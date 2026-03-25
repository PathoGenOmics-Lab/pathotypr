//! Utility functions: file path validation, SSRF protection, and system usage.

use serde::Serialize;
use std::net::{IpAddr, Ipv4Addr, Ipv6Addr, ToSocketAddrs};
use std::path::PathBuf;
use std::process::Command;
use tauri::{AppHandle, Manager};

// ---------------------------------------------------------------------------
// File path validation
// ---------------------------------------------------------------------------

/// Resolve and validate a file path against Tauri's allowed directories.
pub fn resolve_allowed_file_path(app: &AppHandle, path: &str) -> Result<PathBuf, String> {
    let requested = PathBuf::from(path.trim());
    if requested.as_os_str().is_empty() {
        return Err("Path is empty".to_string());
    }
    if !requested.is_absolute() {
        return Err("Path must be absolute".to_string());
    }

    let canonical = requested
        .canonicalize()
        .map_err(|e| format!("Failed to resolve path '{}': {}", requested.display(), e))?;
    if !canonical.is_file() {
        return Err(format!("Path is not a file: {}", canonical.display()));
    }

    let allowed_roots = {
        let path_api = app.path();
        [
            path_api.document_dir().ok(),
            path_api.download_dir().ok(),
            path_api.desktop_dir().ok(),
            path_api.app_data_dir().ok(),
            path_api.temp_dir().ok(),
        ]
        .into_iter()
        .flatten()
        .map(|p: PathBuf| p.canonicalize().unwrap_or(p))
        .collect::<Vec<PathBuf>>()
    };

    let is_allowed = allowed_roots
        .iter()
        .any(|root: &PathBuf| canonical.starts_with(root.as_path()));
    if !is_allowed {
        return Err(format!(
            "Reading files outside user/app directories is not allowed: {}",
            canonical.display()
        ));
    }

    Ok(canonical)
}

// ---------------------------------------------------------------------------
// SSRF protection
// ---------------------------------------------------------------------------

fn is_disallowed_ipv4(ip: Ipv4Addr) -> bool {
    let [a, b, c, _] = ip.octets();
    ip.is_private()
        || ip.is_loopback()
        || ip.is_link_local()
        || ip.is_broadcast()
        || ip.is_unspecified()
        // Carrier-grade NAT
        || (a == 100 && (64..=127).contains(&b))
        // Documentation/testing and benchmark ranges
        || (a == 192 && b == 0 && (c == 0 || c == 2))
        || (a == 198 && (b == 18 || b == 19))
        || (a == 198 && b == 51 && c == 100)
        || (a == 203 && b == 0 && c == 113)
        // Multicast/reserved
        || a >= 224
}

fn is_disallowed_ipv6(ip: Ipv6Addr) -> bool {
    ip.is_loopback()
        || ip.is_unspecified()
        || ip.is_unique_local()
        || ip.is_unicast_link_local()
        || ip.is_multicast()
}

fn is_disallowed_ip(ip: IpAddr) -> bool {
    match ip {
        IpAddr::V4(v4) => is_disallowed_ipv4(v4),
        IpAddr::V6(v6) => is_disallowed_ipv6(v6),
    }
}

fn validate_host_resolution(host: &str, port: u16) -> Result<(), String> {
    let addrs = (host, port)
        .to_socket_addrs()
        .map_err(|e| format!("Could not resolve host '{}': {}", host, e))?;

    let mut resolved_any = false;
    for addr in addrs {
        resolved_any = true;
        if is_disallowed_ip(addr.ip()) {
            return Err("Cannot download from local/reserved addresses".to_string());
        }
    }

    if !resolved_any {
        return Err(format!(
            "Could not resolve host '{}' to any reachable address",
            host
        ));
    }

    Ok(())
}

pub fn validate_download_url(url: &reqwest::Url) -> Result<(), String> {
    if !matches!(url.scheme(), "https") {
        return Err("Only HTTPS URLs are allowed".to_string());
    }

    let host = url
        .host_str()
        .ok_or_else(|| "URL is missing a host".to_string())?
        .trim()
        .to_ascii_lowercase();
    if host == "localhost" || host.ends_with(".localhost") || host.ends_with(".local") {
        return Err("Cannot download from local/reserved addresses".to_string());
    }

    if let Ok(ip) = host.parse::<IpAddr>() {
        if is_disallowed_ip(ip) {
            return Err("Cannot download from local/reserved addresses".to_string());
        }
    }

    let port = url
        .port_or_known_default()
        .ok_or_else(|| "URL has an unknown port".to_string())?;
    validate_host_resolution(&host, port)
}

// ---------------------------------------------------------------------------
// System usage monitoring
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
pub struct SystemUsage {
    pub cpu_percent: f32,
    pub memory_percent: f32,
    pub memory_mb: f32,
    pub memory_rss_bytes: u64,
}

fn parse_decimal_value(token: &str) -> Option<f64> {
    let normalized = token
        .trim()
        .trim_end_matches('%')
        .trim_matches(|c: char| c == ',' || c == ';')
        .replace(',', ".");
    normalized.parse::<f64>().ok()
}

#[cfg(any(target_os = "macos", target_os = "linux"))]
fn parse_ps_usage_output(output: &str) -> Option<SystemUsage> {
    let line = output.lines().find(|line| !line.trim().is_empty())?;
    let mut parts = line.split_whitespace();
    let cpu_percent = parse_decimal_value(parts.next()?)?.max(0.0);
    let memory_percent = parse_decimal_value(parts.next()?)?.max(0.0);
    let rss_kib = parse_decimal_value(parts.next()?)?.max(0.0);
    let memory_rss_bytes = (rss_kib * 1024.0).round().max(0.0) as u64;
    let memory_mb = memory_rss_bytes as f64 / (1024.0 * 1024.0);

    Some(SystemUsage {
        cpu_percent: cpu_percent as f32,
        memory_percent: memory_percent as f32,
        memory_mb: memory_mb as f32,
        memory_rss_bytes,
    })
}

#[cfg(target_os = "windows")]
fn parse_wmic_value(text: &str, key: &str) -> Option<f64> {
    text.lines().find_map(|line| {
        let trimmed = line.trim();
        if !trimmed
            .to_ascii_lowercase()
            .starts_with(&(key.to_ascii_lowercase() + "="))
        {
            return None;
        }
        let value = trimmed.split_once('=')?.1.trim();
        value.parse::<f64>().ok()
    })
}

#[cfg(any(target_os = "macos", target_os = "linux"))]
pub fn read_system_usage() -> Result<SystemUsage, String> {
    let pid = std::process::id().to_string();
    let output = Command::new("ps")
        .arg("-p")
        .arg(&pid)
        .arg("-o")
        .arg("%cpu=")
        .arg("-o")
        .arg("%mem=")
        .arg("-o")
        .arg("rss=")
        .output()
        .map_err(|e| format!("Failed to read app process usage: {}", e))?;

    if !output.status.success() {
        return Err("Failed to read app process usage from ps".to_string());
    }

    let text = String::from_utf8_lossy(&output.stdout);
    parse_ps_usage_output(&text).ok_or_else(|| "Could not parse app process usage".to_string())
}

#[cfg(target_os = "windows")]
pub fn read_system_usage() -> Result<SystemUsage, String> {
    let pid = std::process::id().to_string();
    let filter = format!("IDProcess={}", pid);

    let perf_output = Command::new("wmic")
        .arg("path")
        .arg("Win32_PerfFormattedData_PerfProc_Process")
        .arg("where")
        .arg(&filter)
        .arg("get")
        .arg("PercentProcessorTime,WorkingSetPrivate")
        .arg("/value")
        .output()
        .map_err(|e| format!("Failed to read app process usage: {}", e))?;
    if !perf_output.status.success() {
        return Err("Failed to read app process usage".to_string());
    }
    let perf_text = String::from_utf8_lossy(&perf_output.stdout);
    let cpu_percent = parse_wmic_value(&perf_text, "PercentProcessorTime")
        .unwrap_or(0.0)
        .max(0.0);
    let memory_rss_bytes = parse_wmic_value(&perf_text, "WorkingSetPrivate")
        .or_else(|| parse_wmic_value(&perf_text, "WorkingSet"))
        .unwrap_or(0.0)
        .max(0.0) as u64;

    let mem_output = Command::new("wmic")
        .args(["OS", "get", "TotalVisibleMemorySize", "/value"])
        .output()
        .map_err(|e| format!("Failed to read memory usage: {}", e))?;
    if !mem_output.status.success() {
        return Err("Failed to read memory usage".to_string());
    }
    let mem_text = String::from_utf8_lossy(&mem_output.stdout);
    let total_kib = parse_wmic_value(&mem_text, "TotalVisibleMemorySize").unwrap_or(0.0);
    let total_bytes = (total_kib * 1024.0).max(0.0);
    let memory_percent = if total_bytes > 0.0 {
        (memory_rss_bytes as f64 / total_bytes * 100.0).max(0.0)
    } else {
        0.0
    };
    let memory_mb = memory_rss_bytes as f64 / (1024.0 * 1024.0);

    Ok(SystemUsage {
        cpu_percent: cpu_percent as f32,
        memory_percent: memory_percent as f32,
        memory_mb: memory_mb as f32,
        memory_rss_bytes,
    })
}

#[cfg(not(any(target_os = "macos", target_os = "linux", target_os = "windows")))]
pub fn read_system_usage() -> Result<SystemUsage, String> {
    Err("System usage monitoring is not supported on this platform".to_string())
}
