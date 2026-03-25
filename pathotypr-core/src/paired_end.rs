//! Paired-end FASTQ file detection and grouping.
//!
//! Detects and groups paired-end FASTQ files based on common naming conventions
//! (Illumina `_R1_001`/`_R2_001`, `_R1`/`_R2`, `_1`/`_2`, `.1`/`.2`).

use log::debug;
use regex::Regex;
use std::collections::HashMap;
use std::path::Path;
use std::sync::OnceLock;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Result of paired-end detection for a set of FASTQ files.
#[derive(Debug, Clone)]
pub struct PairedEndResult {
    /// Sample name → FASTQ file paths (2 for paired, 1 for single-end).
    pub samples: HashMap<String, Vec<String>>,
    /// Whether any paired-end files were detected.
    pub is_paired: bool,
    /// Number of paired samples detected.
    pub paired_count: usize,
    /// Number of single-end samples detected.
    pub single_count: usize,
}

// ---------------------------------------------------------------------------
// Detection
// ---------------------------------------------------------------------------

/// Detects and groups paired-end FASTQ files based on naming conventions.
pub fn detect_paired_end_files(paths: &[String]) -> PairedEndResult {
    static PAIRED_PATTERNS: OnceLock<Vec<(Regex, Regex, &'static str)>> = OnceLock::new();
    let patterns = PAIRED_PATTERNS.get_or_init(|| {
        vec![
            (
                Regex::new(r"(.+?)_R1_\d+\.").unwrap(),
                Regex::new(r"(.+?)_R2_\d+\.").unwrap(),
                "_R1_",
            ),
            (
                Regex::new(r"(.+?)_R1[._]").unwrap(),
                Regex::new(r"(.+?)_R2[._]").unwrap(),
                "_R1",
            ),
            (
                Regex::new(r"(.+?)_1\.").unwrap(),
                Regex::new(r"(.+?)_2\.").unwrap(),
                "_1",
            ),
            (
                Regex::new(r"(.+?)\.1\.").unwrap(),
                Regex::new(r"(.+?)\.2\.").unwrap(),
                ".1",
            ),
        ]
    });

    fn extract_base_name(
        path: &str,
        patterns: &[(Regex, Regex, &'static str)],
    ) -> Option<(String, u8)> {
        let filename = Path::new(path)
            .file_name()
            .and_then(|f| f.to_str())
            .unwrap_or("");

        for (r1_pattern, r2_pattern, _) in patterns {
            if let Some(caps) = r1_pattern.captures(filename) {
                if let Some(base) = caps.get(1) {
                    return Some((base.as_str().to_string(), 1));
                }
            }
            if let Some(caps) = r2_pattern.captures(filename) {
                if let Some(base) = caps.get(1) {
                    return Some((base.as_str().to_string(), 2));
                }
            }
        }
        None
    }

    // Group files by their base name.
    let mut grouped: HashMap<String, (Option<String>, Option<String>)> = HashMap::new();
    let mut unmatched: Vec<String> = Vec::new();

    for path in paths {
        if let Some((base_name, read_num)) = extract_base_name(path, patterns) {
            let entry = grouped.entry(base_name).or_insert((None, None));
            match read_num {
                1 => entry.0 = Some(path.clone()),
                2 => entry.1 = Some(path.clone()),
                _ => unmatched.push(path.clone()),
            }
        } else {
            unmatched.push(path.clone());
        }
    }

    // Build the result.
    let mut samples: HashMap<String, Vec<String>> = HashMap::new();
    let mut paired_count = 0;
    let mut single_count = 0;

    for (base_name, (r1, r2)) in grouped {
        match (r1, r2) {
            (Some(r1_path), Some(r2_path)) => {
                samples.insert(base_name, vec![r1_path, r2_path]);
                paired_count += 1;
            }
            (Some(path), None) | (None, Some(path)) => {
                let sample_name = derive_sample_name(&path);
                samples.insert(sample_name, vec![path]);
                single_count += 1;
            }
            (None, None) => {}
        }
    }

    for path in unmatched {
        let mut sample_name = derive_sample_name(&path);
        if samples.contains_key(&sample_name) {
            sample_name = Path::new(&path)
                .file_stem()
                .unwrap_or_default()
                .to_str()
                .unwrap_or(&path)
                .to_string();
        }
        samples.entry(sample_name).or_default().push(path);
        single_count += 1;
    }

    debug!(
        "Paired-end detection: {} paired, {} single-end",
        paired_count, single_count
    );

    PairedEndResult {
        samples,
        is_paired: paired_count > 0,
        paired_count,
        single_count,
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Derives a clean sample name from a file path by stripping extensions and
/// paired-end suffixes (e.g., `_R1_001`, `-R2`, `_1`).
pub(crate) fn derive_sample_name(path: &str) -> String {
    let filename = Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    let mut base = filename.to_string();

    // Handle compressed FASTQ names (e.g., sample_R1.fastq.gz).
    let lower = base.to_ascii_lowercase();
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if lower.ends_with(ext) {
            base.truncate(base.len().saturating_sub(ext.len()));
            break;
        }
    }

    let suffixes = [
        "_R1_001", "_R2_001", "_r1_001", "_r2_001", "_R1", "_R2", "_r1", "_r2", "_1", "_2", ".1",
        ".2", "-R1", "-R2", "-r1", "-r2", "-1", "-2",
    ];
    for suffix in suffixes {
        if let Some(stripped) = base.strip_suffix(suffix) {
            if !stripped.is_empty() {
                base = stripped.to_string();
                break;
            }
        }
    }

    let normalized = base
        .trim_end_matches(|c| c == '_' || c == '-' || c == '.')
        .to_string();
    if normalized.is_empty() {
        Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string()
    } else {
        normalized
    }
}
