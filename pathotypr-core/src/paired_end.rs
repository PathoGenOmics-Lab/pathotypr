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

    // Group files by their base name. Accumulate every R1 and R2 path so that
    // chunked / lane-split output for the same sample (e.g. `_R1_001`,
    // `_R1_002`) is kept together instead of the later file overwriting the
    // earlier one.
    let mut grouped: HashMap<String, (Vec<String>, Vec<String>)> = HashMap::new();
    let mut unmatched: Vec<String> = Vec::new();

    for path in paths {
        if let Some((base_name, read_num)) = extract_base_name(path, patterns) {
            let entry = grouped.entry(base_name).or_default();
            match read_num {
                1 => entry.0.push(path.clone()),
                2 => entry.1.push(path.clone()),
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

    for (base_name, (mut r1, mut r2)) in grouped {
        if !r1.is_empty() && !r2.is_empty() {
            // Paired: all R1 files followed by all R2 files. scan_fastq_with_index
            // scans an arbitrary number of input files for the sample.
            let mut files = Vec::with_capacity(r1.len() + r2.len());
            files.append(&mut r1);
            files.append(&mut r2);
            insert_unique_sample(&mut samples, base_name, files);
            paired_count += 1;
        } else if !r1.is_empty() || !r2.is_empty() {
            // Only one read direction present: treat as single-end.
            let mut files = r1;
            files.append(&mut r2);
            let sample_name = derive_sample_name(&files[0]);
            insert_unique_sample(&mut samples, sample_name, files);
            single_count += 1;
        }
    }

    for path in unmatched {
        let sample_name = derive_sample_name(&path);
        insert_unique_sample(&mut samples, sample_name, vec![path]);
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

/// Inserts a sample under a unique key, appending an incrementing numeric
/// suffix on collision. This guarantees that two distinct inputs never merge
/// into the same sample bucket (which would silently contaminate one sample's
/// read counts with another's).
fn insert_unique_sample(
    samples: &mut HashMap<String, Vec<String>>,
    desired_name: String,
    files: Vec<String>,
) {
    let mut name = desired_name;
    if samples.contains_key(&name) {
        let base = name.clone();
        let mut n = 2;
        loop {
            name = format!("{}_{}", base, n);
            if !samples.contains_key(&name) {
                break;
            }
            n += 1;
        }
    }
    samples.insert(name, files);
}

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
