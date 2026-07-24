//! Paired-end FASTQ file detection and grouping.
//!
//! Detects and groups paired-end FASTQ files based on common naming conventions
//! (Illumina `_R1_001`/`_R2_001`, `_R1`/`_R2`, `_1`/`_2`, `.1`/`.2`).

use log::{debug, warn};
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

    // Group files by their directory *and* base name. Accumulating every R1 and
    // R2 path keeps chunked / lane-split output for the same sample (e.g.
    // `_R1_001`, `_R1_002`) together, while the directory component keeps
    // identically named files from different run folders apart — merging those
    // would silently build a chimera of two distinct isolates.
    let mut grouped: HashMap<(String, String), (Vec<String>, Vec<String>)> = HashMap::new();
    let mut unmatched: Vec<String> = Vec::new();

    for path in paths {
        if let Some((base_name, read_num)) = extract_base_name(path, patterns) {
            let dir = Path::new(path)
                .parent()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_default();
            let entry = grouped.entry((dir, base_name)).or_default();
            match read_num {
                1 => entry.0.push(path.clone()),
                2 => entry.1.push(path.clone()),
                _ => unmatched.push(path.clone()),
            }
        } else {
            unmatched.push(path.clone());
        }
    }

    // Build the result. Sort the groups first: `insert_unique_sample` resolves
    // name collisions by suffixing, so a HashMap's arbitrary iteration order
    // would make the assignment of those suffixes differ between identical
    // runs. Sorting keeps sample naming reproducible.
    let mut grouped: Vec<((String, String), (Vec<String>, Vec<String>))> =
        grouped.into_iter().collect();
    grouped.sort_by(|((d1, b1), _), ((d2, b2), _)| (b1, d1).cmp(&(b2, d2)));

    let mut samples: HashMap<String, Vec<String>> = HashMap::new();
    let mut paired_count = 0;
    let mut single_count = 0;

    for ((_dir, base_name), (mut r1, mut r2)) in grouped {
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
        warn!(
            "Sample name '{}' is already in use; {:?} will be reported as '{}'. \
             Rename the files or use --input-list to control sample names.",
            base, files, name
        );
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

#[cfg(test)]
mod tests {
    use super::*;

    fn v(paths: &[&str]) -> Vec<String> {
        paths.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn groups_r1_r2_into_one_paired_sample() {
        let r = detect_paired_end_files(&v(&["/d/S1_R1.fastq.gz", "/d/S1_R2.fastq.gz"]));
        assert_eq!(r.paired_count, 1);
        assert_eq!(r.single_count, 0);
        assert!(r.is_paired);
        assert_eq!(r.samples["S1"].len(), 2);
    }

    #[test]
    fn keeps_every_chunk_of_a_split_sample() {
        // Lane/chunk-split output must not lose reads: all four files belong
        // to the same sample.
        let r = detect_paired_end_files(&v(&[
            "/d/S_L001_R1_001.fastq.gz",
            "/d/S_L001_R1_002.fastq.gz",
            "/d/S_L001_R2_001.fastq.gz",
            "/d/S_L001_R2_002.fastq.gz",
        ]));
        assert_eq!(r.paired_count, 1);
        let files = r.samples.values().next().unwrap();
        assert_eq!(files.len(), 4, "no chunk may be dropped");
    }

    #[test]
    fn same_basename_in_different_dirs_stays_separate() {
        // Two distinct isolates that happen to share a filename must never be
        // merged into one chimeric sample.
        let r = detect_paired_end_files(&v(&[
            "/runA/S1_R1.fastq.gz",
            "/runA/S1_R2.fastq.gz",
            "/runB/S1_R1.fastq.gz",
            "/runB/S1_R2.fastq.gz",
        ]));
        assert_eq!(r.paired_count, 2, "each directory is its own sample");
        assert_eq!(r.samples.len(), 2);
        for files in r.samples.values() {
            assert_eq!(files.len(), 2);
            let dirs: std::collections::HashSet<_> = files
                .iter()
                .map(|f| Path::new(f).parent().unwrap().to_path_buf())
                .collect();
            assert_eq!(dirs.len(), 1, "a sample must not mix directories");
        }
    }

    #[test]
    fn unmatched_file_does_not_contaminate_an_existing_sample() {
        let r = detect_paired_end_files(&v(&["/d/foo_1.fastq", "/d/foo_2.fastq", "/d/foo.fastq"]));
        assert_eq!(r.samples.len(), 2, "the lone file gets its own sample");
        let paired = r.samples.values().find(|f| f.len() == 2).unwrap();
        assert!(paired.iter().all(|f| f.ends_with("_1.fastq") || f.ends_with("_2.fastq")));
    }

    #[test]
    fn sample_naming_is_deterministic() {
        let paths = v(&[
            "/runA/S1_R1.fastq.gz",
            "/runA/S1_R2.fastq.gz",
            "/runB/S1_R1.fastq.gz",
            "/runB/S1_R2.fastq.gz",
        ]);
        let first: Vec<String> = {
            let mut k: Vec<String> = detect_paired_end_files(&paths).samples.keys().cloned().collect();
            k.sort();
            k
        };
        for _ in 0..8 {
            let mut k: Vec<String> = detect_paired_end_files(&paths).samples.keys().cloned().collect();
            k.sort();
            assert_eq!(k, first, "sample names must not vary between runs");
        }
    }
}
