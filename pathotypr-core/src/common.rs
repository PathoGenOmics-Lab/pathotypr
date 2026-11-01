//! Common utilities and re-exports.
//!
//! This module re-exports types from focused sub-modules for backward
//! compatibility. New code should import directly from the specific module.

use log::{debug, warn};
// Re-export from sub-modules so existing `use crate::common::X` still works.
pub use crate::excel::{write_excel_file, ExcelStreamWriter};
pub use crate::model::{ModelBundle, MODEL_FORMAT_VERSION};
pub use crate::paired_end::detect_paired_end_files;

// ---------------------------------------------------------------------------
// Outputs written over inputs
// ---------------------------------------------------------------------------

/// Refuses to start when an output would be written over one of the inputs.
///
/// Every command creates its outputs before it has finished reading its
/// inputs, so naming an output after an input destroyed the input: `classify
/// -l run.tsv -o run` truncated the sample list into its own detail table,
/// found no sample in it and deleted it while cleaning up, and `match -r
/// refs.fasta -o refs.fasta`, `train -i data.fasta -o data.fasta` and
/// `classify -m panel.tsv -o panel` replaced the reference, the training data
/// and the marker panel with results and exited successfully. Paths are
/// compared once resolved, so `./run.tsv` or a link to it is the same file.
pub fn refuse_output_over_input<'a>(
    inputs: impl IntoIterator<Item = &'a str>,
    outputs: &[String],
) -> Result<(), crate::errors::AppError> {
    let inputs: Vec<(&str, std::path::PathBuf)> = inputs
        .into_iter()
        .filter_map(|input| std::fs::canonicalize(input).ok().map(|resolved| (input, resolved)))
        .collect();
    for output in outputs {
        // An output that does not exist yet cannot be an existing input.
        let Ok(resolved) = std::fs::canonicalize(output) else {
            continue;
        };
        if let Some((input, _)) = inputs.iter().find(|(_, path)| *path == resolved) {
            return Err(crate::errors::AppError::Generic(format!(
                "The output '{}' is the input '{}'. Writing it would destroy that input; choose another output name.",
                output, input
            )));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// K-mer size validation
// ---------------------------------------------------------------------------

/// Validates that a k-mer size is in the range 1–31 (fits in needletail's
/// 2-bit encoding with u64). Call this at CLI/GUI entry points.
pub fn validate_kmer_size(k: usize) -> Result<usize, crate::errors::AppError> {
    if k == 0 || k > 31 {
        Err(crate::errors::AppError::Generic(format!(
            "k-mer size must be between 1 and 31, got {}",
            k
        )))
    } else {
        Ok(k)
    }
}

// ---------------------------------------------------------------------------
// Thread pool
// ---------------------------------------------------------------------------

/// Configures the global rayon thread pool.
///
/// Handles the case where the pool has already been initialized (common in
/// long-running GUI processes) by logging a warning instead of panicking.
pub fn configure_thread_pool(num_threads: Option<usize>) {
    if let Some(n) = num_threads {
        debug!("Setting number of threads to: {}", n);
        match rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
        {
            Ok(()) => debug!("Global thread pool initialized with {} threads.", n),
            Err(e) => {
                warn!(
                    "Could not configure thread pool ({}). Using existing pool configuration.",
                    e
                );
            }
        }
    }
}

#[cfg(test)]
mod output_over_input_tests {
    use super::refuse_output_over_input;

    /// The case that destroyed a sample list: the output resolves to the input,
    /// even written another way.
    #[test]
    fn an_output_that_is_an_input_is_refused() {
        let dir = std::env::temp_dir().join(format!("pathotypr-ooi-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("run.tsv");
        std::fs::write(&input, "s1\tgenome.fasta\n").unwrap();
        let input = input.to_str().unwrap().to_string();
        let same_written_otherwise = format!("{}/./run.tsv", dir.to_str().unwrap());
        assert!(refuse_output_over_input([input.as_str()], &[input.clone()]).is_err());
        assert!(refuse_output_over_input([input.as_str()], &[same_written_otherwise]).is_err());
        // A different file, and one that does not exist yet, are both fine.
        let other = dir.join("run_summary.tsv").to_str().unwrap().to_string();
        assert!(refuse_output_over_input([input.as_str()], &[other]).is_ok());
        std::fs::remove_dir_all(&dir).ok();
    }
}
