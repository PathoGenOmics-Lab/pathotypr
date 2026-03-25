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
