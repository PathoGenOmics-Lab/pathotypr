//! Pathotypr Core Library
//!
//! High-performance toolkit for genome classification and variant genotyping.
//!
//! # Modules
//!
//! - `train` — Train Random Forest models on labeled genome sequences
//! - `predict` — Predict lineages using trained models
//! - `classify` — Marker-based variant calling on assemblies
//! - `classify_split_fastq` — Alignment-free genotyping from FASTQ reads
//! - `match` — K-mer containment scoring against reference databases
//!
//! # Supporting modules
//!
//! - `model` — Model bundle, config, and label encoder
//! - `vectorizer` — Feature hashing (hashing trick) for k-mers
//! - `sparse_tree` — Custom CART decision tree on sparse vectors
//! - `fasta_io` — FASTA file parsing for labeled training data
//! - `excel` — Streaming Excel export alongside TSV output
//! - `paired_end` — Paired-end FASTQ detection and grouping
//! - `errors` — Error types and cancellation token
//! - `split_kmer` — Dynamic variant detection engine
//! - `common` — Thread pool config and backward-compat re-exports

// Configuration defaults
pub mod defaults;

// Core workflow modules
pub mod classify;
pub mod classify_split_fastq;
pub mod predict;
pub mod r#match;
pub mod train;

// Supporting modules
pub mod common;
pub mod errors;
pub mod lineage;
pub mod excel;
pub mod fasta_io;
pub mod model;
pub mod paired_end;
pub mod sparse_tree;
pub mod split_kmer;
pub mod vectorizer;

// ---------------------------------------------------------------------------
// Re-exports for convenience
// ---------------------------------------------------------------------------

pub use common::{configure_thread_pool, validate_kmer_size};
pub use errors::{check_cancelled, AppError, AppResult, CancellationToken, ParallelCancellation};
pub use excel::{parse_tsv_for_excel, write_excel_file, ExcelStreamWriter};
pub use fasta_io::{read_fasta, read_fasta_records};
pub use model::{
    LabelEncoder, ModelBundle, ModelConfig, DEFAULT_HASH_BUCKETS, MODEL_FORMAT_VERSION,
};
pub use paired_end::{detect_paired_end_files, PairedEndResult};
pub use vectorizer::{FeatureHasher, KmerCoord};

// Re-export argument types for Tauri integration
pub use classify::Args as ClassifyArgs;
pub use classify_split_fastq::SplitFastqArgs;
pub use predict::PredictArgs;
pub use r#match::Args as MatchArgs;
pub use train::{Args as TrainArgs, TrainReport};
