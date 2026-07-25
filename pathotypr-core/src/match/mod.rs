//! Match subcommand — k-mer containment scoring against reference databases.
//!
//! Uses streaming mode: processes references in batches for constant memory usage.

use crate::common::{configure_thread_pool, detect_paired_end_files};
use crate::errors::{
    check_cancelled, AppError, AppResult, CancellationToken, ParallelCancellation,
};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::sync::{Arc, Mutex};
use std::time::Duration;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const DEFAULT_EARLY_STOP_CONFIDENCE: f64 = 0.0;
const DEFAULT_EARLY_STOP_MIN_KMERS: u64 = 1_000_000;

#[derive(Parser, Debug)]
pub struct Args {
    /// Path to one or more input FASTQ files (can be gzipped). Use this or --input-list.
    #[arg(short = 'i', long = "input", num_args = 1.., group = "input_method")]
    pub fastqs: Vec<String>,

    /// Path to a TSV file listing FASTQ files to process.
    /// Format: sample_name\t/path/to/reads1.fastq[\t/path/to/reads2.fastq]
    /// Use this or --input.
    #[arg(short = 'l', long = "input-list", group = "input_method")]
    pub input_list: Option<String>,

    /// Force paired-end mode: treats input files as paired-end reads, grouped sequentially in pairs.
    /// By default, paired-end files are auto-detected based on naming conventions (_1/_2, _R1/_R2).
    #[arg(long)]
    pub paired: bool,

    /// Disable automatic paired-end detection based on file naming conventions.
    #[arg(long)]
    pub no_auto_paired: bool,

    /// Path to a single multi-FASTA file containing all reference genomes.
    #[arg(short = 'r', long = "references", required = true)]
    pub references: Option<String>,

    /// Path for the output TSV report. Prints to console if not provided.
    #[arg(short = 'o', long)]
    pub output: Option<String>,

    /// k-mer size for comparison.
    #[arg(short = 'k', long, default_value_t = 31)]
    pub kmer_size: u8,

    /// Number of CPU threads to use. Defaults to all available.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,

    /// Stop exact phase early when confidence is high (0 disables early stop).
    /// Recommended range: 0.98 - 0.999.
    #[arg(
        long = "early-stop-confidence",
        default_value_t = DEFAULT_EARLY_STOP_CONFIDENCE
    )]
    pub early_stop_confidence: f64,

    /// Minimum number of scored query k-mers before early-stop checks can trigger.
    #[arg(
        long = "early-stop-min-kmers",
        default_value_t = DEFAULT_EARLY_STOP_MIN_KMERS
    )]
    pub early_stop_min_kmers: u64,

    /// Reproduce legacy Shared_Kmer_Fraction percentages using weighted query k-mer counts.
    #[arg(
        long = "strict-percentages",
        default_value_t = true
    )]
    pub strict_percentages: bool,

    /// Also generate an Excel (.xlsx) file alongside the TSV output.
    #[arg(long, default_value = "false")]
    pub excel: bool,

    /// Minimum k-mer count to keep when filtering query k-mers.
    /// K-mers with fewer occurrences than this are discarded as sequencing noise.
    /// Set to 1 to disable singleton filtering (keep all k-mers).
    /// Default: 2 (remove singletons). Only applied when total unique k-mers > 100k.
    #[arg(long = "min-kmer-count", default_value_t = 2)]
    pub min_kmer_count: u32,

    /// Cancellation token for stopping the task (GUI only, not CLI).
    #[arg(skip)]
    pub cancel_token: Option<CancellationToken>,
}

/// Reads a list of FASTQ file paths from a TSV file.
fn read_fastq_list_from_tsv(path: &str) -> AppResult<Vec<String>> {
    info!("Reading FASTQ file list from TSV: {}", path);
    let mut fastq_paths = Vec::new();
    let content = fs::read_to_string(path)?;
    for line in content.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let paths: Vec<String> = line
            .split('\t')
            .skip(1)
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(String::from)
            .collect();
        if paths.is_empty() {
            warn!("Skipping malformed line in TSV (no paths found): {}", line);
            continue;
        }
        for file_path in &paths {
            if !std::path::Path::new(file_path).exists() {
                return Err(AppError::Io(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    format!(
                        "FASTQ file listed in '{}' does not exist: {}",
                        path, file_path
                    ),
                )));
            }
        }
        fastq_paths.extend(paths);
    }
    if fastq_paths.is_empty() {
        return Err(AppError::NotEnoughData(format!(
            "No FASTQ files found in the TSV list: {}",
            path
        )));
    }
    Ok(fastq_paths)
}

/// Counts k-mer occurrences from a set of FASTQ files in parallel.
/// Filters out low-count k-mers (likely sequencing noise) based on `min_kmer_count`.
fn count_kmers_from_fastqs(
    paths: &[String],
    k: u8,
    min_kmer_count: u32,
    cancel_token: &Option<CancellationToken>,
) -> AppResult<FxHashMap<u64, u32>> {
    info!("🧬 Counting k-mers from input FASTQ files...");

    check_cancelled(cancel_token)?;

    let pb = ProgressBar::new(paths.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("  {spinner:.blue} Reading FASTQs: [{bar:40.cyan/blue}] {pos}/{len} {msg}")?
            .progress_chars("=>-"),
    );

    let cancellation = ParallelCancellation::new(cancel_token);
    let parse_error: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));

    let mut kmer_counts = paths
        .par_iter()
        .map(|path| {
            let parse_error = Arc::clone(&parse_error);

            if cancellation.is_cancelled() {
                return FxHashMap::default();
            }

            let mut local_counts = FxHashMap::default();
            local_counts.reserve(131_072);
            let mut reader = match parse_fastx_file(path) {
                Ok(r) => r,
                Err(e) => {
                    let mut first_error = parse_error.lock().unwrap_or_else(|e| e.into_inner());
                    if first_error.is_none() {
                        *first_error =
                            Some(format!("Failed to open or parse FASTQ '{}': {}", path, e));
                    }
                    pb.inc(1);
                    return local_counts;
                }
            };

            let mut kmer_count = 0u64;
            'outer: while let Some(record) = reader.next() {
                if cancellation.is_cancelled() {
                    break;
                }

                let seqrec = match record {
                    Ok(r) => r,
                    Err(e) => {
                        let mut first_error = parse_error.lock().unwrap_or_else(|e| e.into_inner());
                        if first_error.is_none() {
                            *first_error =
                                Some(format!("Invalid FASTQ record in '{}': {}", path, e));
                        }
                        break;
                    }
                };
                for (_, kmer, _) in seqrec.sequence().bit_kmers(k, true) {
                    *local_counts.entry(kmer.0).or_insert(0) += 1;

                    kmer_count += 1;
                    if kmer_count % 50000 == 0 && cancellation.is_cancelled() {
                        break 'outer;
                    }
                }
            }
            pb.inc(1);
            local_counts
        })
        .reduce(FxHashMap::default, |mut a, b| {
            a.reserve(b.len());
            for (kmer, count) in b {
                *a.entry(kmer).or_insert(0) += count;
            }
            a
        });

    pb.finish_with_message("FASTQ k-mer counting complete.");

    if let Some(message) = parse_error.lock().unwrap_or_else(|e| e.into_inner()).clone() {
        return Err(AppError::Generic(message));
    }
    cancellation.check_after()?;

    if kmer_counts.is_empty() {
        return Ok(kmer_counts);
    }

    // Filter low-count k-mers (noise reduction).
    // Only applied when dataset is large enough that singletons are likely errors.
    if min_kmer_count > 1 {
        let before_count = kmer_counts.len();
        info!("  Raw k-mer count: {}", before_count);

        if before_count > 100_000 {
            let threshold = min_kmer_count;
            kmer_counts.retain(|_, &mut count| count >= threshold);
            let after_count = kmer_counts.len();
            let reduction = 100.0 * (1.0 - (after_count as f64 / before_count as f64));
            info!(
                "📉 Filtered k-mers with count < {} (noise): {} -> {} unique k-mers (reduced by {:.1}%).",
                threshold, before_count, after_count, reduction
            );
        }
    }

    Ok(kmer_counts)
}

/// Extracts all reference sequences and their headers from a multi-FASTA file.
pub fn read_references_from_multifasta(path: &str) -> AppResult<Vec<(String, Vec<u8>)>> {
    info!("📖 Loading reference genomes from multi-FASTA file...");
    let mut refs = Vec::new();
    let mut reader = parse_fastx_file(path)
        .map_err(|_| AppError::Generic(format!("Cannot open reference file: {}", path)))?;

    while let Some(record) = reader.next() {
        let seqrec = record
            .map_err(|_| AppError::Parsing("Invalid record in reference FASTA.".to_string()))?;
        let header = seqrec.id();
        let seq = seqrec.seq().to_vec();
        refs.push((String::from_utf8_lossy(header).to_string(), seq));
    }
    info!(
        "  Found {} reference sequences to compare against.",
        refs.len()
    );
    Ok(refs)
}

/// Adaptive batch size for streaming reference matching.
/// Balances parallelism (need ≥ num_threads refs to keep all cores busy)
/// against memory (each ref in the batch holds ~34 MB of k-mers for a
/// 4.4 Mb genome). Using exactly num_threads gives full parallelism at
/// minimum memory: batch_size × 34 MB for k-mers + batch_size × genome.
fn adaptive_batch_size() -> usize {
    let cpus = rayon::current_num_threads();
    cpus.clamp(4, 64)
}



/// Main execution logic for the `match` subcommand.
pub fn run(args: Args) -> AppResult<()> {
    configure_thread_pool(args.threads);

    let cancel_token = &args.cancel_token;

    check_cancelled(cancel_token)?;

    // Determine fastq files
    let fastq_files_to_process = if let Some(list_path) = &args.input_list {
        read_fastq_list_from_tsv(list_path)?
    } else if !args.fastqs.is_empty() {
        // Check for paired-end detection
        if args.paired {
            // Manual paired mode: validate even number of files
            if args.fastqs.len() % 2 != 0 {
                return Err(AppError::Generic(format!(
                    "--paired requires an even number of input files. Found {}.",
                    args.fastqs.len()
                )));
            }
            info!(
                "🔗 Using manual paired-end mode: {} file pairs",
                args.fastqs.len() / 2
            );
            args.fastqs.clone()
        } else if !args.no_auto_paired && args.fastqs.len() > 1 {
            // Auto-detect paired-end files
            let detection_result = detect_paired_end_files(&args.fastqs);
            if detection_result.is_paired {
                info!(
                    "🔗 Auto-detected paired-end files: {} paired samples, {} single-end samples",
                    detection_result.paired_count, detection_result.single_count
                );
            }
            // Flatten all files back for processing (Match processes all files together)
            detection_result
                .samples
                .values()
                .flatten()
                .cloned()
                .collect()
        } else {
            args.fastqs.clone()
        }
    } else {
        return Err(AppError::Generic(
            "An input source is required: either --input or --input-list must be provided."
                .to_string(),
        ));
    };

    info!(
        "Analyzing FASTQ files: {}",
        fastq_files_to_process.join(", ")
    );

    crate::common::validate_kmer_size(args.kmer_size as usize)?;
    let kmer_size = args.kmer_size;

    check_cancelled(cancel_token)?;

    // Count query k-mers and match against references using streaming mode
    let references_path = args.references.as_ref().ok_or_else(|| {
        AppError::Generic("--references must be provided.".to_string())
    })?;

    let query_kmer_map =
        count_kmers_from_fastqs(&fastq_files_to_process, kmer_size, args.min_kmer_count, cancel_token)?;
    if query_kmer_map.is_empty() {
        return Err(AppError::NotEnoughData(
            "No k-mers could be extracted (or all were filtered as noise).".to_string(),
        ));
    }
    let total_query_kmers = query_kmer_map.values().map(|&c| c as f64).sum::<f64>();
    debug!("Total retained k-mers in query: {}", total_query_kmers);

    info!("🔬 Comparing sample against references...");
    let best_match = match_with_raw_sequences(
        references_path,
        kmer_size,
        &query_kmer_map,
        total_query_kmers,
        cancel_token,
    )?;

    check_cancelled(cancel_token)?;

    info!("✅ Comparison finished. Writing report...");
    let mut writer: Box<dyn Write> = if let Some(path) = &args.output {
        Box::new(File::create(path)?)
    } else {
        Box::new(std::io::stdout())
    };

    writeln!(
        writer,
        "Query_Files\tBest_Match_Reference\tShared_Kmer_Fraction"
    )?;
    if let Some((best_header, best_score)) = best_match {
        let query_files_str = fastq_files_to_process.join(",");
        writeln!(
            writer,
            "{}\t{}\t{:.4}",
            query_files_str, best_header, best_score
        )?;

        if args.excel {
            if let Some(path) = &args.output {
                let headers = [
                    "Query_Files",
                    "Best_Match_Reference",
                    "Shared_Kmer_Fraction",
                ];
                let rows = vec![vec![
                    query_files_str,
                    best_header,
                    format!("{:.4}", best_score),
                ]];
                if let Err(e) = crate::common::write_excel_file(path, &headers, &rows) {
                    warn!("⚠️ Could not write Excel file: {}", e);
                } else {
                    info!("📊 Excel file written to {}", path);
                }
            }
        }
    } else {
        info!("No matching references found.");
    }

    Ok(())
}

/// Match using raw sequences.
/// OPTIMIZED: Streams references in batches to cap memory usage.
fn match_with_raw_sequences(
    references_path: &str,
    kmer_size: u8,
    query_kmer_map: &FxHashMap<u64, u32>,
    total_query_kmers: f64,
    cancel_token: &Option<CancellationToken>,
) -> AppResult<Option<(String, f64)>> {
    let mut reader = parse_fastx_file(references_path).map_err(|_| {
        AppError::Generic(format!("Cannot open reference file: {}", references_path))
    })?;

    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::with_template(
        "  {spinner:.green} Matching references: {pos} processed",
    )?);
    pb.enable_steady_tick(Duration::from_millis(120));

    let cancellation = ParallelCancellation::new(cancel_token);
    let batch_size = adaptive_batch_size();
    info!("  Streaming batch size: {} references per batch.", batch_size);
    let mut reference_batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(batch_size);
    let mut processed_refs = 0usize;
    let mut best_ref_name: Option<String> = None;
    let mut best_shared_kmers = 0u64;

    let process_batch = |batch: &mut Vec<(String, Vec<u8>)>,
                         processed_refs: &mut usize,
                         best_ref_name: &mut Option<String>,
                         best_shared_kmers: &mut u64|
     -> AppResult<()> {
        if batch.is_empty() {
            return Ok(());
        }

        let batch_best = batch
            .par_iter()
            .enumerate()
            .map(|(idx, (_, seq))| {
                if cancellation.is_cancelled() {
                    return (idx, 0u64);
                }

                let mut ref_kmers: Vec<u64> = seq
                    .bit_kmers(kmer_size, true)
                    .map(|(_, kmer, _)| kmer.0)
                    .collect();
                ref_kmers.sort_unstable();
                ref_kmers.dedup();

                let mut shared_kmer_count = 0u64;
                let mut processed = 0usize;
                for kmer in &ref_kmers {
                    processed += 1;
                    if processed % 100000 == 0 && cancellation.is_cancelled() {
                        return (idx, 0u64);
                    }
                    if let Some(count) = query_kmer_map.get(kmer) {
                        shared_kmer_count += *count as u64;
                    }
                }

                (idx, shared_kmer_count)
            })
            .reduce_with(|a, b| {
                if b.1 > a.1 || (b.1 == a.1 && b.0 < a.0) {
                    b
                } else {
                    a
                }
            });

        cancellation.check_after()?;

        if let Some((best_idx, shared_count)) = batch_best {
            if best_ref_name.is_none() || shared_count > *best_shared_kmers {
                *best_ref_name = Some(batch[best_idx].0.clone());
                *best_shared_kmers = shared_count;
            }
        }

        *processed_refs += batch.len();
        pb.set_position(*processed_refs as u64);
        batch.clear();
        Ok(())
    };

    while let Some(record) = reader.next() {
        if processed_refs % batch_size == 0 {
            check_cancelled(cancel_token)?;
        }

        let seqrec = record
            .map_err(|_| AppError::Parsing("Invalid record in reference FASTA.".to_string()))?;
        let header = String::from_utf8_lossy(seqrec.id()).to_string();
        let seq = seqrec.seq().to_vec();
        reference_batch.push((header, seq));
        if reference_batch.len() >= batch_size {
            process_batch(
                &mut reference_batch,
                &mut processed_refs,
                &mut best_ref_name,
                &mut best_shared_kmers,
            )?;
        }
    }
    process_batch(
        &mut reference_batch,
        &mut processed_refs,
        &mut best_ref_name,
        &mut best_shared_kmers,
    )?;

    pb.finish_with_message("Matching complete.");

    if let Some(best_header) = best_ref_name {
        let score = if total_query_kmers > 0.0 {
            best_shared_kmers as f64 / total_query_kmers
        } else {
            0.0
        };
        Ok(Some((best_header, score)))
    } else {
        Ok(None)
    }
}
