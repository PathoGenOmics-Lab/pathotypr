//! Implements the `SplitFastq` subcommand.
//!
//! This module now supports batch processing of multiple samples, either from
//! direct command-line input or from a TSV file. It automatically determines
//! sample names and generates separate reports for each, plus a final summary.
//! It uses a dynamic k-mer engine to detect SNPs, MNVs, and Indels.
//! MODIFIED: Includes strictly corrected logic for nested lineage classification.

use crate::common::{configure_thread_pool, detect_paired_end_files};
use crate::errors::{check_cancelled, AppError, AppResult, CancellationToken};
use crate::lineage::determine_major_lineage;
use crate::paired_end::derive_sample_name;
use crate::split_kmer;
use clap::Parser;
use log::{debug, info, trace, warn};
use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::UNIX_EPOCH;

#[derive(Parser, Debug)]
pub struct SplitFastqArgs {
    /// One or more FASTQ files to analyze (e.g., -i sample_R1.fq -i sample_R2.fq).
    /// Use this or --input-list.
    #[arg(short = 'i', long = "input", group = "input_method")]
    pub input: Vec<String>,

    /// Path to a TSV file listing samples.
    /// Format: sample_name\tpath/to/R1.fq[\tpath/to/R2.fq].
    /// Use this or --input.
    #[arg(short = 'l', long = "input-list", group = "input_method")]
    pub input_list: Option<String>,

    /// Force paired-end mode: treats input files as paired-end reads, grouped sequentially in pairs.
    /// By default, paired-end files are auto-detected based on naming conventions (_1/_2, _R1/_R2).
    #[arg(long)]
    pub paired: bool,

    /// Disable automatic paired-end detection based on file naming conventions.
    /// Use this if you want to treat all files as single-end samples.
    #[arg(long)]
    pub no_auto_paired: bool,

    /// Reference FASTA file used to define the markers.
    #[arg(short = 'r', long, required = true)]
    pub reference: String,

    /// Path to a TSV file defining markers.
    /// Format: position\tREF\tALT\tlevel1\tlevel2...
    #[arg(short = 'm', long, required = true)]
    pub markers: String,

    /// Number of threads for parallel processing. Defaults to all available cores.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,

    /// Prefix for the output files.
    #[arg(short = 'o', long, default_value = "split")]
    pub output_prefix: String,

    /// Minimum read depth required to call a variant at a marker position.
    #[arg(long, default_value_t = 10)]
    pub min_depth: u32,

    /// Minimum frequency of the alternate allele to call a variant, as a percentage.
    #[arg(long, default_value_t = 95)]
    pub min_alt_percent: u32,

    /// [NEW] Enable nested lineage classification logic.
    /// This requires markers to be defined with multiple columns for each level.
    #[arg(long)]
    pub nested_classification: bool,

    /// Size of diagnostic k-mers (must be odd, between 11 and 31). Default: 31.
    #[arg(short = 'k', long, default_value_t = split_kmer::DEFAULT_MARKER_KMER_LEN)]
    pub kmer_size: usize,

    /// Also generate Excel (.xlsx) files alongside TSV outputs.
    #[arg(long, default_value = "false")]
    pub excel: bool,

    /// Cancellation token for stopping the task (GUI only, not CLI).
    #[arg(skip)]
    pub cancel_token: Option<CancellationToken>,
}

/// Sanitize a sample name for safe file-name usage.
fn sanitize_sample_name_for_path(sample_name: &str) -> String {
    let sanitized: String = sample_name
        .chars()
        .map(|c| match c {
            'A'..='Z' | 'a'..='z' | '0'..='9' | '_' | '-' | '.' => c,
            _ => '_',
        })
        .collect();
    // Trim leading dots and surrounding underscores in a single pass (O(n))
    let sanitized = sanitized
        .trim_start_matches('.')
        .trim_matches('_');
    if sanitized.is_empty() {
        "sample".to_string()
    } else {
        sanitized.to_string()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MarkerCacheFingerprint {
    reference_canonical_path: String,
    reference_len: u64,
    reference_modified_unix_nanos: u128,
    marker_canonical_path: String,
    marker_len: u64,
    marker_modified_unix_nanos: u128,
    kmer_len: usize,
}

#[derive(Clone)]
struct CachedMarkerDatabase {
    fingerprint: MarkerCacheFingerprint,
    markers: Arc<Vec<split_kmer::Marker>>,
    marker_index: Arc<split_kmer::MarkerIndex>,
    bloom_filter: Arc<split_kmer::BloomFilter>,
}

static SPLITFQ_MARKER_CACHE: OnceLock<Mutex<Option<CachedMarkerDatabase>>> = OnceLock::new();

fn splitfq_marker_cache() -> &'static Mutex<Option<CachedMarkerDatabase>> {
    SPLITFQ_MARKER_CACHE.get_or_init(|| Mutex::new(None))
}

fn file_fingerprint(path: &str) -> AppResult<(String, u64, u128)> {
    let canonical_path = fs::canonicalize(path).unwrap_or_else(|_| PathBuf::from(path));
    let metadata = fs::metadata(path)?;
    let modified_unix_nanos = metadata
        .modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    Ok((
        canonical_path.to_string_lossy().to_string(),
        metadata.len(),
        modified_unix_nanos,
    ))
}

fn build_marker_cache_fingerprint(
    reference_path: &str,
    marker_tsv_path: &str,
    kmer_len: usize,
) -> AppResult<MarkerCacheFingerprint> {
    let (reference_canonical_path, reference_len, reference_modified_unix_nanos) =
        file_fingerprint(reference_path)?;
    let (marker_canonical_path, marker_len, marker_modified_unix_nanos) =
        file_fingerprint(marker_tsv_path)?;
    Ok(MarkerCacheFingerprint {
        reference_canonical_path,
        reference_len,
        reference_modified_unix_nanos,
        marker_canonical_path,
        marker_len,
        marker_modified_unix_nanos,
        kmer_len,
    })
}

fn load_or_build_marker_database(
    reference_path: &str,
    marker_tsv_path: &str,
    cancel_token: &Option<CancellationToken>,
    kmer_len: usize,
) -> AppResult<(Arc<Vec<split_kmer::Marker>>, Arc<split_kmer::MarkerIndex>, Arc<split_kmer::BloomFilter>)> {
    let fingerprint = build_marker_cache_fingerprint(reference_path, marker_tsv_path, kmer_len)?;

    if let Some(entry) = splitfq_marker_cache()
        .lock()
        .map_err(|e| AppError::Generic(format!("marker cache lock poisoned: {e}")))?
        .as_ref()
        .filter(|entry| entry.fingerprint == fingerprint)
    {
        info!("♻️ Reusing cached marker database and marker index from memory.");
        return Ok((entry.markers.clone(), entry.marker_index.clone(), entry.bloom_filter.clone()));
    }

    check_cancelled(cancel_token)?;

    info!("▶ Building dynamic marker database...");
    let markers = split_kmer::build_markers(reference_path, marker_tsv_path, cancel_token, kmer_len)
        .map_err(|e| {
            if cancel_token
                .as_ref()
                .map(|t| t.is_cancelled())
                .unwrap_or(false)
                || e.to_string().to_lowercase().contains("cancel")
            {
                AppError::Cancelled
            } else {
                AppError::Generic(format!("Failed to build markers: {}", e))
            }
        })?;
    info!(
        "  Successfully generated {} dynamic markers.",
        markers.len()
    );
    trace!("First marker generated: {:?}", markers.first());

    check_cancelled(cancel_token)?;

    info!(
        "▶ Building reusable marker index for {} markers...",
        markers.len()
    );
    let (marker_index, bloom_filter) = split_kmer::build_marker_index(&markers);

    let markers_arc = Arc::new(markers);
    let marker_index_arc = Arc::new(marker_index);
    let bloom_filter_arc = Arc::new(bloom_filter);

    *splitfq_marker_cache().lock().unwrap_or_else(|e| e.into_inner()) = Some(CachedMarkerDatabase {
        fingerprint,
        markers: markers_arc.clone(),
        marker_index: marker_index_arc.clone(),
        bloom_filter: bloom_filter_arc.clone(),
    });
    info!("🧠 Cached marker database/index in memory.");

    Ok((markers_arc, marker_index_arc, bloom_filter_arc))
}

/// Reads a TSV file to get a map of sample names to their FASTQ file paths.
fn read_sample_list(path: &str) -> AppResult<HashMap<String, Vec<String>>> {
    let mut samples = HashMap::new();
    let reader = fs::read_to_string(path)?;
    for line in reader.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 2 {
            warn!("Skipping malformed line in sample list: {}", line);
            continue;
        }
        let sample_name = fields[0].trim().to_string();
        if sample_name.is_empty() {
            warn!("Skipping sample with empty name in sample list: {}", line);
            continue;
        }
        let fastq_paths: Vec<String> = fields[1..]
            .iter()
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect();
        if fastq_paths.is_empty() {
            warn!("Skipping sample {} with no FASTQ paths", sample_name);
            continue;
        }
        if samples.insert(sample_name.clone(), fastq_paths).is_some() {
            return Err(AppError::Generic(format!(
                "Duplicate sample name '{}' in sample list '{}'. Sample names must be unique.",
                sample_name, path
            )));
        }
        for fq_path in samples.get(&sample_name).into_iter().flatten() {
            if !Path::new(fq_path).exists() {
                return Err(AppError::Io(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    format!(
                        "FASTQ file '{}' for sample '{}' does not exist",
                        fq_path, sample_name
                    ),
                )));
            }
        }
    }
    if samples.is_empty() {
        return Err(AppError::NotEnoughData(format!(
            "No valid samples found in sample list '{}'.",
            path
        )));
    }
    debug!("Read {} samples from list {}", samples.len(), path);
    Ok(samples)
}

fn cleanup_generated_outputs(paths: &[String]) {
    for out_path in paths {
        match fs::remove_file(out_path) {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => warn!("Failed to remove partial output {}: {}", out_path, e),
        }

        // Must match the writer's derivation exactly, or cleanup deletes the
        // wrong file and leaves the real partial .xlsx behind.
        let xlsx_path = crate::excel::excel_path_from_tsv(out_path);
        match fs::remove_file(&xlsx_path) {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => warn!(
                "Failed to remove partial output {}: {}",
                xlsx_path,
                e
            ),
        }
    }
}

fn check_cancelled_with_cleanup(
    cancel_token: &Option<CancellationToken>,
    generated_outputs: &[String],
) -> AppResult<()> {
    match check_cancelled(cancel_token) {
        Ok(()) => Ok(()),
        Err(e) => {
            cleanup_generated_outputs(generated_outputs);
            Err(e)
        }
    }
}

pub fn run(args: SplitFastqArgs) -> AppResult<()> {
    if args.input.is_empty() && args.input_list.is_none() {
        return Err(AppError::Generic(
            "You must provide an input source: either --input or --input-list.".to_string(),
        ));
    }

    check_cancelled(&args.cancel_token)?;

    configure_thread_pool(args.threads);

    crate::common::validate_kmer_size(args.kmer_size)?;
    let kmer_len = args.kmer_size;
    info!("Using k-mer size: {}", kmer_len);

    let (markers, marker_index, bloom_filter_arc) =
        load_or_build_marker_database(&args.reference, &args.markers, &args.cancel_token, kmer_len)?;

    let marker_lineage_paths: Vec<String> = markers
        .iter()
        .map(|marker| marker.lineages.join(";"))
        .collect();
    let marker_lineage_prefixes: Vec<Vec<String>> = markers
        .iter()
        .map(|marker| {
            let mut prefixes = Vec::with_capacity(marker.lineages.len());
            let mut current = String::new();
            for (idx, component) in marker.lineages.iter().enumerate() {
                if idx == 0 {
                    current.push_str(component);
                } else {
                    current.push(';');
                    current.push_str(component);
                }
                prefixes.push(current.clone());
            }
            prefixes
        })
        .collect();
    let marker_annotation_columns: Vec<Option<String>> = markers
        .iter()
        .map(|marker| {
            if marker.annotations.is_empty() {
                None
            } else {
                Some(marker.annotations.join("\t"))
            }
        })
        .collect();

    let mut samples_to_process: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(list_path) = &args.input_list {
        samples_to_process = read_sample_list(list_path)?;
    } else if args.paired {
        // Manual paired mode: files are grouped sequentially in pairs
        if args.input.len() % 2 != 0 {
            return Err(AppError::Generic(format!(
                "--paired requires an even number of input files. Found {}.",
                args.input.len()
            )));
        }
        for chunk in args.input.chunks(2) {
            let r1_path = &chunk[0];
            let r2_path = &chunk[1];
            let sample_name = derive_sample_name(r1_path);
            if sample_name != derive_sample_name(r2_path) {
                warn!(
                    "Paired files may not belong to the same sample: {} and {}",
                    r1_path, r2_path
                );
            }
            if samples_to_process
                .insert(sample_name.clone(), vec![r1_path.clone(), r2_path.clone()])
                .is_some()
            {
                return Err(AppError::Generic(format!(
                    "Duplicate sample name '{}' derived from paired FASTQ files. Use explicit --input-list names or rename files.",
                    sample_name
                )));
            }
        }
    } else if !args.no_auto_paired {
        // Auto-detect paired-end files based on naming conventions (_1/_2, _R1/_R2)
        let detection_result = detect_paired_end_files(&args.input);

        if detection_result.is_paired {
            info!(
                "🔗 Auto-detected paired-end files: {} paired samples, {} single-end samples",
                detection_result.paired_count, detection_result.single_count
            );
        }

        samples_to_process = detection_result.samples;
    } else {
        // No auto-detection: treat each file as a separate sample
        for fq_path in &args.input {
            let sample_name = derive_sample_name(fq_path);
            if samples_to_process
                .insert(sample_name.clone(), vec![fq_path.clone()])
                .is_some()
            {
                return Err(AppError::Generic(format!(
                    "Duplicate sample name '{}' derived from FASTQ input. Use explicit --input-list names or rename files.",
                    sample_name
                )));
            }
        }
    }

    info!("Found {} sample(s) to process.", samples_to_process.len());
    debug!("Samples to process: {:?}", samples_to_process.keys());
    let mut all_summary_lines = Vec::new();
    let mut generated_outputs = Vec::new();

    let mut ordered_samples: Vec<(&String, &Vec<String>)> = samples_to_process.iter().collect();
    ordered_samples.sort_by(|(a, _), (b, _)| a.cmp(b));

    let mut used_output_names: HashMap<String, usize> = HashMap::new();

    for (sample_name, fastq_paths) in ordered_samples {
        check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;

        info!("▶ Processing sample: {}", sample_name);
        let counts = match split_kmer::scan_fastq_with_index(
            fastq_paths,
            marker_index.as_ref(),
            bloom_filter_arc.as_ref(),
            markers.len(),
            &args.cancel_token,
            kmer_len,
        ) {
            Ok(counts) => counts,
            Err(e) => {
                cleanup_generated_outputs(&generated_outputs);
                if args
                    .cancel_token
                    .as_ref()
                    .map(|t| t.is_cancelled())
                    .unwrap_or(false)
                    || e.to_string().to_lowercase().contains("cancel")
                {
                    return Err(AppError::Cancelled);
                }
                return Err(AppError::Generic(format!(
                    "Failed to scan FASTQ for {}: {}",
                    sample_name, e
                )));
            }
        };
        check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;
        info!("  Finished scan for {}. Analyzing results...", sample_name);
        trace!("Raw counts for sample {}: {:?}", sample_name, counts);

        let safe_base = sanitize_sample_name_for_path(sample_name);
        let counter = used_output_names.entry(safe_base.clone()).or_insert(0);
        *counter += 1;
        let safe_output_name = if *counter == 1 {
            safe_base
        } else {
            format!("{}_{}", safe_base, counter)
        };
        let detailed_output_path =
            format!("{}_{}_mutations.tsv", args.output_prefix, safe_output_name);
        let mut detailed_writer = match fs::File::create(&detailed_output_path) {
            Ok(f) => BufWriter::new(f),
            Err(e) => {
                cleanup_generated_outputs(&generated_outputs);
                return Err(AppError::Io(e));
            }
        };
        generated_outputs.push(detailed_output_path.clone());

        let header = "pos\tref_allele\talt_allele\tref_count\talt_count\talt_fraction\tlineage_path\textra_annotations...";
        if let Err(e) = writeln!(detailed_writer, "{}", header) {
            drop(detailed_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }

        // This map will store counts for ALL levels defined by each found marker.
        let mut lineage_counts: HashMap<String, usize> = HashMap::new();

        let detail_headers = [
            "pos",
            "ref_allele",
            "alt_allele",
            "ref_count",
            "alt_count",
            "alt_fraction",
            "lineage_path",
            "annotations",
        ];
        let mut detail_excel_writer = if args.excel {
            match crate::common::ExcelStreamWriter::new(&detailed_output_path, &detail_headers) {
                Ok(mut writer) => {
                    // alt_fraction (col 5): red < 50%, yellow 50-95%, green >= 95%
                    writer.set_conditional_formatting(5, 50.0, 95.0);
                    Some(writer)
                }
                Err(e) => {
                    warn!(
                        "  ⚠️ Could not initialize detailed Excel for {}: {}",
                        sample_name, e
                    );
                    None
                }
            }
        } else {
            None
        };

        for (marker_id, &[ref_count, alt_count]) in counts.iter().enumerate() {
            if marker_id % 2048 == 0 {
                if let Err(e) = check_cancelled(&args.cancel_token) {
                    drop(detailed_writer);
                    cleanup_generated_outputs(&generated_outputs);
                    return Err(e);
                }
            }

            let coverage = ref_count + alt_count;
            if coverage < args.min_depth {
                continue;
            }
            let alt_fraction = if coverage > 0 {
                (alt_count as f32 / coverage as f32) * 100.0
            } else {
                0.0
            };

            if alt_fraction >= args.min_alt_percent as f32 {
                let marker = &markers[marker_id];
                let lineage_path_string = &marker_lineage_paths[marker_id];

                // Write detailed report line
                let mut output_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
                    marker.pos + 1,
                    marker.ref_allele,
                    marker.alt_allele,
                    ref_count,
                    alt_count,
                    alt_fraction,
                    lineage_path_string
                );

                // Collect Excel row
                let mut excel_row = vec![
                    (marker.pos + 1).to_string(),
                    marker.ref_allele.to_string(),
                    marker.alt_allele.to_string(),
                    ref_count.to_string(),
                    alt_count.to_string(),
                    format!("{:.2}", alt_fraction),
                    lineage_path_string.clone(),
                ];

                if let Some(annotation_cols) = marker_annotation_columns[marker_id].as_deref() {
                    output_line.push('\t');
                    output_line.push_str(annotation_cols);
                    excel_row.push(annotation_cols.to_string());
                }
                if let Err(e) = writeln!(detailed_writer, "{}", output_line) {
                    drop(detailed_writer);
                    cleanup_generated_outputs(&generated_outputs);
                    return Err(AppError::Io(e));
                }
                let mut disable_detailed_excel = false;
                if let Some(writer) = detail_excel_writer.as_mut() {
                    if let Err(e) = writer.write_row(&excel_row) {
                        warn!(
                            "  ⚠️ Could not append detailed Excel row for {}: {}",
                            sample_name, e
                        );
                        disable_detailed_excel = true;
                    }
                }
                if disable_detailed_excel {
                    detail_excel_writer = None;
                }

                // Add counts for precomputed lineage prefixes to ensure parents are counted.
                for lineage_prefix in &marker_lineage_prefixes[marker_id] {
                    *lineage_counts.entry(lineage_prefix.clone()).or_default() += 1;
                }
            }
        }
        if let Err(e) = detailed_writer.flush() {
            drop(detailed_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }
        drop(detailed_writer);

        // Write detailed Excel file if requested
        if args.excel {
            check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;
            if let Some(writer) = detail_excel_writer.take() {
                match writer.finish() {
                    Ok(xlsx_path) => info!(
                        "  📊 Detailed Excel for {} written to {}",
                        sample_name, xlsx_path
                    ),
                    Err(e) => warn!(
                        "  ⚠️ Could not write detailed Excel for {}: {}",
                        sample_name, e
                    ),
                }
            }
            check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;
        }

        info!(
            "  Detailed report for {} written to {}",
            sample_name, detailed_output_path
        );

        // Generate summary line
        let mut sorted_lineages: Vec<(String, usize)> = lineage_counts
            .iter()
            .map(|(k, v)| (k.clone(), *v))
            .collect();
        sorted_lineages.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let list_str = sorted_lineages
            .iter()
            .map(|(l, c)| format!("{}:{}", l, c))
            .collect::<Vec<_>>()
            .join(" ");

        let major_lineage =
            determine_major_lineage(&lineage_counts, args.nested_classification);

        let summary_line = format!("{}\t{}\t{}", sample_name, list_str, major_lineage);
        all_summary_lines.push(summary_line);
    }

    check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;

    let summary_path = format!("{}_summary.tsv", args.output_prefix);
    info!(
        "▶ Writing final summary for all samples to {}",
        summary_path
    );
    let mut summary_writer = match fs::File::create(&summary_path) {
        Ok(f) => BufWriter::new(f),
        Err(e) => {
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }
    };
    generated_outputs.push(summary_path.clone());
    if let Err(e) = writeln!(summary_writer, "genome\tlineage:count\tmajor_lineage") {
        drop(summary_writer);
        cleanup_generated_outputs(&generated_outputs);
        return Err(AppError::Io(e));
    }

    let summary_headers = ["genome", "lineage:count", "major_lineage"];
    let mut summary_excel_writer = if args.excel {
        match crate::common::ExcelStreamWriter::new(&summary_path, &summary_headers) {
            Ok(writer) => Some(writer),
            Err(e) => {
                warn!("⚠️ Could not initialize summary Excel writer: {}", e);
                None
            }
        }
    } else {
        None
    };

    for (idx, line) in all_summary_lines.iter().enumerate() {
        if idx % 128 == 0 {
            if let Err(e) = check_cancelled(&args.cancel_token) {
                drop(summary_writer);
                cleanup_generated_outputs(&generated_outputs);
                return Err(e);
            }
        }
        if let Err(e) = writeln!(summary_writer, "{}", line) {
            drop(summary_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }
        let mut disable_summary_excel = false;
        if let Some(writer) = summary_excel_writer.as_mut() {
            let row: Vec<String> = line.split('\t').map(String::from).collect();
            if let Err(e) = writer.write_row(&row) {
                warn!("⚠️ Could not append summary Excel row: {}", e);
                disable_summary_excel = true;
            }
        }
        if disable_summary_excel {
            summary_excel_writer = None;
        }
    }
    if let Err(e) = summary_writer.flush() {
        drop(summary_writer);
        cleanup_generated_outputs(&generated_outputs);
        return Err(AppError::Io(e));
    }
    drop(summary_writer);

    // Write summary Excel file if requested
    if args.excel {
        check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;
        if let Some(writer) = summary_excel_writer.take() {
            match writer.finish() {
                Ok(xlsx_path) => info!("📊 Summary Excel file written to {}", xlsx_path),
                Err(e) => warn!("⚠️ Could not write summary Excel file: {}", e),
            }
        }
        check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;
    }

    info!("✅ Process completed.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::lineage::determine_major_lineage;
    use std::collections::HashMap;

    #[test]
    fn nested_call_is_deterministic_on_ties() {
        let mut counts = HashMap::new();
        counts.insert("L2".to_string(), 2usize);
        counts.insert("L1".to_string(), 2usize);
        assert_eq!(determine_major_lineage(&counts, true), "L1");
    }
}
