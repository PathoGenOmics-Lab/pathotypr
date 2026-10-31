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
use crate::split_kmer::{AnchorTally, Marker};
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

    /// Do not read a marker's locus from its anchors when its k-mers could not.
    /// Faster, at the cost of reporting as unread the markers that a second
    /// change inside the same window hides.
    #[arg(long)]
    pub no_anchor_rescue: bool,

    /// [NEW] Enable nested lineage classification logic.
    /// This requires markers to be defined with multiple columns for each level.
    #[arg(long)]
    pub nested_classification: bool,

    /// Size of diagnostic k-mers (1-31). An odd value centres the variant in the
    /// k-mer, so 31 or 21 are the usual choices. Default: 31.
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

    // Recognised once for the whole run, like classify does.
    let resistance_panel = crate::resistance::is_resistance_panel(
        markers.iter().map(|m| m.lineages.join(";")).collect::<Vec<_>>().iter().map(|s| s.as_str()),
    );
    if resistance_panel {
        info!("  The panel is a drug-resistance catalogue; the summary reports the drugs with a WHO grade 1 or 2 marker.");
    }

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

        let header = "pos\tref_allele\talt_allele\tref_count\talt_count\talt_fraction\tlineage_path\tevidence\textra_annotations...";
        if let Err(e) = writeln!(detailed_writer, "{}", header) {
            drop(detailed_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }

        // Markers the reads could not settle. A marker whose locus carries the
        // reference allele is covered by reference reads, so it is absent and
        // known to be. One with no reads at all has not been looked at: a
        // change next to it stops both k-mers matching, and dropping it
        // alongside the genuinely absent ones reports absence of evidence as
        // evidence of absence.
        let undetermined_path =
            format!("{}_{}_undetermined.tsv", args.output_prefix, safe_output_name);
        let mut undetermined_writer = match fs::File::create(&undetermined_path) {
            Ok(f) => BufWriter::new(f),
            Err(e) => {
                drop(detailed_writer);
                cleanup_generated_outputs(&generated_outputs);
                return Err(AppError::Io(e));
            }
        };
        generated_outputs.push(undetermined_path.clone());
        if let Err(e) = writeln!(
            undetermined_writer,
            "pos\tref_allele\talt_allele\tref_count\talt_count\tlineage_path\treason"
        ) {
            drop(detailed_writer);
            drop(undetermined_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }

        // This map will store counts for ALL levels defined by each found marker.
        let mut lineage_counts: HashMap<String, usize> = HashMap::new();
        let mut undetermined_counts: HashMap<String, usize> = HashMap::new();
        let mut resistance: HashMap<String, crate::resistance::Verdict> = HashMap::new();

        let detail_headers = [
            "pos",
            "ref_allele",
            "alt_allele",
            "ref_count",
            "alt_count",
            "alt_fraction",
            "lineage_path",
            "evidence",
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

        let undetermined_headers = [
            "pos",
            "ref_allele",
            "alt_allele",
            "ref_count",
            "alt_count",
            "lineage_path",
            "reason",
        ];
        let mut undetermined_excel_writer = if args.excel {
            match crate::common::ExcelStreamWriter::new(&undetermined_path, &undetermined_headers) {
                Ok(writer) => Some(writer),
                Err(e) => {
                    warn!(
                        "  ⚠️ Could not initialize undetermined Excel for {}: {}",
                        sample_name, e
                    );
                    None
                }
            }
        } else {
            None
        };

        // Which stretches of the genome the reads actually read. Both k-mers of
        // a marker are cut from the same window and matched exactly, so reads
        // on either one read every base of that window: the bases the window
        // holds are known, and a marker lying inside it is absent rather than
        // unread. Without this, each called variant would turn the other
        // markers of its position (ten of them at rpoB 450 in the WHO
        // catalogue) into markers the run claims not to have looked at.
        let read_windows = merge_windows(
            counts
                .iter()
                .enumerate()
                .filter(|(_, &[ref_count, alt_count])| ref_count + alt_count >= args.min_depth)
                .map(|(marker_id, _)| {
                    let start = markers[marker_id].kmer_start;
                    (start, start + kmer_len)
                }),
        );

        // The markers the k-mers could not settle and no neighbour's window
        // covers. Their own anchors may still read the locus, so they are
        // collected before anything is written and the reads are gone.
        let unsettled: Vec<usize> = counts
            .iter()
            .enumerate()
            .filter(|(marker_id, &[ref_count, alt_count])| {
                let marker = &markers[*marker_id];
                ref_count + alt_count < args.min_depth
                    && !window_covers(
                        &read_windows,
                        marker.pos,
                        marker.pos + marker.ref_allele.len(),
                    )
            })
            .map(|(marker_id, _)| marker_id)
            .collect();

        let anchors = if args.no_anchor_rescue || unsettled.is_empty() {
            None
        } else {
            let index = match split_kmer::build_anchor_index(
                &args.reference,
                &markers,
                &unsettled,
                kmer_len,
            ) {
                Ok(index) => index,
                Err(e) => {
                    drop(detailed_writer);
                    drop(undetermined_writer);
                    cleanup_generated_outputs(&generated_outputs);
                    return Err(AppError::Generic(format!(
                        "Failed to build the marker anchors: {}",
                        e
                    )));
                }
            };
            if index.is_empty() {
                None
            } else {
                match split_kmer::scan_fastq_anchors(
                    fastq_paths,
                    &markers,
                    &index,
                    &args.cancel_token,
                ) {
                    Ok(tallies) => Some((index, tallies)),
                    Err(e) => {
                        drop(detailed_writer);
                        drop(undetermined_writer);
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
                            "Failed to read the unsettled loci from their anchors: {}",
                            e
                        )));
                    }
                }
            }
        };
        let mut rescued = 0usize;

        for (marker_id, &kmer_counts) in counts.iter().enumerate() {
            if marker_id % 2048 == 0 {
                if let Err(e) = check_cancelled(&args.cancel_token) {
                    drop(detailed_writer);
                    drop(undetermined_writer);
                    cleanup_generated_outputs(&generated_outputs);
                    return Err(e);
                }
            }

            let marker = &markers[marker_id];
            let path = &marker_lineage_paths[marker_id];
            let call = call_marker(
                marker,
                kmer_counts,
                &read_windows,
                anchors.as_ref().map(|(index, tallies)| MarkerAnchors {
                    readable: index.can_read(marker_id),
                    tally: &tallies[marker_id],
                }),
                args.min_depth,
                args.min_alt_percent,
            );

            let (ref_count, alt_count, evidence) = match call {
                Call::Absent => continue,
                Call::Unread { reference, alternate, reason } => {
                    if let Err(e) = writeln!(
                        undetermined_writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        marker.pos + 1,
                        marker.ref_allele,
                        marker.alt_allele,
                        reference,
                        alternate,
                        path,
                        reason
                    ) {
                        drop(detailed_writer);
                        drop(undetermined_writer);
                        cleanup_generated_outputs(&generated_outputs);
                        return Err(AppError::Io(e));
                    }
                    *undetermined_counts
                        .entry(path.split(';').next().unwrap_or("").to_string())
                        .or_default() += 1;
                    let mut disable_undetermined_excel = false;
                    if let Some(writer) = undetermined_excel_writer.as_mut() {
                        let row = vec![
                            (marker.pos + 1).to_string(),
                            marker.ref_allele.clone(),
                            marker.alt_allele.clone(),
                            reference.to_string(),
                            alternate.to_string(),
                            path.clone(),
                            reason.to_string(),
                        ];
                        if let Err(e) = writer.write_row(&row) {
                            warn!(
                                "  ⚠️ Could not append undetermined Excel row for {}: {}",
                                sample_name, e
                            );
                            disable_undetermined_excel = true;
                        }
                    }
                    if disable_undetermined_excel {
                        undetermined_excel_writer = None;
                    }
                    continue;
                }
                Call::Present { reference, alternate, evidence } => {
                    (reference, alternate, evidence)
                }
            };
            if evidence != "kmer" {
                rescued += 1;
            }

            let coverage = ref_count + alt_count;
            let alt_fraction = if coverage > 0 {
                (alt_count as f32 / coverage as f32) * 100.0
            } else {
                0.0
            };

            // Write detailed report line
            let mut output_line = format!(
                "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
                marker.pos + 1,
                marker.ref_allele,
                marker.alt_allele,
                ref_count,
                alt_count,
                alt_fraction,
                path,
                evidence
            );

            // Collect Excel row
            let mut excel_row = vec![
                (marker.pos + 1).to_string(),
                marker.ref_allele.to_string(),
                marker.alt_allele.to_string(),
                ref_count.to_string(),
                alt_count.to_string(),
                format!("{:.2}", alt_fraction),
                path.clone(),
                evidence.to_string(),
            ];

            if let Some(annotation_cols) = marker_annotation_columns[marker_id].as_deref() {
                output_line.push('\t');
                output_line.push_str(annotation_cols);
                excel_row.push(annotation_cols.to_string());
            }
            if let Err(e) = writeln!(detailed_writer, "{}", output_line) {
                drop(detailed_writer);
                drop(undetermined_writer);
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
            crate::resistance::fold_marker(&mut resistance, path);
        }
        if rescued > 0 {
            info!(
                "  {} markers of {} were read from their anchors, which their k-mers could not",
                rescued, sample_name
            );
        }
        if let Err(e) = detailed_writer.flush() {
            drop(detailed_writer);
            drop(undetermined_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }
        drop(detailed_writer);
        if let Err(e) = undetermined_writer.flush() {
            drop(undetermined_writer);
            cleanup_generated_outputs(&generated_outputs);
            return Err(AppError::Io(e));
        }
        drop(undetermined_writer);

        // Write detailed Excel file if requested
        if args.excel {
            check_cancelled_with_cleanup(&args.cancel_token, &generated_outputs)?;
            if let Some(writer) = undetermined_excel_writer.take() {
                match writer.finish() {
                    Ok(xlsx_path) => info!(
                        "  📊 Undetermined Excel for {} written to {}",
                        sample_name, xlsx_path
                    ),
                    Err(e) => warn!(
                        "  ⚠️ Could not write undetermined Excel for {}: {}",
                        sample_name, e
                    ),
                }
            }
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

        // The WHO catalogue is not a hierarchy, so tallying its paths answers
        // which drug has the most markers, and most of that catalogue is
        // graded as not associated with resistance.
        let major_lineage = if resistance_panel {
            crate::resistance::call(&resistance)
        } else {
            determine_major_lineage(&lineage_counts, args.nested_classification)
        };
        let mut sorted_undetermined: Vec<(&String, &usize)> = undetermined_counts.iter().collect();
        sorted_undetermined.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        let undetermined_total: usize = undetermined_counts.values().sum();
        let undetermined_str = sorted_undetermined
            .iter()
            .map(|(l, c)| format!("{}:{}", l, c))
            .collect::<Vec<_>>()
            .join(" ");
        if undetermined_total > 0 {
            info!(
                "  {} markers of {} could be read neither way; they are listed in {}",
                undetermined_total, sample_name, undetermined_path
            );
        }

        let summary_line = format!(
            "{}\t{}\t{}\t{}\t{}",
            sample_name, list_str, major_lineage, undetermined_total, undetermined_str
        );
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
    if let Err(e) = writeln!(
        summary_writer,
        "genome\tlineage:count\tmajor_lineage\tundetermined\tundetermined:count"
    ) {
        drop(summary_writer);
        cleanup_generated_outputs(&generated_outputs);
        return Err(AppError::Io(e));
    }

    let summary_headers = [
        "genome",
        "lineage:count",
        "major_lineage",
        "undetermined",
        "undetermined:count",
    ];
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

/// What one marker turned out to be in one sample.
enum Call {
    /// Read, and not carrying the variant.
    Absent,
    /// Carrying it, with the counts behind the call and where they came from.
    Present { reference: u32, alternate: u32, evidence: &'static str },
    /// Not read, with the reason for the undetermined table.
    Unread { reference: u32, alternate: u32, reason: &'static str },
}

/// What the anchor pass has on one marker.
struct MarkerAnchors<'a> {
    /// Whether any of its anchors is unique in the reference. False means none
    /// of them can say which copy of the locus a read came from.
    readable: bool,
    /// Its left and right anchors, in that order.
    tally: &'a [AnchorTally; 2],
}

/// Whether reads on one anchor say the variant is there.
fn anchor_says_present(tally: &AnchorTally, min_alt_percent: u32) -> bool {
    let depth = tally.depth();
    depth > 0 && (tally.alternate as f32 / depth as f32) * 100.0 >= min_alt_percent as f32
}

/// Decides one marker from the evidence: first the k-mers, then the windows the
/// reads read, then the marker's own anchors.
///
/// A k-mer is an exact match, so it settles a marker whenever it matches at
/// all. When neither k-mer does, the locus may still have been read: by a
/// neighbour whose window covers it, or from one of this marker's anchors,
/// which survives a change on the other side of the variant.
fn call_marker(
    marker: &Marker,
    kmer_counts: [u32; 2],
    read_windows: &[(usize, usize)],
    anchors: Option<MarkerAnchors>,
    min_depth: u32,
    min_alt_percent: u32,
) -> Call {
    let [reference, alternate] = kmer_counts;
    let coverage = reference + alternate;
    if coverage >= min_depth {
        let fraction = (alternate as f32 / coverage as f32) * 100.0;
        return if fraction >= min_alt_percent as f32 {
            Call::Present { reference, alternate, evidence: "kmer" }
        } else {
            Call::Absent
        };
    }

    if window_covers(read_windows, marker.pos, marker.pos + marker.ref_allele.len()) {
        return Call::Absent;
    }

    let depth_reason = if coverage == 0 { "not_found" } else { "low_depth" };
    let Some(anchors) = anchors else {
        return Call::Unread { reference, alternate, reason: depth_reason };
    };
    if !anchors.readable {
        return Call::Unread { reference, alternate, reason: "no_unique_anchor" };
    }
    let tally = anchors.tally;

    // A read spanning the whole locus matches both anchors, so the two sides
    // cannot be added up. The deeper one answers, and the other has to agree if
    // it was read deeply enough to have an opinion.
    let (chosen, other, side) = if tally[0].depth() >= tally[1].depth() {
        (&tally[0], &tally[1], "left_anchor")
    } else {
        (&tally[1], &tally[0], "right_anchor")
    };
    if chosen.depth() < min_depth {
        let reason = if chosen.depth() == 0 { depth_reason } else { "low_depth" };
        return Call::Unread { reference: chosen.reference, alternate: chosen.alternate, reason };
    }
    let present = anchor_says_present(chosen, min_alt_percent);
    let evidence = if other.depth() >= min_depth {
        if anchor_says_present(other, min_alt_percent) != present {
            return Call::Unread {
                reference: chosen.reference,
                alternate: chosen.alternate,
                reason: "anchors_disagree",
            };
        }
        "both_anchors"
    } else {
        side
    };
    if present {
        Call::Present { reference: chosen.reference, alternate: chosen.alternate, evidence }
    } else {
        Call::Absent
    }
}

/// The stretches of the reference the reads read, as sorted disjoint
/// half-open intervals. Windows that touch or overlap are joined, because a
/// base is known if any one window holds it.
fn merge_windows(windows: impl Iterator<Item = (usize, usize)>) -> Vec<(usize, usize)> {
    let mut windows: Vec<(usize, usize)> = windows.collect();
    windows.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(windows.len());
    for (start, end) in windows {
        match merged.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    merged
}

/// Whether every base of `start .. end` falls inside one read window.
fn window_covers(merged: &[(usize, usize)], start: usize, end: usize) -> bool {
    // The last window that begins at or before the span.
    match merged.partition_point(|&(window_start, _)| window_start <= start) {
        0 => false,
        n => merged[n - 1].1 >= end,
    }
}


#[cfg(test)]
mod tests {
    use super::{
        call_marker, merge_windows, window_covers, AnchorTally, Call, Marker, MarkerAnchors,
    };
    use crate::lineage::determine_major_lineage;
    use std::collections::HashMap;

    #[test]
    fn nested_call_is_deterministic_on_ties() {
        let mut counts = HashMap::new();
        counts.insert("L2".to_string(), 2usize);
        counts.insert("L1".to_string(), 2usize);
        assert_eq!(determine_major_lineage(&counts, true), "L1");
    }

    /// A marker with nothing but its position and alleles, which is all the
    /// decision looks at.
    fn marker_at(pos: usize, reference: &str, alt: &str) -> Marker {
        Marker {
            pos,
            ref_allele: reference.to_string(),
            alt_allele: alt.to_string(),
            lineages: vec!["RIF".to_string()],
            alt_kmer: Vec::new(),
            ref_kmer: Vec::new(),
            kmer_start: pos,
            annotations: Vec::new(),
        }
    }

    fn tally(reference: u32, alternate: u32) -> AnchorTally {
        AnchorTally { reference, alternate, other: 0 }
    }

    /// A k-mer is exact, so when one matches at all it settles the marker and
    /// the anchors are never asked.
    #[test]
    fn the_kmers_settle_a_marker_they_could_read() {
        let marker = marker_at(75, "G", "A");
        // Enough reads on the alternate k-mer: present, on the k-mer's word.
        let call = call_marker(&marker, [0, 20], &[], Some(MarkerAnchors { readable: true, tally: &[tally(0, 0); 2] }), 5, 95);
        assert!(matches!(
            call,
            Call::Present { reference: 0, alternate: 20, evidence: "kmer" }
        ));
        // Enough reads, all on the reference: absent, and not reported at all.
        let call = call_marker(&marker, [20, 0], &[], Some(MarkerAnchors { readable: true, tally: &[tally(0, 20); 2] }), 5, 95);
        assert!(matches!(call, Call::Absent), "the anchors must not overrule a read k-mer");
    }

    /// Read through a neighbour's window: absent, without asking the anchors.
    #[test]
    fn a_marker_inside_a_window_the_reads_read_is_absent() {
        let marker = marker_at(75, "G", "A");
        let windows = [(60usize, 100usize)];
        let call = call_marker(&marker, [0, 0], &windows, Some(MarkerAnchors { readable: true, tally: &[tally(0, 30); 2] }), 5, 95);
        assert!(matches!(call, Call::Absent));
    }

    /// Without the anchor pass the depth is all there is to go on, and the
    /// marker is reported unread rather than absent.
    #[test]
    fn without_anchors_an_unread_marker_says_so() {
        let marker = marker_at(75, "G", "A");
        assert!(matches!(
            call_marker(&marker, [0, 0], &[], None, 5, 95),
            Call::Unread { reason: "not_found", .. }
        ));
        assert!(matches!(
            call_marker(&marker, [2, 1], &[], None, 5, 95),
            Call::Unread { reason: "low_depth", .. }
        ));
    }

    /// The anchor reads the locus the k-mers could not: this is the call that
    /// used to be missed, and the one that used to be reported as absent.
    #[test]
    fn an_anchor_calls_a_marker_the_kmers_missed() {
        let marker = marker_at(75, "G", "A");
        // Only the left anchor was read, and it carries the variant.
        let call = call_marker(
            &marker,
            [0, 0],
            &[],
            Some(MarkerAnchors { readable: true, tally: &[tally(0, 112), tally(0, 0)] }),
            5,
            95,
        );
        assert!(matches!(
            call,
            Call::Present { reference: 0, alternate: 112, evidence: "left_anchor" }
        ));
        // The right anchor alone, reading the reference: absent, not unread.
        let call = call_marker(
            &marker,
            [0, 0],
            &[],
            Some(MarkerAnchors { readable: true, tally: &[tally(0, 0), tally(40, 0)] }),
            5,
            95,
        );
        assert!(matches!(call, Call::Absent));
        // Both anchors, agreeing. A read spanning the locus matches both, so
        // the deeper one answers and the counts are not added up.
        let call = call_marker(
            &marker,
            [0, 0],
            &[],
            Some(MarkerAnchors { readable: true, tally: &[tally(0, 30), tally(0, 25)] }),
            5,
            95,
        );
        assert!(matches!(
            call,
            Call::Present { reference: 0, alternate: 30, evidence: "both_anchors" }
        ));
    }

    /// Two anchors that disagree are two anchors nobody should believe.
    #[test]
    fn anchors_that_disagree_leave_the_marker_unread() {
        let marker = marker_at(75, "G", "A");
        let call = call_marker(
            &marker,
            [0, 0],
            &[],
            Some(MarkerAnchors { readable: true, tally: &[tally(0, 30), tally(25, 0)] }),
            5,
            95,
        );
        assert!(matches!(call, Call::Unread { reason: "anchors_disagree", .. }));
        // The shallow side has no opinion, so it does not create a conflict.
        let call = call_marker(
            &marker,
            [0, 0],
            &[],
            Some(MarkerAnchors { readable: true, tally: &[tally(0, 30), tally(2, 0)] }),
            5,
            95,
        );
        assert!(matches!(call, Call::Present { evidence: "left_anchor", .. }));
    }

    /// An anchor the reference repeats cannot say which copy a read came from,
    /// and the marker is reported unread with that as the reason.
    #[test]
    fn a_marker_whose_anchors_are_repeated_says_which_reason() {
        let marker = marker_at(75, "G", "A");
        let call = call_marker(&marker, [0, 0], &[], Some(MarkerAnchors { readable: false, tally: &[tally(0, 0); 2] }), 5, 95);
        assert!(matches!(call, Call::Unread { reason: "no_unique_anchor", .. }));
    }

    /// Too few reads on the anchor is still too few reads, and the count
    /// reported is the one the decision rests on.
    #[test]
    fn a_shallow_anchor_reports_its_own_depth() {
        let marker = marker_at(75, "G", "A");
        let call = call_marker(
            &marker,
            [0, 0],
            &[],
            Some(MarkerAnchors { readable: true, tally: &[tally(1, 2), tally(0, 0)] }),
            5,
            95,
        );
        assert!(matches!(
            call,
            Call::Unread { reference: 1, alternate: 2, reason: "low_depth" }
        ));
    }

    /// A base is known if any read window holds it, so touching windows join.
    #[test]
    fn read_windows_are_joined_and_sorted() {
        let merged = merge_windows([(100, 131), (40, 71), (120, 151), (71, 90)].into_iter());
        assert_eq!(merged, vec![(40, 90), (100, 151)]);
        assert!(merge_windows(std::iter::empty()).is_empty());
        // A window inside another one must not shorten it.
        assert_eq!(
            merge_windows([(100, 131), (105, 120)].into_iter()),
            vec![(100, 131)]
        );
    }

    /// The case this exists for: rpoB 450 carries ten markers in the WHO
    /// catalogue. One of them is called, which reads the window for all of
    /// them, so the other nine are absent and not reported as unread. A marker
    /// hanging over the edge of the window is still unread.
    #[test]
    fn a_marker_inside_a_read_window_is_not_reported_as_unread() {
        let merged = merge_windows([(761140, 761171)].into_iter());
        assert!(window_covers(&merged, 761154, 761155), "the called SNP");
        assert!(window_covers(&merged, 761154, 761156), "its two-base neighbour");
        assert!(window_covers(&merged, 761154, 761165), "the eleven-base marker");
        assert!(window_covers(&merged, 761168, 761171), "ends on the last base read");
        assert!(window_covers(&merged, 761140, 761143), "starts on the first base read");
        assert!(!window_covers(&merged, 761154, 761180), "runs past the window");
        assert!(!window_covers(&merged, 761168, 761172), "one base past the last");
        assert!(!window_covers(&merged, 761130, 761145), "starts before it");
        assert!(!window_covers(&merged, 900000, 900001), "nowhere near it");
        assert!(!window_covers(&[], 761154, 761155), "nothing was read");
    }

    /// Two windows side by side do not vouch for a span that crosses the gap
    /// between them.
    #[test]
    fn a_gap_between_windows_leaves_a_span_unread() {
        let merged = merge_windows([(100, 131), (140, 171)].into_iter());
        assert!(window_covers(&merged, 125, 130));
        assert!(!window_covers(&merged, 128, 145));
        assert!(!window_covers(&merged, 131, 141));
    }
}
