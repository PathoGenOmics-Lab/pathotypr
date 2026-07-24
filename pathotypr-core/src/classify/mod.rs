//! Classify subcommand — marker-based variant calling on assemblies.
//!
//! Submodules:
//! - `markers` — Marker TSV parsing, k-mer generation, genome scanning
//! - `annotation` — GFF parsing and variant-to-amino-acid translation
//! - `masking` — FASTA sequence masking at marker positions

pub mod annotation;
pub mod markers;
pub mod masking;

// Re-export public types for backward compatibility with Tauri and other modules.
pub use annotation::{Gene, Strand, parse_gff_and_build_tree};
pub use markers::{
    find_markers, generate_markerkmer, get_positions, MarkerIndex, MarkerKmerEntry, MarkerMatch,
    MarkerVariant,
};
pub use masking::{collect_mask_positions, mask_sequence, write_masked_fasta};

use annotation::translate_variant_info;

use crate::common::configure_thread_pool;
use crate::errors::{
    check_cancelled, AppError, AppResult, CancellationToken, ParallelCancellation,
};
use crate::lineage::determine_major_lineage;
use needletail::Sequence;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use rayon::prelude::*;
use rust_lapper::Lapper;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

/// Command-line arguments for the classify subcommand.
#[derive(Parser, Debug)]
pub struct Args {
    /// Path to a TSV file defining markers.
    #[arg(short = 'm', long = "markers")]
    pub tsv_pos: String,

    /// Reference FASTA file used to define the markers.
    #[arg(short = 'r', long = "reference")]
    pub ref_fasta: String,

    /// Path to a TSV file listing samples. Format: sample_name\tfasta_path[\tgff_path].
    #[arg(short = 'l', long = "input-list", required_unless_present = "fasta_genomes")]
    pub tsv_genomes: Option<String>,

    /// Single FASTA file to analyze.
    #[arg(short = 'i', long = "input", required_unless_present_any = ["tsv_genomes", "fasta_files"])]
    pub fasta_genomes: Option<String>,

    /// Multiple FASTA files to analyze (for GUI batch mode).
    #[arg(long = "input-files", required_unless_present_any = ["tsv_genomes", "fasta_genomes"])]
    pub fasta_files: Option<Vec<String>>,

    /// Optional GFF file for annotation when using --input.
    #[arg(long = "gff", requires = "fasta_genomes")]
    pub gff_file: Option<String>,

    /// Multiple GFF files for annotation when using --input-files (matched by filename).
    #[arg(long = "gff-files")]
    pub gff_files: Option<Vec<String>>,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "output-prefix")]
    pub ofile: String,

    /// k-mer size (default is 31)
    #[arg(long, default_value_t = 31)]
    pub kmer_size: usize,

    /// Number of threads for parallel processing.
    #[arg(short = 't', long = "threads")]
    pub num_cpu: Option<usize>,

    /// Enable nested lineage classification logic.
    #[arg(long)]
    pub nested_classification: bool,

    /// Minimum flanking bases on each side of the allele in the marker k-mer.
    #[arg(long, default_value_t = 10)]
    pub min_flank_bases: usize,

    /// Output masked FASTA files with marker positions replaced by N.
    #[arg(long)]
    pub output_masked_fasta: bool,

    /// Also generate Excel (.xlsx) files alongside TSV outputs.
    #[arg(long, default_value = "false")]
    pub excel: bool,

    /// Cancellation token for stopping the task (GUI only, not CLI).
    #[arg(skip)]
    pub cancel_token: Option<CancellationToken>,
}

// ---------------------------------------------------------------------------
// Reverse complement
// ---------------------------------------------------------------------------

pub fn reverse_complement_sequence(sequence: &str) -> String {
    static COMPLEMENT: [u8; 256] = {
        let mut table = [b'N'; 256];
        table[b'A' as usize] = b'T';
        table[b'T' as usize] = b'A';
        table[b'C' as usize] = b'G';
        table[b'G' as usize] = b'C';
        table[b'a' as usize] = b'T';
        table[b't' as usize] = b'A';
        table[b'c' as usize] = b'G';
        table[b'g' as usize] = b'C';
        table
    };
    let bytes = sequence.as_bytes();
    let mut result = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        result.push(COMPLEMENT[b as usize]);
    }
    String::from_utf8(result).expect("complement table produced invalid UTF-8")
}

// ---------------------------------------------------------------------------
// Reference reading
// ---------------------------------------------------------------------------

pub fn get_ref(fasta_file: &str) -> AppResult<String> {
    let mut reader = needletail::parse_fastx_file(fasta_file).map_err(|e| {
        AppError::Generic(format!("Failed to open FASTA file {}: {}", fasta_file, e))
    })?;
    let record = reader.next().ok_or_else(|| {
        AppError::NotEnoughData("No record found in the reference FASTA file.".to_string())
    })?.map_err(|e| {
        AppError::Generic(format!("Failed to read record from {}: {}", fasta_file, e))
    })?;
    let seq_bytes = record.normalize(true).into_owned();
    if reader.next().is_some() {
        return Err(AppError::Generic(format!(
            "Reference FASTA '{}' contains multiple records; provide a single-record FASTA.",
            fasta_file
        )));
    }
    String::from_utf8(seq_bytes).map_err(|e| {
        AppError::Parsing(format!("Invalid UTF-8 in FASTA: {}", e))
    })
}

// ---------------------------------------------------------------------------
// Genome paths from TSV
// ---------------------------------------------------------------------------

fn get_genomepaths(tsv_file: &str) -> AppResult<HashMap<String, (String, Option<String>)>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_path(tsv_file)?;
    let mut genome_paths = HashMap::new();
    for result in rdr.records() {
        let record = result?;
        if record.len() < 2 { continue; }
        let col0 = record.get(0).map(str::trim).unwrap_or("");
        let col1 = record.get(1).map(str::trim).unwrap_or("");
        if col0.starts_with('#') || col0.is_empty() { continue; }
        if (col0.eq_ignore_ascii_case("sample") || col0.eq_ignore_ascii_case("sample_name") || col0.eq_ignore_ascii_case("genome"))
            && (col1.contains("fasta") || col1.contains("path"))
        {
            continue;
        }
        let genome_name = record[0].trim().to_string();
        if genome_name.is_empty() {
            warn!("Skipping row with empty genome/sample name in {}", tsv_file);
            continue;
        }
        let fasta_path = record[1].trim().to_string();
        let gff_path = record.get(2).and_then(|p| {
            let trimmed = p.trim();
            if trimmed.is_empty() { None } else { Some(trimmed.to_string()) }
        });
        if !Path::new(&fasta_path).exists() {
            return Err(AppError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("The FASTA file '{}' for sample '{}' does not exist", fasta_path, genome_name),
            )));
        }
        if let Some(gff) = &gff_path {
            if !Path::new(gff).exists() {
                return Err(AppError::Io(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    format!("The GFF file '{}' for sample '{}' does not exist", gff, genome_name),
                )));
            }
        }
        if genome_paths.insert(genome_name.clone(), (fasta_path, gff_path)).is_some() {
            return Err(AppError::Generic(format!(
                "Duplicate sample/genome name '{}' in input list '{}'. Sample names must be unique.",
                genome_name, tsv_file
            )));
        }
    }
    if genome_paths.is_empty() {
        return Err(AppError::NotEnoughData(format!(
            "No valid genome entries found in '{}'.", tsv_file
        )));
    }
    Ok(genome_paths)
}

pub fn get_genomes_from_fasta(fasta_file: &str) -> AppResult<Vec<(String, String)>> {
    let mut reader = needletail::parse_fastx_file(fasta_file).map_err(|e| {
        AppError::Generic(format!("Failed to open FASTA file {}: {}", fasta_file, e))
    })?;
    let mut genomes = Vec::new();
    let mut seen_ids = HashSet::new();
    while let Some(result) = reader.next() {
        let record = result.map_err(|e| {
            AppError::Generic(format!("Failed to read record from {}: {}", fasta_file, e))
        })?;
        let genome_id = String::from_utf8_lossy(record.id()).to_string();
        if !seen_ids.insert(genome_id.clone()) {
            return Err(AppError::Generic(format!(
                "Duplicate sequence ID '{}' found in FASTA '{}'. Sequence IDs must be unique.",
                genome_id, fasta_file
            )));
        }
        let seq = String::from_utf8(record.normalize(true).into_owned())
            .map_err(|e| AppError::Parsing(format!("Invalid UTF-8 in FASTA: {}", e)))?;
        genomes.push((genome_id, seq));
    }
    Ok(genomes)
}

// ---------------------------------------------------------------------------
// Genome analysis
// ---------------------------------------------------------------------------

fn analyze_genome(
    genome_name: &str,
    fasta_path: &str,
    gff_path: &Option<String>,
    shared_annotations: Option<&Lapper<usize, Gene>>,
    marker_index: &MarkerIndex,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> AppResult<Vec<String>> {
    if !Path::new(fasta_path).exists() {
        return Err(AppError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("The FASTA file '{}' does not exist", fasta_path),
        )));
    }
    info!("Analyzing genome {} ({})", genome_name, fasta_path);
    // Read every record so multi-contig draft assemblies (the common case for
    // bacterial WGS) are supported. `get_ref` is reserved for the single-record
    // reference genome; sample assemblies routinely contain several contigs.
    let contigs = get_genomes_from_fasta(fasta_path)?;
    if contigs.is_empty() {
        return Err(AppError::NotEnoughData(format!(
            "No records found in FASTA '{}'.",
            fasta_path
        )));
    }

    let owned_annotations;
    let annotations: Option<&Lapper<usize, Gene>> = if shared_annotations.is_some() {
        shared_annotations
    } else if let Some(gff) = gff_path {
        owned_annotations = parse_gff_and_build_tree(gff)?;
        Some(&owned_annotations)
    } else {
        None
    };

    Ok(collect_marker_lines(
        genome_name,
        contigs.iter().map(|(_, seq)| seq.as_str()),
        marker_index,
        &annotations,
        k,
        ref_seq,
        ref_seq_rc,
    ))
}

/// Collects the detail lines for one genome, scanning every contig on **both
/// strands** and reporting each marker only once.
///
/// Contig orientation in a draft assembly is arbitrary, so a marker sitting on
/// a reverse-oriented contig only matches the reverse complement; scanning the
/// forward strand alone silently loses it (the FASTQ path already indexes both
/// orientations). Matches are deduplicated by marker k-mer, so a marker seen on
/// several contigs — or on both strands of one contig — counts as a single
/// observation rather than inflating the lineage tally.
fn collect_marker_lines<'a>(
    genome_name: &str,
    contigs: impl Iterator<Item = &'a str>,
    marker_index: &MarkerIndex,
    annotations: &Option<&Lapper<usize, Gene>>,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> Vec<String> {
    let mut result = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for contig_seq in contigs {
        let contig_rc = reverse_complement_sequence(contig_seq);
        for seq in [contig_seq, contig_rc.as_str()] {
            for m in find_markers(seq, marker_index, k) {
                if !seen.insert(m.kmer.clone()) {
                    continue;
                }
                // `seq` is the strand that produced the match, so the positions
                // and the ALT bases sliced out of it stay self-consistent.
                result.extend(format_marker_match(
                    genome_name,
                    &m,
                    annotations,
                    seq,
                    ref_seq,
                    ref_seq_rc,
                ));
            }
        }
    }
    if result.is_empty() {
        result.push(format!("{}\t\t\t\t\t\t\t\t\t\t\t\t\n", genome_name));
    }
    result
}

pub fn analyze_genome_seq(
    genome_name: &str,
    genome_seq: &str,
    marker_index: &MarkerIndex,
    annotations: &Option<Lapper<usize, Gene>>,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> Vec<String> {
    info!("Analyzing genome {} (from provided sequence)", genome_name);
    let ann_ref: Option<&Lapper<usize, Gene>> = annotations.as_ref();
    collect_marker_lines(
        genome_name,
        std::iter::once(genome_seq),
        marker_index,
        &ann_ref,
        k,
        ref_seq,
        ref_seq_rc,
    )
}

fn format_marker_match(
    genome_name: &str,
    m: &MarkerMatch,
    annotations: &Option<&Lapper<usize, Gene>>,
    genome_seq: &str,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> Vec<String> {
    let ref_pos_0based = m.ref_position.saturating_sub(1);
    let alt_len = m.alt_allele_len;
    let ref_len = m.ref_allele_len;
    let variant_start = m.genome_position + m.left_flank_len;
    let variant_end = (variant_start + alt_len).min(genome_seq.len());
    let alt_bases = &genome_seq[variant_start..variant_end];

    let (gff_gene_id, gene_start, gene_end, gff_aa_pos, gff_aa_change) = if let Some(tree) = annotations {
        if let Some(gene_interval) = tree.find(ref_pos_0based, ref_pos_0based + ref_len).next() {
            let (gene_id, aa_pos, aa_change) = translate_variant_info(
                gene_interval, ref_pos_0based, alt_bases, ref_len, ref_seq, ref_seq_rc,
            );
            (
                gene_id,
                Some((gene_interval.start + 1).to_string()),
                Some(gene_interval.stop.to_string()),
                aa_pos,
                aa_change,
            )
        } else {
            (None, None, None, None, None)
        }
    } else {
        (None, None, None, None, None)
    };

    let gene_id = gff_gene_id.or_else(|| m.gene.clone());
    let aa_pos = gff_aa_pos;
    let aa_change = gff_aa_change.or_else(|| m.mutation.clone());

    vec![format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        genome_name,
        m.kmer,
        m.genome_position + 1,
        variant_start + 1,
        m.ref_position,
        m.ref_allele,
        m.alt_allele,
        m.lineage,
        gene_id.as_deref().unwrap_or(""),
        gene_start.as_deref().unwrap_or(""),
        gene_end.as_deref().unwrap_or(""),
        aa_pos.as_deref().unwrap_or(""),
        aa_change.as_deref().unwrap_or("")
    )]
}

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------

fn check_input_files(args: &Args) -> AppResult<()> {
    let all_files = vec![
        (args.tsv_pos.as_str(), "TSV positions file"),
        (args.ref_fasta.as_str(), "Reference FASTA file"),
        (args.tsv_genomes.as_deref().unwrap_or(""), "TSV genomes file"),
        (args.fasta_genomes.as_deref().unwrap_or(""), "FASTA genomes file"),
    ];
    for (path, description) in all_files.into_iter().filter(|(p, _)| !p.is_empty()) {
        if !Path::new(path).exists() {
            return Err(AppError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("The {} {} does not exist", description, path),
            )));
        }
    }
    Ok(())
}

fn find_matching_gff<'a>(fasta_file: &str, gff_files: &'a [String]) -> Option<&'a String> {
    let fasta_stem = Path::new(fasta_file).file_stem().and_then(|s| s.to_str())?;
    gff_files.iter().find(|gff_file| {
        Path::new(gff_file)
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|gff_stem| gff_stem == fasta_stem)
            .unwrap_or(false)
    })
}

// ---------------------------------------------------------------------------
// Genome processing orchestration
// ---------------------------------------------------------------------------

const PARALLEL_PROCESS_CHUNK_SIZE: usize = 256;

fn process_genomes<F>(
    args: &Args,
    marker_index: &MarkerIndex,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
    mut emit_line: F,
) -> AppResult<()>
where
    F: FnMut(&str) -> AppResult<()>,
{
    let cancel_token = &args.cancel_token;

    if let Some(tsv_genomes) = &args.tsv_genomes {
        let genome_paths_vec: Vec<(String, (String, Option<String>))> =
            get_genomepaths(tsv_genomes)?.into_iter().collect();
        debug!("Processing {} genomes from input list.", genome_paths_vec.len());
        let pb = ProgressBar::new(genome_paths_vec.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")?
                .progress_chars("=>-"),
        );
        let cancellation = ParallelCancellation::new(cancel_token);
        let fallback_gff = &args.gff_file;
        let shared_gff_tree = if let Some(gff) = fallback_gff {
            Some(parse_gff_and_build_tree(gff)?)
        } else {
            None
        };
        for chunk in genome_paths_vec.chunks(PARALLEL_PROCESS_CHUNK_SIZE) {
            let res_nested: AppResult<Vec<Vec<String>>> = chunk
                .par_iter()
                .map(|(genome_name, (fasta_path, gff_path))| {
                    if cancellation.is_cancelled() { return Ok(Vec::new()); }
                    let (effective_gff, shared_ann) = if gff_path.is_some() {
                        (gff_path, None)
                    } else {
                        (&None, shared_gff_tree.as_ref())
                    };
                    let lines = analyze_genome(genome_name, fasta_path, effective_gff, shared_ann, marker_index, k, ref_seq, ref_seq_rc)?;
                    pb.inc(1);
                    Ok(lines)
                })
                .collect();
            cancellation.check_after()?;
            for genome_lines in res_nested? {
                for line in genome_lines {
                    emit_line(&line)?;
                }
            }
        }
        pb.finish_with_message("Done!");
        cancellation.check_after()?;
    } else if let Some(fasta_genomes) = &args.fasta_genomes {
        debug!("Processing genomes from single FASTA input: {}", fasta_genomes);
        let annotations = if let Some(gff_path) = &args.gff_file {
            Some(parse_gff_and_build_tree(gff_path)?)
        } else {
            None
        };
        let mut reader = needletail::parse_fastx_file(fasta_genomes).map_err(|e| {
            AppError::Generic(format!("Failed to open FASTA file {}: {}", fasta_genomes, e))
        })?;
        let mut seen_ids: HashSet<String> = HashSet::new();
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::with_template("[{elapsed_precise}] {spinner:.blue} Genomes processed: {pos}")?);
        pb.enable_steady_tick(std::time::Duration::from_millis(120));
        let cancellation = ParallelCancellation::new(cancel_token);
        loop {
            let mut batch: Vec<(String, String)> = Vec::with_capacity(PARALLEL_PROCESS_CHUNK_SIZE);
            for _ in 0..PARALLEL_PROCESS_CHUNK_SIZE {
                match reader.next() {
                    Some(Ok(record)) => {
                        let genome_id = String::from_utf8_lossy(record.id()).to_string();
                        if !seen_ids.insert(genome_id.clone()) {
                            return Err(AppError::Generic(format!(
                                "Duplicate sequence ID '{}' found in FASTA '{}'. Sequence IDs must be unique.",
                                genome_id, fasta_genomes
                            )));
                        }
                        let seq = String::from_utf8(record.normalize(true).into_owned())
                            .map_err(|e| AppError::Parsing(format!("Invalid UTF-8 in FASTA: {}", e)))?;
                        batch.push((genome_id, seq));
                    }
                    Some(Err(e)) => {
                        return Err(AppError::Generic(format!(
                            "Failed to read record from {}: {}", fasta_genomes, e
                        )));
                    }
                    None => break,
                }
            }
            if batch.is_empty() { break; }
            let chunk_results: Vec<Vec<String>> = batch
                .par_iter()
                .map(|(genome_name, genome_seq)| {
                    if cancellation.is_cancelled() { return Vec::new(); }
                    let lines = analyze_genome_seq(genome_name, genome_seq, marker_index, &annotations, k, ref_seq, ref_seq_rc);
                    pb.inc(1);
                    lines
                })
                .collect();
            cancellation.check_after()?;
            for genome_lines in chunk_results {
                for line in genome_lines {
                    emit_line(&line)?;
                }
            }
        }
        pb.finish_with_message("Done!");
        cancellation.check_after()?;
    } else if let Some(fasta_files) = &args.fasta_files {
        info!("📂 Processing {} FASTA files...", fasta_files.len());
        let gff_list: Vec<String> = args.gff_files.clone().unwrap_or_default();
        if !gff_list.is_empty() {
            info!("📎 {} GFF annotation file(s) provided for matching", gff_list.len());
        }
        let shared_gff: Option<Arc<Option<Lapper<usize, Gene>>>> =
            if gff_list.len() == 1 {
                info!("📎 Pre-parsing shared GFF: {}", gff_list[0]);
                Some(Arc::new(Some(parse_gff_and_build_tree(&gff_list[0])?)))
            } else if gff_list.is_empty() {
                if let Some(gff_path) = &args.gff_file {
                    info!("📎 Pre-parsing single GFF: {}", gff_path);
                    Some(Arc::new(Some(parse_gff_and_build_tree(gff_path)?)))
                } else {
                    None
                }
            } else {
                None
            };
        let cancellation = ParallelCancellation::new(cancel_token);
        let pb = ProgressBar::new(fasta_files.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")?
                .progress_chars("=>-"),
        );
        for fasta_file in fasta_files {
            check_cancelled(cancel_token)?;
            let annotations: Arc<Option<Lapper<usize, Gene>>> =
                if let Some(gff_path) = find_matching_gff(fasta_file, &gff_list) {
                    info!("📎 Matched GFF {} for FASTA {}", gff_path, fasta_file);
                    Arc::new(Some(parse_gff_and_build_tree(gff_path)?))
                } else if let Some(ref shared) = shared_gff {
                    info!("📎 Using shared GFF for FASTA {}", fasta_file);
                    shared.clone()
                } else {
                    info!("📎 No matching GFF for FASTA {} (no annotation)", fasta_file);
                    Arc::new(None)
                };
            let genomes = get_genomes_from_fasta(fasta_file)?;
            debug!("Found {} sequences in {}", genomes.len(), fasta_file);
            let filename = std::path::Path::new(fasta_file)
                .file_name().and_then(|n| n.to_str()).unwrap_or(fasta_file);
            let file_genomes: Vec<(String, String)> = genomes
                .into_iter()
                .map(|(name, seq)| (format!("[{}] {}", filename, name), seq))
                .collect();
            for chunk in file_genomes.chunks(PARALLEL_PROCESS_CHUNK_SIZE) {
                let chunk_results: Vec<Vec<String>> = chunk
                    .par_iter()
                    .map(|(genome_name, genome_seq)| {
                        if cancellation.is_cancelled() { return Vec::new(); }
                        analyze_genome_seq(genome_name, genome_seq, marker_index, annotations.as_ref(), k, ref_seq, ref_seq_rc)
                    })
                    .collect();
                cancellation.check_after()?;
                for genome_lines in chunk_results {
                    for line in genome_lines {
                        emit_line(&line)?;
                    }
                }
            }
            pb.inc(1);
        }
        pb.finish_with_message("Done!");
        cancellation.check_after()?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Summary output
// ---------------------------------------------------------------------------

fn update_summary_from_detail_line(
    lineage_count_map: &mut HashMap<String, HashMap<String, usize>>,
    line: &str,
) {
    let mut fields = line.split('\t');
    let genome_name = match fields.next() {
        Some(value) if !value.is_empty() => value,
        _ => return,
    };
    for _ in 0..6 { if fields.next().is_none() { return; } }
    // Register the genome before inspecting the lineage field: a genome with no
    // marker hits emits a detail row with empty columns, and returning early
    // here would drop it from the summary entirely, leaving the user with fewer
    // summary rows than input genomes and no clue which ones were missing.
    let genome_entry = lineage_count_map.entry(genome_name.to_string()).or_default();
    let lineage_path = match fields.next() {
        Some(value) if !value.is_empty() => value,
        _ => return,
    };
    let mut current_path_part = String::with_capacity(lineage_path.len());
    for (i, component) in lineage_path.split(';').enumerate() {
        if i > 0 { current_path_part.push(';'); }
        current_path_part.push_str(component);
        *genome_entry.entry(current_path_part.clone()).or_default() += 1;
    }
}

fn write_summary(
    summary_out: &mut BufWriter<File>,
    genome: String,
    lineage_counts: &HashMap<String, usize>,
    nested_classification: bool,
) -> AppResult<()> {
    let mut sorted_lineages: Vec<(&String, &usize)> = lineage_counts.iter().collect();
    sorted_lineages.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    let lineage_count_str = sorted_lineages
        .iter()
        .map(|(lin, cnt)| format!("{}:{}", lin, cnt))
        .collect::<Vec<_>>()
        .join(" ");
    let majority_lineage = determine_major_lineage(lineage_counts, nested_classification);
    writeln!(summary_out, "{}\t{}\t{}", genome, lineage_count_str, majority_lineage)?;
    Ok(())
}

fn cleanup_generated_outputs(paths: &[String]) {
    for out_path in paths {
        match std::fs::remove_file(out_path) {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => warn!("Failed to remove partial output {}: {}", out_path, e),
        }
        // Must match the writer's derivation exactly, or cleanup deletes the
        // wrong file and leaves the real partial .xlsx behind.
        let xlsx_path = crate::excel::excel_path_from_tsv(out_path);
        match std::fs::remove_file(&xlsx_path) {
            Ok(_) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => warn!("Failed to remove partial output {}: {}", xlsx_path, e),
        }
    }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

pub fn run(args: Args) -> AppResult<()> {
    configure_thread_pool(args.num_cpu);
    let cancel_token = &args.cancel_token;
    check_cancelled(cancel_token)?;
    check_input_files(&args)?;

    let k = crate::common::validate_kmer_size(args.kmer_size)?;
    debug!("k-mer size set to: {}", k);

    info!("Reading reference and marker positions...");
    let ref_seq = get_ref(&args.ref_fasta)?;
    check_cancelled(cancel_token)?;
    let ref_seq_rc = reverse_complement_sequence(&ref_seq);
    debug!("Reference sequence loaded, length: {}", ref_seq.len());
    let marker_variants = get_positions(&args.tsv_pos)?;
    debug!("Found {} marker entries.", marker_variants.len());

    // A k-mer needs room for the allele between the two flanks; otherwise every
    // marker is skipped and the run silently reports every genome as
    // Unclassified instead of failing.
    if 2 * args.min_flank_bases >= k {
        return Err(AppError::Generic(format!(
            "--min-flank-bases {} leaves no room for an allele in a {}-mer; use a value below {}.",
            args.min_flank_bases,
            k,
            (k + 1) / 2
        )));
    }

    info!("Generating marker k-mers...");
    let marker_index = generate_markerkmer(&marker_variants, &ref_seq, k, args.min_flank_bases);
    if marker_index.len() == 0 {
        return Err(AppError::NotEnoughData(format!(
            "No usable marker k-mers could be built from '{}' (k={}, --min-flank-bases={}). \
             Every marker was skipped, so no genome could be classified.",
            args.tsv_pos, k, args.min_flank_bases
        )));
    }
    debug!("Generated {} unique marker k-mers using {} index storage.", marker_index.len(), marker_index.storage_name());
    check_cancelled(cancel_token)?;

    info!("Processing genomes...");
    let mut generated_outputs: Vec<String> = Vec::new();
    let detailed_file = if args.ofile.ends_with(".tsv") { args.ofile.clone() } else { format!("{}.tsv", args.ofile) };
    generated_outputs.push(detailed_file.clone());

    info!("Writing detailed output to: {}", &detailed_file);
    let mut outfile = match File::create(&detailed_file) {
        Ok(f) => BufWriter::new(f),
        Err(e) => { cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e)); }
    };
    if let Err(e) = writeln!(outfile, "genome\tk-mer\tk-merPOS\tSNPgenome\tSNPreference\tREF\tALT\tlineage\tGene\tGene_Start\tGene_End\tAA_Pos\tAA_Change") {
        drop(outfile);
        cleanup_generated_outputs(&generated_outputs);
        return Err(AppError::Io(e));
    }

    let detail_headers = ["genome", "k-mer", "k-merPOS", "SNPgenome", "SNPreference", "REF", "ALT", "lineage", "Gene", "Gene_Start", "Gene_End", "AA_Pos", "AA_Change"];
    let mut detail_excel_writer = if args.excel {
        match crate::common::ExcelStreamWriter::new(&detailed_file, &detail_headers) {
            Ok(writer) => Some(writer),
            Err(e) => { warn!("⚠️ Could not initialize detailed Excel writer: {}", e); None }
        }
    } else {
        None
    };
    let mut lineage_count_map: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut emitted_lines = 0usize;
    let process_result = process_genomes(&args, &marker_index, k, &ref_seq, &ref_seq_rc, |line| {
        emitted_lines += 1;
        if emitted_lines % 2048 == 0 { check_cancelled(cancel_token)?; }
        write!(outfile, "{}", line).map_err(AppError::Io)?;
        let trimmed = line.trim_end();
        if !trimmed.is_empty() {
            update_summary_from_detail_line(&mut lineage_count_map, trimmed);
        }
        let mut disable_detail_excel = false;
        if let Some(writer) = detail_excel_writer.as_mut() {
            if !trimmed.is_empty() {
                let row: Vec<String> = trimmed.split('\t').map(String::from).collect();
                if let Err(e) = writer.write_row(&row) {
                    warn!("⚠️ Could not append detailed Excel row: {}", e);
                    disable_detail_excel = true;
                }
            }
        }
        if disable_detail_excel { detail_excel_writer = None; }
        Ok(())
    });
    if let Err(e) = process_result {
        drop(outfile);
        cleanup_generated_outputs(&generated_outputs);
        return Err(e);
    }
    if let Err(e) = outfile.flush() {
        drop(outfile); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e));
    }
    drop(outfile);
    if let Err(e) = check_cancelled(cancel_token) {
        cleanup_generated_outputs(&generated_outputs); return Err(e);
    }
    if args.excel {
        if let Some(writer) = detail_excel_writer.take() {
            match writer.finish() {
                Ok(xlsx_path) => info!("📊 Detailed Excel file written to {}", xlsx_path),
                Err(e) => warn!("⚠️ Could not write detailed Excel file: {}", e),
            }
        }
    }

    let base_name = if args.ofile.ends_with(".tsv") { args.ofile.trim_end_matches(".tsv").to_string() } else { args.ofile.clone() };
    let summary_file = format!("{}_summary.tsv", base_name);
    generated_outputs.push(summary_file.clone());
    info!("Writing summary output to: {}", &summary_file);
    let mut summary_out = match File::create(&summary_file) {
        Ok(f) => BufWriter::new(f),
        Err(e) => { cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e)); }
    };
    if let Err(e) = writeln!(summary_out, "genome\tlineage:count\tmajor_lineage") {
        drop(summary_out); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e));
    }

    let mut ordered_lineages: Vec<(String, HashMap<String, usize>)> = lineage_count_map.into_iter().collect();
    ordered_lineages.sort_by(|(a, _), (b, _)| a.cmp(b));
    let summary_headers = ["genome", "lineage:count", "major_lineage"];
    let mut summary_excel_writer = if args.excel {
        match crate::common::ExcelStreamWriter::new(&summary_file, &summary_headers) {
            Ok(writer) => Some(writer),
            Err(e) => { warn!("⚠️ Could not initialize summary Excel writer: {}", e); None }
        }
    } else {
        None
    };

    for (idx, (genome, lineage_map)) in ordered_lineages.into_iter().enumerate() {
        if idx % 128 == 0 {
            if let Err(e) = check_cancelled(cancel_token) {
                drop(summary_out); cleanup_generated_outputs(&generated_outputs); return Err(e);
            }
        }
        // Match the TSV exactly: sorted by count (then name) and space-joined.
        // Iterating the HashMap directly made the Excel column disagree with
        // its own TSV and vary between otherwise identical runs.
        let mut sorted_lineages: Vec<(&String, &usize)> = lineage_map.iter().collect();
        sorted_lineages.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
        let lineage_counts_str: String = sorted_lineages.iter()
            .map(|(lin, count)| format!("{}:{}", lin, count))
            .collect::<Vec<_>>().join(" ");
        let major_lineage = determine_major_lineage(&lineage_map, args.nested_classification);

        let mut disable_summary_excel = false;
        if let Some(writer) = summary_excel_writer.as_mut() {
            let row = vec![genome.clone(), lineage_counts_str.clone(), major_lineage.clone()];
            if let Err(e) = writer.write_row(&row) {
                warn!("⚠️ Could not append summary Excel row: {}", e); disable_summary_excel = true;
            }
        }
        if disable_summary_excel { summary_excel_writer = None; }

        if let Err(e) = write_summary(&mut summary_out, genome, &lineage_map, args.nested_classification) {
            drop(summary_out); cleanup_generated_outputs(&generated_outputs); return Err(e);
        }
    }
    if let Err(e) = summary_out.flush() {
        drop(summary_out); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e));
    }
    drop(summary_out);
    if let Err(e) = check_cancelled(cancel_token) {
        cleanup_generated_outputs(&generated_outputs); return Err(e);
    }
    if args.excel {
        if let Some(writer) = summary_excel_writer.take() {
            match writer.finish() {
                Ok(xlsx_path) => info!("📊 Summary Excel file written to {}", xlsx_path),
                Err(e) => warn!("⚠️ Could not write summary Excel file: {}", e),
            }
        }
    }

    if args.output_masked_fasta {
        info!("Generating masked FASTA files...");
        let mask_ranges = collect_mask_positions(&marker_variants);
        info!("Masking {} genomic regions ({} total bases) across all marker positions.",
            mask_ranges.len(), mask_ranges.iter().map(|(s, e)| e - s).sum::<usize>());

        let fasta_sources: Vec<String> = if let Some(tsv_genomes) = &args.tsv_genomes {
            get_genomepaths(tsv_genomes)?.into_values().map(|(path, _)| path).collect()
        } else if let Some(fasta_genomes) = &args.fasta_genomes {
            vec![fasta_genomes.clone()]
        } else if let Some(fasta_files) = &args.fasta_files {
            fasta_files.clone()
        } else {
            Vec::new()
        };

        let mut used_masked_names: HashSet<String> = HashSet::new();
        for fasta_path in &fasta_sources {
            if let Err(e) = check_cancelled(cancel_token) {
                cleanup_generated_outputs(&generated_outputs);
                return Err(e);
            }
            let stem = Path::new(fasta_path).file_stem().unwrap_or_default().to_string_lossy();
            let out_dir = Path::new(&base_name).parent().unwrap_or(Path::new("."));
            // Inputs from different folders routinely share a file stem (the
            // usual `<sample>/assembly.fasta` layout). Without disambiguation
            // the second masked FASTA silently overwrote the first.
            let mut name = format!("{}_masked.fasta", stem);
            let mut n = 2;
            while !used_masked_names.insert(name.clone()) {
                name = format!("{}_{}_masked.fasta", stem, n);
                n += 1;
            }
            let masked_path = out_dir.join(&name);
            let masked_path_str = masked_path.to_string_lossy().to_string();
            // Register the path before writing so a partially written file is
            // cleaned up if the write fails or the run is cancelled.
            generated_outputs.push(masked_path_str.clone());
            if let Err(e) = write_masked_fasta(fasta_path, &masked_path_str, &mask_ranges) {
                cleanup_generated_outputs(&generated_outputs);
                return Err(e);
            }
        }
    }

    info!("Classification complete.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use crate::lineage::get_final_lineage_call;
    use std::collections::HashMap;

    #[test]
    fn nested_call_is_deterministic_on_ties() {
        let mut counts = HashMap::new();
        counts.insert("L2".to_string(), 2usize);
        counts.insert("L1".to_string(), 2usize);
        assert_eq!(get_final_lineage_call(&counts), "L1");
    }

    #[test]
    fn gff_parses_cds_features_correctly() {
        let gff_path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../TEST_PATHOTYPR/sequence_L4.gff3");
        if !std::path::Path::new(gff_path).exists() {
            eprintln!("Skipping GFF test: file not found at {}", gff_path);
            return;
        }
        let tree = super::parse_gff_and_build_tree(gff_path).expect("Failed to parse GFF");
        let hits: Vec<_> = tree.find(500, 501).collect();
        assert!(!hits.is_empty(), "Expected CDS hit at position 500 (within dnaA)");
        let gene = &hits[0].val;
        eprintln!("Gene ID at pos 500: {}", gene.id);
        assert!(!gene.id.is_empty());
        assert_ne!(gene.id, "Unknown");
    }

    #[test]
    fn nested_call_uses_lexical_tiebreak_for_same_depth() {
        let mut counts = HashMap::new();
        counts.insert("L1".to_string(), 3usize);
        counts.insert("L1;B".to_string(), 2usize);
        counts.insert("L1;A".to_string(), 2usize);
        assert_eq!(get_final_lineage_call(&counts), "L1;A");
    }

    #[test]
    fn markers_are_found_on_reverse_oriented_contigs() {
        use crate::classify::markers::{generate_markerkmer, MarkerVariant};
        use crate::classify::{analyze_genome_seq, reverse_complement_sequence};

        // 60 bp reference; the marker sits at 1-based position 30.
        let ref_seq = concat!(
            "ACGTTGCAAG", "GCTTAACCGG", "ATCGATTCAG", "CTAGCCATGG", "TACGTTAACG", "GCATTGCAGT",
        );
        assert_eq!(ref_seq.len(), 60);
        assert_eq!(&ref_seq[29..30], "G", "reference allele at the marker position");

        let markers = vec![MarkerVariant {
            pos: 30,
            ref_allele: "G".to_string(),
            alt_allele: "A".to_string(),
            lineage: "L9.9".to_string(),
            gene: None,
            mutation: None,
        }];
        let k = 11;
        let index = generate_markerkmer(&markers, ref_seq, k, 5);
        let ref_rc = reverse_complement_sequence(ref_seq);

        // A genome carrying the ALT allele, in the same orientation as the reference.
        let mut forward = ref_seq.to_string();
        forward.replace_range(29..30, "A");
        let hit_fwd = analyze_genome_seq("fwd", &forward, &index, &None, k, ref_seq, &ref_rc);
        assert!(
            hit_fwd.iter().any(|l| l.contains("L9.9")),
            "marker must be found on a forward-oriented contig"
        );

        // The very same contig stored in the opposite orientation, which is
        // arbitrary in a draft assembly, must yield the same call.
        let reverse = reverse_complement_sequence(&forward);
        let hit_rev = analyze_genome_seq("rev", &reverse, &index, &None, k, ref_seq, &ref_rc);
        assert!(
            hit_rev.iter().any(|l| l.contains("L9.9")),
            "marker must also be found when the contig is reverse-oriented"
        );
    }

    #[test]
    fn a_marker_present_on_several_contigs_is_counted_once() {
        use crate::classify::markers::{generate_markerkmer, MarkerVariant};
        use crate::classify::{analyze_genome_seq, reverse_complement_sequence};

        let ref_seq = concat!(
            "ACGTTGCAAG", "GCTTAACCGG", "ATCGATTCAG", "CTAGCCATGG", "TACGTTAACG", "GCATTGCAGT",
        );
        let markers = vec![MarkerVariant {
            pos: 30,
            ref_allele: "G".to_string(),
            alt_allele: "A".to_string(),
            lineage: "L9.9".to_string(),
            gene: None,
            mutation: None,
        }];
        let k = 11;
        let index = generate_markerkmer(&markers, ref_seq, k, 5);
        let ref_rc = reverse_complement_sequence(ref_seq);

        let mut forward = ref_seq.to_string();
        forward.replace_range(29..30, "A");
        // Same sequence twice over: the marker is one observation, not two.
        let doubled = format!("{}{}", forward, forward);
        let lines = analyze_genome_seq("dup", &doubled, &index, &None, k, ref_seq, &ref_rc);
        let hits = lines.iter().filter(|l| l.contains("L9.9")).count();
        assert_eq!(hits, 1, "a repeated marker must not inflate the lineage tally");
    }

}
