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
    find_markers, generate_markerkmer, get_positions, scan_anchors, scan_markers, AnchorEvidence,
    AnchorIndex, AnchorRead, AnchorScan, MarkerIndex, MarkerKmerEntry, MarkerMatch, MarkerScan,
    MarkerVariant,
};
pub use masking::{collect_mask_positions, mask_sequence, write_masked_fasta};

use annotation::{translate_variant_info, GenomeContext};

/// Bases kept either side of a variant so its codon can be read from the sample.
const CODON_MARGIN: usize = 2;

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
    // `--input-files` is a full input source too, so it must satisfy this
    // requirement as well; listing only `--input` made `--input-files` on its
    // own fail with "the following required arguments were not provided".
    #[arg(
        short = 'l',
        long = "input-list",
        required_unless_present_any = ["fasta_genomes", "fasta_files"]
    )]
    pub tsv_genomes: Option<String>,

    /// Single FASTA file to analyze.
    #[arg(short = 'i', long = "input", required_unless_present_any = ["tsv_genomes", "fasta_files"])]
    pub fasta_genomes: Option<String>,

    /// Multiple FASTA files to analyze (for GUI batch mode).
    #[arg(
        long = "input-files",
        num_args = 1..,
        required_unless_present_any = ["tsv_genomes", "fasta_genomes"]
    )]
    pub fasta_files: Option<Vec<String>>,

    /// Optional GFF file for annotation when using --input.
    #[arg(long = "gff", requires = "fasta_genomes")]
    pub gff_file: Option<String>,

    /// Multiple GFF files for annotation when using --input-files (matched by filename).
    #[arg(long = "gff-files", num_args = 1..)]
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

/// Columns of the detailed output, one row per marker whose alternative allele
/// was found.
const DETAIL_HEADERS: [&str; 15] = [
    "genome", "k-mer", "k-merPOS", "SNPgenome", "SNPreference", "REF", "ALT", "lineage",
    "Gene", "Gene_Start", "Gene_End", "AA_Pos", "AA_Change", "alleles_found", "evidence",
];

/// Columns of the undetermined output, one row per marker that could be read
/// neither way.
const UNDETERMINED_HEADERS: [&str; 8] = [
    "genome", "SNPreference", "REF", "ALT", "lineage", "Gene", "mutation", "reason",
];

/// Output rows for one genome.
///
/// A marker ends up in one of four states. Its alternative k-mer found:
/// a `detail` row, with `alleles_found` saying `ALT`, or `ALT,REF` when the
/// reference k-mer was found as well (two copies of the region that disagree,
/// or a k-mer that also occurs elsewhere in the genome). Its reference k-mer
/// found alone: no row, because the genome carries the reference allele.
/// Neither found: the anchors are consulted, and the marker is reported
/// present, absent or, when they cannot read it either, as `undetermined`.
#[derive(Clone, Debug, Default)]
pub struct GenomeLines {
    pub detail: Vec<String>,
    pub undetermined: Vec<String>,
}

fn analyze_genome(
    genome_name: &str,
    fasta_path: &str,
    gff_path: &Option<String>,
    shared_annotations: Option<&Lapper<usize, Gene>>,
    marker_index: &MarkerIndex,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> AppResult<GenomeLines> {
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

/// Collects the output rows for one genome, scanning every contig on **both
/// strands** and reporting each marker only once.
///
/// Contig orientation in a draft assembly is arbitrary, so a marker sitting on
/// a reverse-oriented contig only matches the reverse complement; scanning the
/// forward strand alone silently loses it (the FASTQ path already indexes both
/// orientations). Matches are deduplicated per marker, so a marker seen on
/// several contigs, or on both strands of one contig, counts as a single
/// observation rather than inflating the lineage tally.
///
/// A marker's state is only known once every contig has been scanned, since
/// its two alleles can turn up on different contigs, so rows are written at
/// the end.
fn collect_marker_lines<'a>(
    genome_name: &str,
    contigs: impl Iterator<Item = &'a str>,
    marker_index: &MarkerIndex,
    annotations: &Option<&Lapper<usize, Gene>>,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> GenomeLines {
    // Collected because the anchors need a second look at the same contigs.
    let contigs: Vec<String> = contigs.map(str::to_string).collect();
    let mut alt_found = vec![false; marker_index.len()];
    let mut reference_found = vec![false; marker_index.site_count()];
    // A site is read when any k-mer of it turned up: its reference allele, or
    // the alternative allele of any marker there. Reading one allele of a
    // site rules out the others at the same position.
    let mut site_read = vec![false; marker_index.site_count()];
    // How many times each site's reference window turned up, summed over every
    // contig and strand. A window the reference carries twice is proved absent
    // by finding it twice. Overlapping contigs of a draft assembly can inflate
    // this, which only takes the marker back to being reported absent.
    let mut site_hit_counts = vec![0u32; marker_index.site_count()];
    let mut observed: HashMap<usize, (MarkerMatch, String, GenomeContext)> = HashMap::new();

    for contig_seq in &contigs {
        let contig_rc = reverse_complement_sequence(contig_seq);
        for seq in [contig_seq.as_str(), contig_rc.as_str()] {
            let scan = scan_markers(seq, marker_index, k);
            for site in scan.reference_sites {
                reference_found[site] = true;
                site_read[site] = true;
            }
            for (site, hits) in scan.site_hits.iter().enumerate() {
                site_hit_counts[site] = site_hit_counts[site].saturating_add(*hits);
            }
            for m in scan.matches {
                if alt_found[m.marker_id] {
                    continue;
                }
                alt_found[m.marker_id] = true;
                if let Some(site) = marker_index.site_of(m.marker_id) {
                    site_read[site] = true;
                }
                // `seq` is the strand that produced the match, so the positions
                // and the ALT bases sliced out of it stay self-consistent.
                let variant_start = m.genome_position + m.left_flank_len;
                let variant_end = (variant_start + m.alt_allele_len).min(seq.len());
                let alt_bases = seq[variant_start..variant_end].to_string();
                // Two bases either side put the whole codon in view, and the
                // k-mer's own flanks guarantee they are inside the match.
                let from = variant_start.saturating_sub(CODON_MARGIN);
                let to = (variant_end + CODON_MARGIN).min(seq.len());
                let context = GenomeContext {
                    bases: seq[from..to].as_bytes().to_vec(),
                    before: variant_start - from,
                };
                observed.insert(m.marker_id, (m, alt_bases, context));
            }
        }
    }

    // The k-mers alone cannot say *where* they matched, so a copy of either
    // allele anywhere in the genome looks like the marker's own locus. The
    // anchors can: each occurs once or it is not used. So wherever an anchor
    // can read this locus, it decides, and the k-mers are the fallback.
    //
    // This matters in both directions. A paralogous copy of the ALT k-mer made
    // two markers of the published resistance panel "present" in every genome,
    // the reference included. A paralogous copy of the REF k-mer did the
    // opposite: it made a genome carrying rpoB S450L read as susceptible.
    let anchor_scan = (!marker_index.anchors().is_empty()).then(|| {
        let mut scan = AnchorScan::new(marker_index.anchors());
        for contig_seq in &contigs {
            for seq in [contig_seq.as_str(), reverse_complement_sequence(contig_seq).as_str()] {
                scan_anchors(seq, marker_index.anchors(), &mut scan);
            }
        }
        scan
    });

    let mut lines = GenomeLines::default();
    for id in 0..marker_index.len() {
        let entry = &marker_index.markers()[id];
        let kmer_alt = alt_found[id];
        let kmer_reference = marker_index.site_of(id).is_some_and(|site| reference_found[site]);
        let alleles_found = if kmer_alt && kmer_reference { "ALT,REF" } else { "ALT" };
        let read = anchor_scan
            .as_ref()
            .map(|scan| scan.allele_of(marker_index.anchors(), id, entry.alt_allele_len))
            .unwrap_or(AnchorRead::Unusable);

        let reason = match read {
            // The anchor read this locus and it carries the variant.
            AnchorRead::Allele(allele, context, evidence) if allele == entry.alt_allele => {
                match observed.get(&id).filter(|_| !kmer_reference) {
                    // The k-mer matched too and nothing contradicts it, so keep
                    // its row: it carries the genome coordinates of the match.
                    Some((m, alt_bases, ctx)) => lines.detail.push(format_marker_match(
                        genome_name, m, annotations, alt_bases, Some(ctx), alleles_found, "kmer", ref_seq, ref_seq_rc,
                    )),
                    None => lines.detail.push(format_anchor_match(
                        genome_name, entry, annotations, Some(&context), alleles_found, evidence, ref_seq, ref_seq_rc,
                    )),
                }
                continue;
            }
            // The locus carries the reference allele. This overrules a k-mer
            // hit, which can only have come from somewhere else.
            // Anything else: the locus was read and does not carry this
            // marker's allele, so the marker is absent. What it carries
            // instead is another marker's business, and finding variants
            // nobody catalogued is a variant caller's job, not a panel's.
            AnchorRead::Allele(..) => continue,
            AnchorRead::Conflict => "anchors_disagree",
            // No anchor could read it, so fall back to what the k-mers saw.
            AnchorRead::Unusable => {
                if entry.alt_in_reference && kmer_alt {
                    // The reference itself carries this alternative k-mer, so
                    // finding it says nothing: every genome has it. Only an
                    // anchor could have told the real locus from the copy, and
                    // none could be read.
                    "alt_kmer_in_reference"
                } else if kmer_alt && !kmer_reference {
                    match observed.get(&id) {
                        Some((m, alt_bases, ctx)) => {
                            lines.detail.push(format_marker_match(
                                genome_name, m, annotations, alt_bases, Some(ctx), "ALT", "kmer", ref_seq, ref_seq_rc,
                            ));
                            continue;
                        }
                        None => "not_found",
                    }
                } else if kmer_alt {
                    // Both alleles turned up and nothing can say which copy
                    // this locus is. Reporting it present would be a guess.
                    "ambiguous_both_alleles"
                } else if marker_index
                    .site_of(id)
                    .is_some_and(|site| site_read[site] && site_hit_counts[site] < entry.ref_copies)
                {
                    // The reference holds this window more than once and the
                    // genome did not show every copy, so the hit may be the
                    // other copy and cannot rule the marker out. Showing them
                    // all does rule it out: every copy was read as reference.
                    "ref_kmer_repeated"
                } else if marker_index.site_of(id).is_some_and(|site| site_read[site]) {
                    continue; // read as the reference allele, so absent
                } else if marker_index.site_of(id).is_some() {
                    "not_found"
                } else {
                    "no_reference_kmer"
                }
            }
        };
        lines.undetermined.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            genome_name,
            entry.ref_position,
            entry.ref_allele,
            entry.alt_allele,
            entry.lineage,
            entry.gene.as_deref().unwrap_or(""),
            entry.mutation.as_deref().unwrap_or(""),
            reason
        ));
    }

    if lines.detail.is_empty() {
        lines.detail.push(format!("{}{}\n", genome_name, "\t".repeat(DETAIL_HEADERS.len() - 1)));
    }
    lines
}

pub fn analyze_genome_seq(
    genome_name: &str,
    genome_seq: &str,
    marker_index: &MarkerIndex,
    annotations: &Option<Lapper<usize, Gene>>,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> GenomeLines {
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

/// Gene and amino-acid columns for one variant, from the GFF when there is one
/// and from the marker's own fields otherwise.
fn annotate(
    annotations: &Option<&Lapper<usize, Gene>>,
    ref_position: usize,
    ref_len: usize,
    alt_bases: &str,
    context: Option<&GenomeContext>,
    marker_gene: Option<&str>,
    marker_mutation: Option<&str>,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> (String, String, String, String, String) {
    let ref_pos_0based = ref_position.saturating_sub(1);
    let (gff_gene_id, gene_start, gene_end, gff_aa_pos, gff_aa_change) = if let Some(tree) = annotations {
        if let Some(gene_interval) = tree.find(ref_pos_0based, ref_pos_0based + ref_len).next() {
            let (gene_id, aa_pos, aa_change) = translate_variant_info(
                gene_interval, ref_pos_0based, alt_bases, ref_len, context, ref_seq, ref_seq_rc,
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

    let gene_id = gff_gene_id.or_else(|| marker_gene.map(str::to_string));
    let aa_change = gff_aa_change.or_else(|| marker_mutation.map(str::to_string));
    (
        gene_id.unwrap_or_default(),
        gene_start.unwrap_or_default(),
        gene_end.unwrap_or_default(),
        gff_aa_pos.unwrap_or_default(),
        aa_change.unwrap_or_default(),
    )
}

fn format_marker_match(
    genome_name: &str,
    m: &MarkerMatch,
    annotations: &Option<&Lapper<usize, Gene>>,
    alt_bases: &str,
    context: Option<&GenomeContext>,
    alleles_found: &str,
    evidence: &str,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> String {
    let variant_start = m.genome_position + m.left_flank_len;
    let (gene, gene_start, gene_end, aa_pos, aa_change) = annotate(
        annotations,
        m.ref_position,
        m.ref_allele_len,
        alt_bases,
        context,
        m.gene.as_deref(),
        m.mutation.as_deref(),
        ref_seq,
        ref_seq_rc,
    );
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        genome_name,
        m.kmer,
        m.genome_position + 1,
        variant_start + 1,
        m.ref_position,
        m.ref_allele,
        m.alt_allele,
        m.lineage,
        gene,
        gene_start,
        gene_end,
        aa_pos,
        aa_change,
        alleles_found,
        evidence
    )
}

/// A marker the anchors read as present. The k-mer never matched, so there are
/// no genome coordinates to report: those columns stay empty rather than
/// carrying a position no match established.
fn format_anchor_match(
    genome_name: &str,
    entry: &MarkerKmerEntry,
    annotations: &Option<&Lapper<usize, Gene>>,
    context: Option<&GenomeContext>,
    alleles_found: &str,
    evidence: AnchorEvidence,
    ref_seq: &str,
    ref_seq_rc: &str,
) -> String {
    let (gene, gene_start, gene_end, aa_pos, aa_change) = annotate(
        annotations,
        entry.ref_position,
        entry.ref_allele_len,
        &entry.alt_allele,
        context,
        entry.gene.as_deref(),
        entry.mutation.as_deref(),
        ref_seq,
        ref_seq_rc,
    );
    format!(
        "{}\t\t\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        genome_name,
        entry.ref_position,
        entry.ref_allele,
        entry.alt_allele,
        entry.lineage,
        gene,
        gene_start,
        gene_end,
        aa_pos,
        aa_change,
        alleles_found,
        evidence.as_str()
    )
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

/// One output row, tagged with the file it belongs to.
enum OutputRow<'a> {
    Detail(&'a str),
    Undetermined(&'a str),
}

fn emit_genome<F>(lines: &GenomeLines, emit_line: &mut F) -> AppResult<()>
where
    F: FnMut(OutputRow<'_>) -> AppResult<()>,
{
    for line in &lines.detail {
        emit_line(OutputRow::Detail(line))?;
    }
    for line in &lines.undetermined {
        emit_line(OutputRow::Undetermined(line))?;
    }
    Ok(())
}

fn process_genomes<F>(
    args: &Args,
    marker_index: &MarkerIndex,
    k: usize,
    ref_seq: &str,
    ref_seq_rc: &str,
    mut emit_line: F,
) -> AppResult<()>
where
    F: FnMut(OutputRow<'_>) -> AppResult<()>,
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
            let res_nested: AppResult<Vec<GenomeLines>> = chunk
                .par_iter()
                .map(|(genome_name, (fasta_path, gff_path))| {
                    if cancellation.is_cancelled() { return Ok(GenomeLines::default()); }
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
                emit_genome(&genome_lines, &mut emit_line)?;
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
            let chunk_results: Vec<GenomeLines> = batch
                .par_iter()
                .map(|(genome_name, genome_seq)| {
                    if cancellation.is_cancelled() { return GenomeLines::default(); }
                    let lines = analyze_genome_seq(genome_name, genome_seq, marker_index, &annotations, k, ref_seq, ref_seq_rc);
                    pb.inc(1);
                    lines
                })
                .collect();
            cancellation.check_after()?;
            for genome_lines in chunk_results {
                emit_genome(&genome_lines, &mut emit_line)?;
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
                let chunk_results: Vec<GenomeLines> = chunk
                    .par_iter()
                    .map(|(genome_name, genome_seq)| {
                        if cancellation.is_cancelled() { return GenomeLines::default(); }
                        analyze_genome_seq(genome_name, genome_seq, marker_index, annotations.as_ref(), k, ref_seq, ref_seq_rc)
                    })
                    .collect();
                cancellation.check_after()?;
                for genome_lines in chunk_results {
                    emit_genome(&genome_lines, &mut emit_line)?;
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

/// The strongest verdict per drug, per genome, for a resistance panel.
type ResistanceMap = HashMap<String, HashMap<String, crate::resistance::Verdict>>;

/// Folds one detail row into the per-drug verdicts of its genome.
fn update_resistance_from_detail_line(resistance_map: &mut ResistanceMap, line: &str) {
    let mut fields = line.split('\t');
    let genome_name = match fields.next() {
        Some(value) if !value.is_empty() => value,
        _ => return,
    };
    // Register the genome even when it matched nothing, so it still gets a
    // call rather than dropping out of the summary.
    let entry = resistance_map.entry(genome_name.to_string()).or_default();
    if let Some(path) = fields.nth(6).filter(|p| !p.is_empty()) {
        crate::resistance::fold_marker(entry, path);
    }
}

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

/// Tallies an undetermined row under the first level of its marker's path:
/// the drug for a resistance panel, the main lineage for a lineage panel.
fn update_undetermined_from_line(
    undetermined_count_map: &mut HashMap<String, HashMap<String, usize>>,
    line: &str,
) {
    let mut fields = line.split('\t');
    let genome_name = match fields.next() {
        Some(value) if !value.is_empty() => value,
        _ => return,
    };
    let label = fields
        .nth(3)
        .and_then(|path| path.split(';').next())
        .unwrap_or("");
    *undetermined_count_map
        .entry(genome_name.to_string())
        .or_default()
        .entry(label.to_string())
        .or_default() += 1;
}

/// `label:count` pairs, most frequent first and then by name, space-joined.
fn format_label_counts(counts: &HashMap<String, usize>) -> String {
    let mut sorted: Vec<(&String, &usize)> = counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    sorted
        .iter()
        .map(|(label, count)| format!("{}:{}", label, count))
        .collect::<Vec<_>>()
        .join(" ")
}

/// The two undetermined summary columns: the total and the per-label tally.
fn undetermined_columns(counts: Option<&HashMap<String, usize>>) -> (usize, String) {
    match counts {
        Some(counts) => (counts.values().sum(), format_label_counts(counts)),
        None => (0, String::new()),
    }
}

/// The call for one genome: a lineage, or the resistance verdict when the
/// panel is the WHO catalogue, where tallying paths answers the wrong question.
fn genome_call(
    lineage_counts: &HashMap<String, usize>,
    resistance_counts: Option<&HashMap<String, crate::resistance::Verdict>>,
    nested_classification: bool,
) -> String {
    match resistance_counts {
        Some(per_drug) => crate::resistance::call(per_drug),
        None => determine_major_lineage(lineage_counts, nested_classification),
    }
}

fn write_summary(
    summary_out: &mut BufWriter<File>,
    genome: String,
    lineage_counts: &HashMap<String, usize>,
    undetermined_counts: Option<&HashMap<String, usize>>,
    resistance_counts: Option<&HashMap<String, crate::resistance::Verdict>>,
    nested_classification: bool,
) -> AppResult<()> {
    let majority_lineage = genome_call(lineage_counts, resistance_counts, nested_classification);
    let (undetermined_total, undetermined_by_label) = undetermined_columns(undetermined_counts);
    writeln!(
        summary_out,
        "{}\t{}\t{}\t{}\t{}",
        genome,
        format_label_counts(lineage_counts),
        majority_lineage,
        undetermined_total,
        undetermined_by_label
    )?;
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

    // A WHO resistance catalogue is not a hierarchy, so tallying its paths the
    // lineage way answers "which drug has the most markers", and most markers
    // of the catalogue are graded as not associated with resistance. Recognise
    // it once here and answer it on its own terms.
    let resistance_panel =
        crate::resistance::is_resistance_panel(marker_variants.iter().map(|m| m.lineage.as_str()));
    if resistance_panel {
        info!("The panel is a drug-resistance catalogue; the summary reports the drugs with a WHO grade 1 or 2 marker.");
    }

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
    if marker_index.is_empty() {
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
    let base_name = if args.ofile.ends_with(".tsv") { args.ofile.trim_end_matches(".tsv").to_string() } else { args.ofile.clone() };
    generated_outputs.push(detailed_file.clone());

    info!("Writing detailed output to: {}", &detailed_file);
    let mut outfile = match File::create(&detailed_file) {
        Ok(f) => BufWriter::new(f),
        Err(e) => { cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e)); }
    };
    if let Err(e) = writeln!(outfile, "{}", DETAIL_HEADERS.join("\t")) {
        drop(outfile);
        cleanup_generated_outputs(&generated_outputs);
        return Err(AppError::Io(e));
    }

    // Markers that could be read neither way. Without this file a marker the
    // genome hides behind a nearby change looks exactly like one it lacks.
    let undetermined_file = format!("{}_undetermined.tsv", base_name);
    generated_outputs.push(undetermined_file.clone());
    info!("Writing undetermined markers to: {}", &undetermined_file);
    let mut undetermined_out = match File::create(&undetermined_file) {
        Ok(f) => BufWriter::new(f),
        Err(e) => { drop(outfile); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e)); }
    };
    if let Err(e) = writeln!(undetermined_out, "{}", UNDETERMINED_HEADERS.join("\t")) {
        drop(outfile);
        drop(undetermined_out);
        cleanup_generated_outputs(&generated_outputs);
        return Err(AppError::Io(e));
    }

    let mut detail_excel_writer = if args.excel {
        match crate::common::ExcelStreamWriter::new(&detailed_file, &DETAIL_HEADERS) {
            Ok(writer) => Some(writer),
            Err(e) => { warn!("⚠️ Could not initialize detailed Excel writer: {}", e); None }
        }
    } else {
        None
    };
    let mut undetermined_excel_writer = if args.excel {
        match crate::common::ExcelStreamWriter::new(&undetermined_file, &UNDETERMINED_HEADERS) {
            Ok(writer) => Some(writer),
            Err(e) => { warn!("⚠️ Could not initialize undetermined Excel writer: {}", e); None }
        }
    } else {
        None
    };
    let mut lineage_count_map: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut undetermined_count_map: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let mut resistance_map: ResistanceMap = HashMap::new();
    let mut emitted_lines = 0usize;
    let process_result = process_genomes(&args, &marker_index, k, &ref_seq, &ref_seq_rc, |row| {
        emitted_lines += 1;
        if emitted_lines % 2048 == 0 { check_cancelled(cancel_token)?; }
        match row {
            OutputRow::Detail(line) => {
                write!(outfile, "{}", line).map_err(AppError::Io)?;
                // Trim the newline only, not the trailing tabs: a genome with no marker
                // hits emits a row of empty columns, and trimming those away leaves a
                // bare name that the summary parser rejects as malformed, which is how
                // such genomes used to vanish from the summary entirely.
                let trimmed = line.trim_end_matches(['\n', '\r']);
                if !trimmed.trim().is_empty() {
                    update_summary_from_detail_line(&mut lineage_count_map, trimmed);
                    if resistance_panel {
                        update_resistance_from_detail_line(&mut resistance_map, trimmed);
                    }
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
            }
            OutputRow::Undetermined(line) => {
                write!(undetermined_out, "{}", line).map_err(AppError::Io)?;
                let trimmed = line.trim_end_matches(['\n', '\r']);
                update_undetermined_from_line(&mut undetermined_count_map, trimmed);
                let mut disable_undetermined_excel = false;
                if let Some(writer) = undetermined_excel_writer.as_mut() {
                    let row: Vec<String> = trimmed.split('\t').map(String::from).collect();
                    if let Err(e) = writer.write_row(&row) {
                        warn!("⚠️ Could not append undetermined Excel row: {}", e);
                        disable_undetermined_excel = true;
                    }
                }
                if disable_undetermined_excel { undetermined_excel_writer = None; }
            }
        }
        Ok(())
    });
    if let Err(e) = process_result {
        drop(outfile);
        drop(undetermined_out);
        cleanup_generated_outputs(&generated_outputs);
        return Err(e);
    }
    if let Err(e) = outfile.flush() {
        drop(outfile); drop(undetermined_out); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e));
    }
    drop(outfile);
    if let Err(e) = undetermined_out.flush() {
        drop(undetermined_out); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e));
    }
    drop(undetermined_out);
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
        if let Some(writer) = undetermined_excel_writer.take() {
            match writer.finish() {
                Ok(xlsx_path) => info!("📊 Undetermined Excel file written to {}", xlsx_path),
                Err(e) => warn!("⚠️ Could not write undetermined Excel file: {}", e),
            }
        }
    }

    let undetermined_total: usize = undetermined_count_map
        .values()
        .map(|counts| counts.values().sum::<usize>())
        .sum();
    if undetermined_total > 0 {
        info!(
            "{} markers could be read neither way, in {} of {} genomes: some other change inside the \
             k-mer window, a missing region or a contig break hides them. They are listed in {}.",
            undetermined_total,
            undetermined_count_map.len(),
            lineage_count_map.len(),
            &undetermined_file
        );
    } else {
        info!("Every marker was read in every genome, as present or as reference.");
    }

    let summary_file = format!("{}_summary.tsv", base_name);
    generated_outputs.push(summary_file.clone());
    info!("Writing summary output to: {}", &summary_file);
    let mut summary_out = match File::create(&summary_file) {
        Ok(f) => BufWriter::new(f),
        Err(e) => { cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e)); }
    };
    let summary_headers = ["genome", "lineage:count", "major_lineage", "undetermined", "undetermined:count"];
    if let Err(e) = writeln!(summary_out, "{}", summary_headers.join("\t")) {
        drop(summary_out); cleanup_generated_outputs(&generated_outputs); return Err(AppError::Io(e));
    }

    let mut ordered_lineages: Vec<(String, HashMap<String, usize>)> = lineage_count_map.into_iter().collect();
    ordered_lineages.sort_by(|(a, _), (b, _)| a.cmp(b));
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
        let lineage_counts_str = format_label_counts(&lineage_map);
        let resistance_counts = resistance_panel.then(|| resistance_map.get(&genome)).flatten();
        let major_lineage = genome_call(&lineage_map, resistance_counts, args.nested_classification);
        let undetermined_counts = undetermined_count_map.get(&genome);
        let (undetermined_total, undetermined_by_label) = undetermined_columns(undetermined_counts);

        let mut disable_summary_excel = false;
        if let Some(writer) = summary_excel_writer.as_mut() {
            let row = vec![
                genome.clone(),
                lineage_counts_str,
                major_lineage,
                undetermined_total.to_string(),
                undetermined_by_label,
            ];
            if let Err(e) = writer.write_row(&row) {
                warn!("⚠️ Could not append summary Excel row: {}", e); disable_summary_excel = true;
            }
        }
        if disable_summary_excel { summary_excel_writer = None; }

        if let Err(e) = write_summary(&mut summary_out, genome.clone(), &lineage_map, undetermined_counts, resistance_counts, args.nested_classification) {
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
            match write_masked_fasta(fasta_path, &masked_path_str, &mask_ranges, ref_seq.len()) {
                Ok(true) => {}
                // Not in reference coordinates: nothing was written (the reason
                // is logged), so drop the reservation and free the name again.
                Ok(false) => {
                    generated_outputs.pop();
                    used_masked_names.remove(&name);
                }
                Err(e) => {
                    cleanup_generated_outputs(&generated_outputs);
                    return Err(e);
                }
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
            hit_fwd.detail.iter().any(|l| l.contains("L9.9")),
            "marker must be found on a forward-oriented contig"
        );

        // The very same contig stored in the opposite orientation, which is
        // arbitrary in a draft assembly, must yield the same call.
        let reverse = reverse_complement_sequence(&forward);
        let hit_rev = analyze_genome_seq("rev", &reverse, &index, &None, k, ref_seq, &ref_rc);
        assert!(
            hit_rev.detail.iter().any(|l| l.contains("L9.9")),
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
        let hits = lines.detail.iter().filter(|l| l.contains("L9.9")).count();
        assert_eq!(hits, 1, "a repeated marker must not inflate the lineage tally");
    }


    #[test]
    fn a_genome_with_no_markers_still_reaches_the_summary() {
        use std::collections::HashMap;
        // The exact row analyze_genome emits when nothing matched: the name
        // followed by empty columns. It must register the genome, or the sample
        // silently disappears from the summary.
        let line = format!("{}{}", "sampleX", "\t".repeat(super::DETAIL_HEADERS.len() - 1));
        let mut map: HashMap<String, HashMap<String, usize>> = HashMap::new();
        super::update_summary_from_detail_line(&mut map, &line);
        assert!(map.contains_key("sampleX"), "the genome must appear in the summary");
        assert!(map["sampleX"].is_empty(), "with no lineage counts");
        assert_eq!(
            crate::lineage::determine_major_lineage(&map["sampleX"], false),
            "Unclassified"
        );
    }

    // -----------------------------------------------------------------------
    // The reference k-mer: present, absent, both, and undetermined
    // -----------------------------------------------------------------------

    use crate::classify::markers::{generate_markerkmer, MarkerVariant};
    use crate::classify::{analyze_genome_seq, reverse_complement_sequence, GenomeLines};

    /// 60 bp reference shared by the tests below. With k = 11 the window of a
    /// single-base marker at position 30 covers positions 25 to 35.
    const REF60: &str = concat!(
        "ACGTTGCAAG", "GCTTAACCGG", "ATCGATTCAG", "CTAGCCATGG", "TACGTTAACG", "GCATTGCAGT",
    );

    fn snp_marker(alt: &str) -> MarkerVariant {
        MarkerVariant {
            pos: 30,
            ref_allele: "G".to_string(),
            alt_allele: alt.to_string(),
            lineage: format!("L9.9;{}", alt),
            gene: None,
            mutation: None,
        }
    }

    /// The reference with single-base edits at 1-based positions.
    fn edited(base: &str, edits: &[(usize, &str)]) -> String {
        let mut genome = base.to_string();
        for (pos, bases) in edits {
            genome.replace_range(pos - 1..*pos, bases);
        }
        genome
    }

    fn classify_one(genome: &str, markers: &[MarkerVariant]) -> GenomeLines {
        let k = 11;
        let index = generate_markerkmer(markers, REF60, k, 5);
        let ref_rc = reverse_complement_sequence(REF60);
        analyze_genome_seq("g", genome, &index, &None, k, REF60, &ref_rc)
    }

    /// Column of a row, by its name in the header.
    fn col<'a>(headers: &[&str], row: &'a str, name: &str) -> &'a str {
        let i = headers.iter().position(|h| *h == name).expect("column exists");
        row.trim_end_matches('\n').split('\t').nth(i).unwrap_or("")
    }

    fn alleles_found(row: &str) -> &str {
        col(&super::DETAIL_HEADERS, row, "alleles_found")
    }

    fn evidence(row: &str) -> &str {
        col(&super::DETAIL_HEADERS, row, "evidence")
    }

    fn reason(row: &str) -> &str {
        col(&super::UNDETERMINED_HEADERS, row, "reason")
    }

    #[test]
    fn a_genome_like_the_reference_confirms_the_marker_absent() {
        assert_eq!(&REF60[29..30], "G");
        let lines = classify_one(REF60, &[snp_marker("A")]);
        assert!(!lines.detail.iter().any(|l| l.contains("L9.9")));
        assert!(
            lines.undetermined.is_empty(),
            "the reference k-mer was found, so the marker is absent, not undetermined"
        );
    }

    #[test]
    fn a_genome_with_the_variant_reports_only_the_alternative_allele() {
        let genome = edited(REF60, &[(30, "A")]);
        let lines = classify_one(&genome, &[snp_marker("A")]);
        let row = lines.detail.iter().find(|l| l.contains("L9.9")).expect("marker found");
        assert_eq!(alleles_found(row), "ALT");
        assert!(lines.undetermined.is_empty());
    }



    #[test]
    fn a_change_outside_the_window_changes_nothing() {
        assert_eq!(&REF60[39..40], "G");
        let lines = classify_one(&edited(REF60, &[(40, "C")]), &[snp_marker("A")]);
        assert!(lines.undetermined.is_empty());
        let lines = classify_one(&edited(REF60, &[(30, "A"), (40, "C")]), &[snp_marker("A")]);
        assert!(lines.detail.iter().any(|l| l.contains("L9.9")));
        assert!(lines.undetermined.is_empty());
    }

    #[test]
    fn two_copies_that_disagree_are_not_a_variant_call() {
        // One contig with each allele. In a haploid genome that is a question,
        // not a call: which copy is the marker's own locus? With a 60 bp
        // reference both anchors occur twice, so nothing can answer it.
        let k = 11;
        let index = generate_markerkmer(&[snp_marker("A")], REF60, k, 5);
        let ref_rc = reverse_complement_sequence(REF60);
        let variant = edited(REF60, &[(30, "A")]);
        let lines = super::collect_marker_lines(
            "g",
            [REF60, variant.as_str()].into_iter(),
            &index,
            &None,
            k,
            REF60,
            &ref_rc,
        );
        assert!(
            !lines.detail.iter().any(|l| l.contains("L9.9")),
            "an ambiguous marker must not be reported present"
        );
        assert_eq!(reason(&lines.undetermined[0]), "ambiguous_both_alleles");
    }

    #[test]
    fn a_paralogous_reference_kmer_does_not_hide_a_real_variant() {
        // The mirror of the case above, and the one that matters clinically: a
        // copy of the marker's REFERENCE k-mer elsewhere used to read as "the
        // locus carries the reference", so a genome actually carrying the
        // variant was reported as not carrying it, with no row anywhere.
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        let k = 31;
        let index = generate_markerkmer(std::slice::from_ref(&marker), reference.as_str(), k, 10);
        let ref_rc = reverse_complement_sequence(&reference);
        let ref_kmer = index.markers()[0]
            .ref_kmer
            .clone()
            .expect("a substitution marker has a reference k-mer");

        // The variant, a neighbour that hides both k-mers at the locus, and a
        // copy of the reference k-mer pasted far away.
        let neighbour = other_base(&reference[304..305]);
        let mut genome = edited(&reference, &[(300, alt), (305, neighbour)]);
        genome.replace_range(500..500 + ref_kmer.len(), &ref_kmer);

        let lines = analyze_genome_seq("g", &genome, &index, &None, k, &reference, &ref_rc);
        assert!(
            lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")),
            "the anchors read the real locus, so the variant must still be reported: {:?}",
            lines.detail
        );
    }

    /// The WHO catalogue has five markers whose alternative k-mer H37Rv itself
    /// carries, two of them in PPE35, so the word turns up in every genome and
    /// the marker read as present everywhere, the ancestor included. When no
    /// anchor can tell the real locus from the copy, the only honest answer is
    /// that the marker was not read.
    #[test]
    fn an_alternative_kmer_the_reference_carries_is_not_evidence() {
        let base = long_ref();
        let pos = 300;
        let alt = other_base(&base[pos - 1..pos]);
        // A stretch holding the locus with the variant already in it, and both
        // of the marker's anchors, copied elsewhere in the reference. The
        // alternative k-mer is now in the reference and neither anchor is
        // unique, so nothing can say which copy a hit came from.
        let mut stretch = base[250..350].to_string();
        stretch.replace_range(pos - 1 - 250..pos - 250, alt);
        let mut reference = base.clone();
        reference.replace_range(450..550, &stretch);

        let marker = marker_at(&reference, pos, alt);
        let k = 31;
        let index = generate_markerkmer(std::slice::from_ref(&marker), reference.as_str(), k, 10);
        assert!(
            index.markers()[0].alt_in_reference,
            "the fixture needs the copy to be in the reference"
        );

        // A genome identical to the reference: its locus carries the reference
        // allele, so the marker is not there.
        let ref_rc = reverse_complement_sequence(&reference);
        let lines = analyze_genome_seq("g", &reference, &index, &None, k, &reference, &ref_rc);
        assert!(
            !lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")),
            "a word the reference carries must not be read as the marker being present: {:?}",
            lines.detail
        );
        assert!(
            lines.undetermined.iter().any(|l| l.contains("alt_kmer_in_reference")),
            "and the row must say why it could not be read: {:?}",
            lines.undetermined
        );
    }

    /// Sixty-seven markers of the catalogue have a second copy of their
    /// reference window. A reference hit may be that copy, so it cannot rule
    /// the marker out: a genome whose real locus carries the variant behind
    /// another change used to be reported as not carrying it at all.
    #[test]
    fn a_repeated_reference_kmer_cannot_rule_a_marker_out() {
        let base = long_ref();
        let pos = 300;
        let alt = other_base(&base[pos - 1..pos]);
        // The window and both anchors, copied as they are: the reference word
        // now occurs twice and no anchor is unique.
        let mut reference = base.clone();
        let stretch = base[250..350].to_string();
        reference.replace_range(450..550, &stretch);

        let marker = marker_at(&reference, pos, alt);
        let k = 31;
        let index = generate_markerkmer(std::slice::from_ref(&marker), reference.as_str(), k, 10);
        assert_eq!(index.markers()[0].ref_copies, 2, "the fixture needs the repeated window");
        assert!(!index.markers()[0].alt_in_reference, "but not the alternative word");

        // The genome carries the variant, with a second change seven bases
        // along that hides the alternative k-mer. The reference k-mer still
        // matches, at the copy.
        let mut genome = reference.clone();
        genome.replace_range(pos - 1..pos, alt);
        let neighbour = other_base(&base[pos + 6..pos + 7]);
        genome.replace_range(pos + 6..pos + 7, neighbour);

        let ref_rc = reverse_complement_sequence(&reference);
        let lines = analyze_genome_seq("g", &genome, &index, &None, k, &reference, &ref_rc);
        assert!(
            lines.undetermined.iter().any(|l| l.contains("ref_kmer_repeated")),
            "the copy must not be read as this locus holding the reference allele: {:?}",
            lines.undetermined
        );

        // A genome that shows both copies as the reference has proved the
        // marker absent: there is no third copy for it to be hiding in. This
        // is the ordinary case, and the one that would drown the undetermined
        // table if a repeated window were treated as unreadable by itself.
        let lines = analyze_genome_seq("g", &reference, &index, &None, k, &reference, &ref_rc);
        assert!(
            lines.undetermined.is_empty(),
            "a genome showing every copy as the reference must simply be absent: {:?}",
            lines.undetermined
        );
        assert!(
            !lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")),
            "and nothing is called either: {:?}",
            lines.detail
        );

        // The same with the neighbour on the other side, so each anchor in
        // turn is the one left looking unique while the real copy is broken.
        let mut genome = reference.clone();
        genome.replace_range(pos - 1..pos, alt);
        let neighbour = other_base(&base[pos - 8..pos - 7]);
        genome.replace_range(pos - 8..pos - 7, neighbour);
        let lines = analyze_genome_seq("g", &genome, &index, &None, k, &reference, &ref_rc);
        assert!(
            lines.undetermined.iter().any(|l| l.contains("ref_kmer_repeated")),
            "nor with the neighbour on the left: {:?}",
            lines.undetermined
        );
    }

    #[test]
    fn a_paralogous_alt_kmer_does_not_make_the_marker_present() {
        // The marker's ALT k-mer also occurs far away, as happens in the
        // resistance panel where a PE/PPE copy carries it. The genome holds
        // the reference allele at the marker's own locus, and an anchor, being
        // unique there, says so: the marker is absent, not present.
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        let k = 31;
        let index = generate_markerkmer(std::slice::from_ref(&marker), reference.as_str(), k, 10);
        let ref_rc = reverse_complement_sequence(&reference);
        // The reference genome, with the marker's own ALT k-mer pasted over a
        // distant stretch: the ALT k-mer now occurs once, far from the marker,
        // while the marker's own locus and its anchors are untouched.
        let paralogue = index.markers()[0].kmer.clone();
        let mut genome = reference.clone();
        genome.replace_range(500..500 + paralogue.len(), &paralogue);
        assert_eq!(&genome[299..300], &reference[299..300], "the locus itself is unchanged");

        let lines = analyze_genome_seq("g", &genome, &index, &None, k, &reference, &ref_rc);
        assert!(
            !lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")),
            "a paralogous hit must not be read as the marker being present"
        );
        assert!(
            lines.undetermined.is_empty(),
            "the anchors read the real locus, so the marker is absent, not unresolved: {:?}",
            lines.undetermined
        );
    }

    #[test]
    fn reading_one_allele_of_a_position_rules_out_the_others() {
        // Two markers at one position share their reference k-mer. A genome
        // carrying one of them has been read there, so the other is absent.
        let genome = edited(REF60, &[(30, "A")]);
        let lines = classify_one(&genome, &[snp_marker("A"), snp_marker("T")]);
        assert!(lines.detail.iter().any(|l| l.contains("L9.9;A")));
        assert!(!lines.detail.iter().any(|l| l.contains("L9.9;T")));
        assert!(lines.undetermined.is_empty());
    }

    #[test]
    fn the_anchors_read_the_same_allele_on_a_reverse_oriented_contig() {
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        let neighbour = other_base(&reference[304..305]);
        let forward = edited(&reference, &[(300, alt), (305, neighbour)]);
        let reverse = reverse_complement_sequence(&forward);

        for genome in [forward.as_str(), reverse.as_str()] {
            let lines = classify_long(genome, &reference, std::slice::from_ref(&marker));
            assert!(
                lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")),
                "contig orientation must not change the anchor read"
            );
        }
    }

    #[test]
    fn a_marker_without_a_usable_reference_kmer_is_never_called_absent() {
        // A reference allele the 2-bit index cannot hold leaves the marker with
        // an alternative k-mer only, so it can be seen present but never
        // confirmed absent from the k-mers. It still must not pass as absent.
        let marker = MarkerVariant {
            pos: 30,
            ref_allele: "N".to_string(),
            alt_allele: "A".to_string(),
            lineage: "L9.9".to_string(),
            gene: None,
            mutation: None,
        };
        let lines = classify_one(REF60, &[marker]);
        assert!(
            !lines.detail.iter().any(|l| l.contains("L9.9")),
            "a marker whose reference allele cannot be matched must not be called present"
        );
        // The anchors read the locus and it does not carry the marker's ALT,
        // so the marker is absent even though its own REF is unmatchable.
        assert!(lines.undetermined.is_empty());
    }

    // -----------------------------------------------------------------------
    // The anchors: reading the allele the k-mer could not
    // -----------------------------------------------------------------------

    /// A long reference, so a 15-base anchor is unique in it the way it is in
    /// a real genome. The marker sits at 1-based position 300.
    fn long_ref() -> String {
        // Deterministic pseudo-random sequence: a linear congruential walk.
        let mut seq = String::with_capacity(600);
        let mut x: u64 = 12345;
        for _ in 0..600 {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            seq.push(['A', 'C', 'G', 'T'][(x >> 33) as usize % 4]);
        }
        seq
    }

    fn classify_long(genome: &str, reference: &str, markers: &[MarkerVariant]) -> GenomeLines {
        let k = 31;
        let index = generate_markerkmer(markers, reference, k, 10);
        let ref_rc = reverse_complement_sequence(reference);
        analyze_genome_seq("g", genome, &index, &None, k, reference, &ref_rc)
    }

    fn marker_at(reference: &str, pos: usize, alt: &str) -> MarkerVariant {
        MarkerVariant {
            pos,
            ref_allele: reference[pos - 1..pos].to_string(),
            alt_allele: alt.to_string(),
            lineage: "RIF;RIF-R;test".to_string(),
            gene: None,
            mutation: None,
        }
    }

    /// A base other than the one already there.
    fn other_base(current: &str) -> &'static str {
        match current {
            "A" => "C",
            "C" => "G",
            "G" => "T",
            _ => "A",
        }
    }

    #[test]
    fn the_anchors_read_a_variant_the_kmer_cannot_see() {
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        // The variant, plus a neighbour five bases away, inside the k-mer
        // window but outside the left anchor.
        let neighbour = other_base(&reference[304..305]);
        let genome = edited(&reference, &[(300, alt), (305, neighbour)]);

        let lines = classify_long(&genome, &reference, std::slice::from_ref(&marker));
        let row = lines
            .detail
            .iter()
            .find(|l| l.contains("RIF;RIF-R;test"))
            .expect("the anchors must find the variant the k-mer missed");
        assert_eq!((alleles_found(row), evidence(row)), ("ALT", "left_anchor"));
        assert!(lines.undetermined.is_empty());
    }

    #[test]
    fn a_neighbour_on_the_left_leaves_the_right_anchor_to_read_it() {
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        let neighbour = other_base(&reference[294..295]);
        let genome = edited(&reference, &[(295, neighbour), (300, alt)]);

        let lines = classify_long(&genome, &reference, std::slice::from_ref(&marker));
        let row = lines.detail.iter().find(|l| l.contains("RIF;RIF-R;test")).expect("found");
        assert_eq!((alleles_found(row), evidence(row)), ("ALT", "right_anchor"));
    }

    #[test]
    fn the_anchors_rule_a_marker_out_instead_of_leaving_it_unread() {
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        // Only the neighbour: the genome keeps the reference allele, so the
        // marker is absent, not undetermined.
        let neighbour = other_base(&reference[304..305]);
        let genome = edited(&reference, &[(305, neighbour)]);

        let lines = classify_long(&genome, &reference, std::slice::from_ref(&marker));
        assert!(!lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")));
        assert!(
            lines.undetermined.is_empty(),
            "the anchors read the reference allele, so nothing is undetermined: {:?}",
            lines.undetermined
        );
    }

    #[test]
    fn a_locus_read_as_something_else_leaves_the_marker_absent() {
        // A third base at the marker's own position, which the panel does not
        // list. The locus was read, so the marker is settled: it is absent.
        // Finding what is there instead is a variant caller's job.
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        let third = other_base(alt);
        let neighbour = other_base(&reference[304..305]);
        let genome = edited(&reference, &[(300, third), (305, neighbour)]);

        let lines = classify_long(&genome, &reference, std::slice::from_ref(&marker));
        assert!(!lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")));
        assert!(
            lines.undetermined.is_empty(),
            "the locus was read, so nothing is undetermined: {:?}",
            lines.undetermined
        );
    }

    #[test]
    fn a_repeated_anchor_is_not_trusted() {
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let marker = marker_at(&reference, 300, alt);
        let neighbour = other_base(&reference[304..305]);
        let genome = edited(&reference, &[(300, alt), (305, neighbour)]);
        // The same contig twice over: every anchor now occurs twice, so none of
        // them says which copy was read.
        let doubled = format!("{genome}{genome}");

        let lines = classify_long(&doubled, &reference, std::slice::from_ref(&marker));
        assert!(!lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")));
        assert_eq!(reason(&lines.undetermined[0]), "not_found");
    }

    #[test]
    fn an_indel_keeps_the_kmer_answer() {
        let reference = long_ref();
        let marker = MarkerVariant {
            pos: 300,
            ref_allele: reference[299..301].to_string(),
            alt_allele: reference[299..300].to_string(),
            lineage: "RIF;RIF-R;del".to_string(),
            gene: None,
            mutation: None,
        };
        // A neighbour hides both k-mers. The allele's own length is what is in
        // question, so the anchors must not guess: the marker stays unread.
        let neighbour = other_base(&reference[309..310]);
        let genome = edited(&reference, &[(310, neighbour)]);
        let lines = classify_long(&genome, &reference, std::slice::from_ref(&marker));
        assert_eq!(reason(&lines.undetermined[0]), "not_found");
    }

    #[test]
    fn a_genome_read_whole_by_the_kmers_never_runs_the_anchor_pass() {
        // Nothing to resolve, so `collect_marker_lines` must not pay for the
        // second pass. Asserted through the result it would change: a clean
        // genome has no undetermined rows either way.
        let reference = long_ref();
        let alt = other_base(&reference[299..300]);
        let lines = classify_long(&reference, &reference, &[marker_at(&reference, 300, alt)]);
        assert!(lines.undetermined.is_empty());
        assert!(!lines.detail.iter().any(|l| l.contains("RIF;RIF-R;test")));
    }

    #[test]
    fn undetermined_rows_are_tallied_under_the_first_level_of_the_path() {
        let mut map: HashMap<String, HashMap<String, usize>> = HashMap::new();
        for line in [
            "g1\t761155\tC\tT\tRIF;RIF-R;rpoB_S450L\trpoB\tS450L\tnot_found",
            "g1\t761150\tA\tC\tRIF;RIF-R;rpoB_R448S\trpoB\tR448S\tnot_found",
            "g1\t1673425\tC\tT\tINH;INH-R;fabG1_c.-15C>T\tfabG1\tc.-15C>T\tnot_found",
        ] {
            super::update_undetermined_from_line(&mut map, line);
        }
        let (total, by_label) = super::undetermined_columns(map.get("g1"));
        assert_eq!(total, 3);
        assert_eq!(by_label, "RIF:2 INH:1");
        assert_eq!(super::undetermined_columns(map.get("g2")), (0, String::new()));
    }
}
