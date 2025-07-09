//! This module handles the `classify` subcommand.
//!
//! It performs lineage typing by identifying known genetic markers within assembled
//! genome sequences. The process involves:
//! 1. Reading a list of variants (SNPs, MNVs, Indels) and their corresponding lineages.
//! 2. Dynamically generating a unique diagnostic k-mer for each variant's alternate allele.
//! 3. Scanning input genomes (from a list of FASTA files or a single multifasta)
//!    for the presence of these diagnostic marker k-mers in parallel.
//! 4. Reporting all found markers for each genome and its annotations.
//! 5. Generating a summary file that lists the most likely lineage(s) for each genome.

use bio::io::fasta;
use bio::io::gff::{GffType, Reader};
use bio_types::strand::Strand as BioStrand; // Renamed to avoid conflict
use clap::Parser;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{error, info, warn}; // Import 'warn'
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

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

    /// One or more FASTA files to analyze.
    #[arg(short = 'i', long = "input", required_unless_present = "tsv_genomes")]
    pub fasta_genomes: Option<String>,
    
    /// Optional GFF file for annotation when using --input.
    #[arg(long = "gff", requires = "fasta_genomes")]
    pub gff_file: Option<String>,

    /// Prefix for the output files.
    #[arg(short = 'o', long = "output-prefix")]
    pub ofile: String,

    /// k-mer size (default is 21)
    #[arg(long, default_value_t = 21)]
    pub kmer_size: usize,

    /// Number of threads for parallel processing.
    #[arg(short = 't', long = "threads")]
    pub num_cpu: Option<usize>,
}

// Local Strand enum that derives Eq, required by rust-lapper
#[derive(Clone, Debug, PartialEq, Eq)]
enum Strand {
    Forward,
    Reverse,
    Unknown,
}

/// Structure to hold relevant gene information.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Gene {
    id: String,
    strand: Strand,
}

/// Standard genetic code map (3-letter).
fn genetic_code_3_letter() -> HashMap<&'static str, &'static str> {
    [
        ("GCA", "Ala"), ("GCC", "Ala"), ("GCG", "Ala"), ("GCT", "Ala"),
        ("CGA", "Arg"), ("CGC", "Arg"), ("CGG", "Arg"), ("CGT", "Arg"), ("AGA", "Arg"), ("AGG", "Arg"),
        ("AAC", "Asn"), ("AAT", "Asn"),
        ("GAC", "Asp"), ("GAT", "Asp"),
        ("TGC", "Cys"), ("TGT", "Cys"),
        ("GAA", "Glu"), ("GAG", "Glu"),
        ("CAA", "Gln"), ("CAG", "Gln"),
        ("GGA", "Gly"), ("GGC", "Gly"), ("GGG", "Gly"), ("GGT", "Gly"),
        ("CAC", "His"), ("CAT", "His"),
        ("ATA", "Ile"), ("ATC", "Ile"), ("ATT", "Ile"),
        ("CTA", "Leu"), ("CTC", "Leu"), ("CTG", "Leu"), ("CTT", "Leu"), ("TTA", "Leu"), ("TTG", "Leu"),
        ("AAA", "Lys"), ("AAG", "Lys"),
        ("ATG", "Met"),
        ("TTC", "Phe"), ("TTT", "Phe"),
        ("CCA", "Pro"), ("CCC", "Pro"), ("CCG", "Pro"), ("CCT", "Pro"),
        ("TCA", "Ser"), ("TCC", "Ser"), ("TCG", "Ser"), ("TCT", "Ser"), ("AGC", "Ser"), ("AGT", "Ser"),
        ("ACA", "Thr"), ("ACC", "Thr"), ("ACG", "Thr"), ("ACT", "Thr"),
        ("TGG", "Trp"),
        ("TAC", "Tyr"), ("TAT", "Tyr"),
        ("GTA", "Val"), ("GTC", "Val"), ("GTG", "Val"), ("GTT", "Val"),
        ("TAA", "Stp"), ("TAG", "Stp"), ("TGA", "Stp"),
    ]
    .iter().cloned().collect()
}

/// Generates all k-mers of length `k` from the given sequence.
fn generate_kmers(sequence: &str, k: usize) -> HashMap<String, usize> {
    let mut kmers = HashMap::new();
    let seq_len = sequence.len();
    if seq_len < k {
        return kmers;
    }
    for i in 0..=(seq_len - k) {
        let kmer = &sequence[i..i + k];
        kmers.insert(kmer.to_string(), i);
    }
    kmers
}

/// Finds markers in the genome sequence by comparing k-mers.
/// For each marker found, returns a tuple of (genome k-mer starting position, reference position, lineage).
fn find_markers(
    genome_sequence: &str,
    markers_kmers: &HashMap<String, (usize, String)>,
    k: usize,
) -> HashMap<String, (usize, usize, String)> {
    let genome_kmers = generate_kmers(genome_sequence, k);
    let mut matched_markers = HashMap::new();
    for (marker_kmer, (ref_position, lineage)) in markers_kmers.iter() {
        if let Some(&genome_position) = genome_kmers.get(marker_kmer) {
            matched_markers.insert(
                marker_kmer.clone(),
                (genome_position, *ref_position, lineage.clone()),
            );
        }
    }
    matched_markers
}

/// Reads the marker positions TSV file and returns two maps.
fn get_positions(
    tsv_file: &str,
) -> Result<(HashMap<usize, String>, HashMap<usize, String>), Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file)?;
    let mut reference_positions = HashMap::new();
    let mut markers_lineage = HashMap::new();
    for result in rdr.records() {
        let record = result?;
        if record.len() < 3 {
            continue;
        }
        let pos: usize = record[0].parse()?;
        let alt_base = record[1].to_string();
        let lineage = record[2].to_string();
        reference_positions.insert(pos, alt_base);
        markers_lineage.insert(pos, lineage);
    }
    Ok((reference_positions, markers_lineage))
}

/// Reads the first record from the reference FASTA file.
fn get_ref(fasta_file: &str) -> Result<String, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(fasta_file)?;
    let mut records = reader.records();
    if let Some(result) = records.next() {
        let record = result?;
        let seq = String::from_utf8(record.seq().to_vec())?;
        Ok(seq.to_uppercase())
    } else {
        Err("No record found in the reference FASTA file.".into())
    }
}

/// Generates a marker k-mer for each marker by extracting a window around the marker position in the reference,
/// replacing the middle base with the alternative base.
fn generate_markerkmer(
    reference_positions: &HashMap<usize, String>,
    ref_seq: &str,
    markers_lineage: &HashMap<usize, String>,
    k: usize,
) -> HashMap<String, (usize, String)> {
    let mut markers_kmers = HashMap::new();
    let seq_len = ref_seq.len();
    let half = k / 2;
    for (&pos, alt_base) in reference_positions.iter() {
        if pos > half && pos < seq_len - half {
            let start = pos - half - 1;
            let end = pos + half;
            if end <= seq_len {
                let kmer = ref_seq[start..end].to_string();
                let mut kmer_chars: Vec<char> = kmer.chars().collect();
                if half < kmer_chars.len() {
                    kmer_chars[half] = alt_base.chars().next().unwrap();
                }
                let new_kmer: String = kmer_chars.into_iter().collect();
                if let Some(lineage) = markers_lineage.get(&pos) {
                    markers_kmers.insert(new_kmer, (pos, lineage.clone()));
                }
            }
        }
    }
    markers_kmers
}

/// Reads the genomes TSV file and returns a map from genome name to its FASTA and optional GFF path.
fn get_genomepaths(tsv_file: &str) -> Result<HashMap<String, (String, Option<String>)>, Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file)?;
    let mut genome_paths = HashMap::new();
    for result in rdr.records() {
        let record = result?;
        if record.len() < 2 {
            continue;
        }
        let genome_name = record[0].to_string();
        let fasta_path = record[1].to_string();
        let gff_path = record
            .get(2)
            .and_then(|p| if p.is_empty() { None } else { Some(p.to_string()) });

        if !Path::new(&fasta_path).exists() {
            error!("The file {} does not exist", fasta_path);
            continue;
        }
        if let Some(gff) = &gff_path {
            if !Path::new(gff).exists() {
                error!("The GFF file {} does not exist", gff);
                continue;
            }
        }
        genome_paths.insert(genome_name, (fasta_path, gff_path));
    }
    Ok(genome_paths)
}

/// Reads a FASTA file (single or multi-FASTA) and returns a map from record ID to its sequence.
fn get_genomes_from_fasta(fasta_file: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(fasta_file)?;
    let mut genomes = HashMap::new();
    for result in reader.records() {
        let record = result?;
        let seq = String::from_utf8(record.seq().to_vec())?;
        genomes.insert(record.id().to_string(), seq.to_uppercase());
    }
    Ok(genomes)
}

/// Analyzes a single genome from a FASTA file, returning lines with marker details.
fn analyze_genome(
    genome_name: &str,
    fasta_path: &str,
    gff_path: &Option<String>,
    markers_kmers: &HashMap<String, (usize, String)>,
    k: usize,
) -> Vec<String> {
    let mut result = Vec::new();
    if !Path::new(fasta_path).exists() {
        error!("The file {} does not exist", fasta_path);
        return result;
    }
    info!("Analyzing genome {} ({})", genome_name, fasta_path);
    let genome_seq = match get_ref(fasta_path) {
        Ok(seq) => seq,
        Err(e) => {
            error!("Error reading FASTA file {}: {}", fasta_path, e);
            return result;
        }
    };

    // Parse GFF and build interval tree if available
    let annotations = if let Some(gff) = gff_path {
        match parse_gff_and_build_tree(gff) {
            Ok(tree) => Some(tree),
            Err(e) => {
                error!("Error processing GFF file {}: {}", gff, e);
                None
            }
        }
    } else {
        None
    };

    let matched_markers = find_markers(&genome_seq, markers_kmers, k);
    if matched_markers.is_empty() {
        result.push(format!("{}\t\t\t\t\t\t\t\t\n", genome_name)); // Add extra tabs for new columns
    } else {
        for (kmer, (position, ref_position, lineage)) in matched_markers {
            let snp_position = position + k / 2;
            let (mut gene_id, mut aa_pos, mut aa_change) = (None, None, None);

            if let Some(tree) = &annotations {
                // Find gene overlapping with the SNP
                for gene_interval in tree.find(snp_position, snp_position + 1) {
                    let alt_base = &genome_seq[snp_position..snp_position + 1];
                    let (gid, pos, change) =
                        translate_snp_info(gene_interval, snp_position, alt_base, &genome_seq);
                    gene_id = gid;
                    aa_pos = pos;
                    aa_change = change;
                    break; // Take the first annotation found
                }
            }

            result.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                genome_name, kmer, position, snp_position, ref_position, lineage,
                gene_id.as_deref().unwrap_or(""),
                aa_pos.as_deref().unwrap_or(""),
                aa_change.as_deref().unwrap_or("")
            ));
        }
    }
    result
}

/// Analyzes a single genome from a provided sequence (in memory), returning lines with marker details.
fn analyze_genome_seq(
    genome_name: &str,
    genome_seq: &str,
    markers_kmers: &HashMap<String, (usize, String)>,
    annotations: &Option<Lapper<usize, Gene>>,
    k: usize,
) -> Vec<String> {
    let mut result = Vec::new();
    info!("Analyzing genome {} (from provided sequence)", genome_name);
    let matched_markers = find_markers(genome_seq, markers_kmers, k);
    if matched_markers.is_empty() {
        result.push(format!("{}\t\t\t\t\t\t\t\t\n", genome_name));
    } else {
        for (kmer, (position, ref_position, lineage)) in matched_markers {
            let snp_position = position + k / 2;
            let (mut gene_id, mut aa_pos, mut aa_change) = (None, None, None);

            if let Some(tree) = annotations {
                for gene_interval in tree.find(snp_position, snp_position + 1) {
                    let alt_base = &genome_seq[snp_position..snp_position + 1];
                    let (gid, pos, change) =
                        translate_snp_info(gene_interval, snp_position, alt_base, genome_seq);
                    gene_id = gid;
                    aa_pos = pos;
                    aa_change = change;
                    break;
                }
            }
            
            result.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                genome_name, kmer, position, snp_position, ref_position, lineage,
                gene_id.as_deref().unwrap_or(""),
                aa_pos.as_deref().unwrap_or(""),
                aa_change.as_deref().unwrap_or("")
            ));
        }
    }
    result
}

/// Parses a GFF file and builds an Interval Tree with gene features ('CDS').
fn parse_gff_and_build_tree(gff_file: &str) -> Result<Lapper<usize, Gene>, Box<dyn Error>> {
    let mut reader = Reader::new(File::open(gff_file)?, GffType::GFF3);
    let mut intervals = Vec::new();
    for result in reader.records() {
        match result {
            Ok(record) => {
                // We focus on coding sequences (CDS)
                if record.feature_type() == "CDS" {
                    let id = record
                        .attributes()
                        .get("ID")
                        .or_else(|| record.attributes().get("Name"))
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| "Unknown".to_string());
                    let gene = Gene {
                        id,
                        strand: match record.strand() {
                            Some(BioStrand::Forward) => Strand::Forward,
                            Some(BioStrand::Reverse) => Strand::Reverse,
                            _ => Strand::Unknown,
                        },
                    };
                    intervals.push(Interval {
                        start: *record.start() as usize - 1, // GFF is 1-based, lapper is 0-based
                        stop: *record.end() as usize,
                        val: gene,
                    });
                }
            }
            Err(e) => {
                // If a line is malformed, print a warning and skip it.
                warn!("Skipping malformed GFF record: {}", e);
            }
        }
    }
    let lapper = Lapper::new(intervals);
    Ok(lapper)
}

/// Translates SNP information into an amino acid change. Returns (GeneID, AAPosition, AAChange).
fn translate_snp_info(
    gene_interval: &Interval<usize, Gene>,
    snp_pos: usize,
    alt_base: &str,
    sequence: &str,
) -> (Option<String>, Option<String>, Option<String>) {
    let gene = &gene_interval.val;
    let code = genetic_code_3_letter();

    let (_ref_codon_seq, alt_codon_seq, aa_pos) = match gene.strand {
        Strand::Forward => {
            let frame_start = (snp_pos - gene_interval.start) % 3;
            let codon_start = snp_pos - frame_start;
            let aa_pos = (snp_pos - gene_interval.start) / 3 + 1;

            if codon_start + 3 > sequence.len() {
                return (None, None, None);
            }

            let mut ref_codon: Vec<char> = sequence[codon_start..codon_start + 3].chars().collect();
            let mut alt_codon = ref_codon.clone();
            alt_codon[frame_start] = alt_base.chars().next().unwrap_or('N');

            (
                ref_codon.into_iter().collect::<String>(),
                alt_codon.into_iter().collect::<String>(),
                aa_pos,
            )
        }
        Strand::Reverse => {
            // For the reverse strand, we first get the reverse complement
            let sequence_rc = sequence
                .chars()
                .rev()
                .map(|b| match b {
                    'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N',
                })
                .collect::<String>();

            // Coordinates on the RC strand are inverted
            let snp_pos_rc = sequence.len() - snp_pos - 1;
            let start_rc = sequence.len() - gene_interval.stop;

            let frame_start = (snp_pos_rc - start_rc) % 3;
            let codon_start = snp_pos_rc - frame_start;
            let aa_pos = (snp_pos_rc - start_rc) / 3 + 1;

            if codon_start + 3 > sequence_rc.len() {
                return (None, None, None);
            }

            let mut ref_codon: Vec<char> = sequence_rc[codon_start..codon_start + 3].chars().collect();
            let mut alt_codon = ref_codon.clone();
            let alt_base_rc = match alt_base.chars().next().unwrap_or('N') {
                'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N',
            };
            alt_codon[frame_start] = alt_base_rc;

            (
                ref_codon.into_iter().collect::<String>(),
                alt_codon.into_iter().collect::<String>(),
                aa_pos,
            )
        }
        _ => return (None, None, None), // Strand::Unknown
    };

    let alt_aa_3l = code.get(alt_codon_seq.as_str()).unwrap_or(&"???");

    let gene_id = Some(gene.id.clone());
    let aa_pos_str = Some(aa_pos.to_string());
    let aa_change_str = Some(alt_aa_3l.to_string());

    (gene_id, aa_pos_str, aa_change_str)
}

/// Checks that all required input files exist.
fn check_input_files(args: &Args) -> Result<(), Box<dyn Error>> {
    let all_files = vec![
        (args.tsv_pos.as_str(), "TSV positions file"),
        (args.ref_fasta.as_str(), "Reference FASTA file"),
        (
            args.tsv_genomes.as_deref().unwrap_or(""),
            "TSV genomes file",
        ),
        (
            args.fasta_genomes.as_deref().unwrap_or(""),
            "FASTA genomes file",
        ),
    ];
    for (path, description) in all_files.into_iter().filter(|(p, _)| !p.is_empty()) {
        if !Path::new(path).exists() {
            error!("The {} {} does not exist", description, path);
            return Err(format!("File {} does not exist", path).into());
        }
    }
    Ok(())
}

/// Processes genomes and returns lines with marker matches for each genome.
fn process_genomes(
    args: &Args,
    markers_kmers: &HashMap<String, (usize, String)>,
    k: usize,
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut results = Vec::new();

    // If we have a TSV of genome names -> paths
    if let Some(tsv_genomes) = &args.tsv_genomes {
        let genome_paths = get_genomepaths(tsv_genomes)?;
        let pb = ProgressBar::new(genome_paths.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("=>-"),
        );
        let res: Vec<String> = genome_paths
            .par_iter()
            .map(|(genome_name, (fasta_path, gff_path))| {
                let lines = analyze_genome(genome_name, fasta_path, gff_path, markers_kmers, k);
                pb.inc(1);
                lines
            })
            .flatten()
            .collect();
        pb.finish_with_message("Done!");
        results.extend(res);
    }
    // Otherwise, we have a multi-FASTA with possibly multiple genomes
    else if let Some(fasta_genomes) = &args.fasta_genomes {
        let annotations = if let Some(gff_path) = &args.gff_file {
            match parse_gff_and_build_tree(gff_path) {
                Ok(tree) => Some(tree),
                Err(e) => {
                    error!("Error processing GFF file {}: {}", gff_path, e);
                    None
                }
            }
        } else {
            None
        };

        let genomes = get_genomes_from_fasta(fasta_genomes)?;
        let pb = ProgressBar::new(genomes.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("=>-"),
        );
        let res: Vec<String> = genomes
            .par_iter()
            .map(|(genome_name, genome_seq)| {
                let lines = analyze_genome_seq(genome_name, genome_seq, markers_kmers, &annotations, k);
                pb.inc(1);
                lines
            })
            .flatten()
            .collect();
        pb.finish_with_message("Done!");
        results.extend(res);
    }
    Ok(results)
}

/// Generates a summary of lineage counts per genome.
fn generate_summary(results: &[String]) -> HashMap<String, HashMap<String, usize>> {
    let mut lineage_count_map = HashMap::new();

    for line in results {
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 6 || fields[5].is_empty() {
            continue;
        }
        let genome_name = fields[0].to_string();
        let lineage = fields[5].to_string();

        lineage_count_map
            .entry(genome_name)
            .or_insert_with(HashMap::new)
            .entry(lineage)
            .and_modify(|c| *c += 1)
            .or_insert(1);
    }

    lineage_count_map
}

/// Writes a single summary line for a genome.
fn write_summary(
    summary_out: &mut BufWriter<File>,
    genome: String,
    lineage_counts: Vec<(String, usize)>,
) -> Result<(), Box<dyn Error>> {
    let lineage_count_str = lineage_counts
        .iter()
        .map(|(lin, cnt)| format!("{}:{}", lin, cnt))
        .collect::<Vec<_>>()
        .join(",");

    let majority_lineage = if lineage_counts.is_empty() {
        String::new()
    } else if lineage_counts.len() == 1 {
        lineage_counts[0].0.clone()
    } else if lineage_counts[0].1 > lineage_counts[1].1 {
        lineage_counts[0].0.clone()
    } else {
        lineage_counts
            .iter()
            .map(|(lin, _)| lin.clone())
            .collect::<Vec<_>>()
            .join(",")
    };

    writeln!(
        summary_out,
        "{}\t{}\t{}",
        genome, lineage_count_str, majority_lineage
    )?;
    Ok(())
}

/// Main function for the classify subcommand.
pub fn run(args: Args) -> Result<(), Box<dyn Error>> {
    // env_logger::init(); // This should be called in main.rs
    if let Some(n) = args.num_cpu {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap();
    }
    check_input_files(&args)?;

    let k = args.kmer_size;
    let ref_seq = get_ref(&args.ref_fasta)?;
    let (reference_positions, markers_lineage) = get_positions(&args.tsv_pos)?;
    let markers_kmers = generate_markerkmer(&reference_positions, &ref_seq, &markers_lineage, k);
    let results = process_genomes(&args, &markers_kmers, k)?;

    let mut outfile = BufWriter::new(File::create(&args.ofile)?);
    writeln!(
        outfile,
        "genome\tk-mer\tk-merPOS\tSNPgenome\tSNPreference\tlineage\tGene\tAA_Pos\tAA_Change"
    )?;
    for line in &results {
        write!(outfile, "{}", line)?;
    }

    let summary_file = format!("{}_summary.tsv", args.ofile);
    let mut summary_out = BufWriter::new(File::create(&summary_file)?);
    writeln!(summary_out, "genome\tlineage:count\tmajor_lineage")?;

    let lineage_count_map = generate_summary(&results);

    for (genome, lineage_map) in lineage_count_map {
        let mut lineage_counts: Vec<(String, usize)> = lineage_map.into_iter().collect();
        lineage_counts.sort_by(|a, b| b.1.cmp(&a.1));
        write_summary(&mut summary_out, genome, lineage_counts)?;
    }

    Ok(())
}
    info!("✅ Process completed.");
    Ok(())
}
