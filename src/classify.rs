use clap::Parser;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{error, info};
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use bio::io::fasta;

/// Command-line arguments for the classify subcommand.
#[derive(Parser, Debug)]
pub struct Args {
    /// TSV file with marker positions
    #[arg(short = 'p', long = "tsv_pos")]
    pub tsv_pos: String,

    /// Reference genome in FASTA format
    #[arg(short = 'r', long = "ref_fasta")]
    pub ref_fasta: String,

    /// TSV file with genome names and paths to FASTA files (required unless --fasta_genomes is provided)
    #[arg(short = 'g', long = "tsv_genomes", required_unless_present = "fasta_genomes")]
    pub tsv_genomes: Option<String>,

    /// FASTA file with one or multiple genomes (required unless --tsv_genomes is provided)
    #[arg(short = 'f', long = "fasta_genomes", required_unless_present = "tsv_genomes")]
    pub fasta_genomes: Option<String>,

    /// Base name for the output file (the main per-marker file and a second summary TSV)
    #[arg(short = 'o', long = "output")]
    pub ofile: String,

    /// k-mer size (default is 21)
    #[arg(long, default_value_t = 21)]
    pub kmer_size: usize,

    /// Number of CPUs to use (optional)
    #[arg(short = 'c', long = "num_cpu")]
    pub num_cpu: Option<usize>,
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
            matched_markers.insert(marker_kmer.clone(), (genome_position, *ref_position, lineage.clone()));
        }
    }
    matched_markers
}

/// Reads the marker positions TSV file and returns two maps.
fn get_positions(tsv_file: &str) -> Result<(HashMap<usize, String>, HashMap<usize, String>), Box<dyn Error>> {
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

/// Reads the genomes TSV file and returns a map from genome name to its FASTA path.
fn get_genomepaths(tsv_file: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file)?;
    let mut genome_paths = HashMap::new();
    for result in rdr.records() {
        let record = result?;
        if record.len() < 2 {
            continue;
        }
        let genome_name = record[0].to_string();
        let fasta_path = record[1].to_string();
        if !Path::new(&fasta_path).exists() {
            eprintln!("The file {} does not exist", fasta_path);
            continue;
        }
        genome_paths.insert(genome_name, fasta_path);
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
    let matched_markers = find_markers(&genome_seq, markers_kmers, k);
    if matched_markers.is_empty() {
        result.push(format!("{}\t\t\t\t\t\n", genome_name));
    } else {
        for (kmer, (position, ref_position, lineage)) in matched_markers {
            let snp_position = position + k / 2;
            result.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                genome_name, kmer, position, snp_position, ref_position, lineage
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
    k: usize,
) -> Vec<String> {
    let mut result = Vec::new();
    info!("Analyzing genome {} (provided sequence)", genome_name);
    let matched_markers = find_markers(genome_seq, markers_kmers, k);
    if matched_markers.is_empty() {
        result.push(format!("{}\t\t\t\t\t\n", genome_name));
    } else {
        for (kmer, (position, ref_position, lineage)) in matched_markers {
            let snp_position = position + k / 2;
            result.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                genome_name, kmer, position, snp_position, ref_position, lineage
            ));
        }
    }
    result
}

/// Checks that all required input files exist.
fn check_input_files(args: &Args) -> Result<(), Box<dyn Error>> {
    let all_files = vec![
        (args.tsv_pos.as_str(), "TSV positions file"),
        (args.ref_fasta.as_str(), "Reference FASTA file"),
        (args.tsv_genomes.as_deref().unwrap_or(""), "TSV genomes file"),
        (args.fasta_genomes.as_deref().unwrap_or(""), "FASTA genomes file"),
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
    k: usize
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
            .map(|(genome_name, fasta_path)| {
                let lines = analyze_genome(genome_name, fasta_path, markers_kmers, k);
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
                let lines = analyze_genome_seq(genome_name, genome_seq, markers_kmers, k);
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

/// Main function for the classify subcommand.
fn generate_summary(results: &[String]) -> HashMap<String, HashMap<String, usize>> {
    let mut lineage_count_map = HashMap::new();
    
    for line in results {
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 6 {
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

fn write_summary(summary_out: &mut BufWriter<File>, genome: String, lineage_counts: Vec<(String, usize)>) -> Result<(), Box<dyn Error>> {
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

pub fn run(args: Args) -> Result<(), Box<dyn Error>> {
    env_logger::init();
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
        "genome\tk-mer\tk-merPOS\tSNPgenome\tSNPreference\tlineage"
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
