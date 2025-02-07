use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use clap::Parser;
use csv::ReaderBuilder;
use log::{error, info};
use rayon::prelude::*;
use bio::io::fasta;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author = "Paula Ruiz Rodriguez <paula.ruiz-rodriguez@csic.es>",
          version = "0.1.0",
          about = "Classify genomes with markers into different lineages")]
struct Args {
    /// File with positions that are markers (TSV file)
    #[arg(short = 'p', long = "tsv_pos")]
    tsv_pos: String,

    /// Reference genome in FASTA format
    #[arg(short = 'r', long = "ref_fasta")]
    ref_fasta: String,

    /// TSV file with genome names and paths to FASTA files (required unless --fasta_genomes is provided)
    #[arg(short = 'g', long = "tsv_genomes", required_unless_present = "fasta_genomes")]
    tsv_genomes: Option<String>,

    /// FASTA file with one or multiple genomes (required unless --tsv_genomes is provided)
    #[arg(short = 'f', long = "fasta_genomes", required_unless_present = "tsv_genomes")]
    fasta_genomes: Option<String>,

    /// Output file
    #[arg(short = 'o', long = "output")]
    ofile: String,

    /// Number of CPUs to use (optional)
    #[arg(short = 'c', long = "num_cpu")]
    num_cpu: Option<usize>,
}

/// Generates all kmers of length `k` from the provided sequence.
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

/// Finds markers in the genome sequence by comparing kmers.
/// For each marker found, it returns a tuple of (genome kmer starting position, reference position, lineage).
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

/// Reads the markers TSV file and returns two maps:
/// - reference_positions: marker position (1-indexed) -> alternative base
/// - markers_lineage: marker position -> lineage
fn get_positions(tsv_file: &str) -> Result<(HashMap<usize, String>, HashMap<usize, String>), Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tsv_file)?;
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

/// Reads the first record from a FASTA file and returns the sequence (uppercased).
fn get_ref(fasta_file: &str) -> Result<String, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(fasta_file)?;
    let mut records = reader.records();
    if let Some(result) = records.next() {
        let record = result?;
        let seq = String::from_utf8(record.seq().to_vec())?;
        Ok(seq.to_uppercase())
    } else {
        Err("No record found in FASTA file".into())
    }
}

/// Generates a kmer for each marker by extracting a window around the marker position in the reference
/// and replacing the middle base with the alternative base. The marker positions are assumed to be 1-indexed.
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

/// Reads the genomes TSV file (genome name and FASTA path) and returns a map from genome name to its FASTA path.
fn get_genomepaths(tsv_file: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(tsv_file)?;
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

/// Analyzes a single genome from a file: reads its FASTA file, generates kmers,
/// and finds any markers that match. Returns a vector of formatted output lines.
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
        for (kmer, (position, ref_position, lineage)) in matched_markers.into_iter() {
            let snp_position = position + k / 2;
            result.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                genome_name, kmer, position, snp_position, ref_position, lineage
            ));
        }
    }
    result
}

/// Analyzes a single genome from a provided sequence.
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
        for (kmer, (position, ref_position, lineage)) in matched_markers.into_iter() {
            let snp_position = position + k / 2;
            result.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                genome_name, kmer, position, snp_position, ref_position, lineage
            ));
        }
    }
    result
}

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::init();
    let args = Args::parse();
    let tsv_pos = args.tsv_pos;
    let ref_fasta = args.ref_fasta;
    let output_file = args.ofile;
    let num_cpu = args.num_cpu;

    // Optionally set the number of threads for Rayon.
    if let Some(n) = num_cpu {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap();
    }

    // Check that required input files exist.
    if !Path::new(&tsv_pos).exists() {
        error!("The file {} does not exist", tsv_pos);
        return Err(format!("File {} does not exist", tsv_pos).into());
    }
    if !Path::new(&ref_fasta).exists() {
        error!("The file {} does not exist", ref_fasta);
        return Err(format!("File {} does not exist", ref_fasta).into());
    }
    if let Some(ref tsv_genomes) = args.tsv_genomes {
        if !Path::new(tsv_genomes).exists() {
            error!("The file {} does not exist", tsv_genomes);
            return Err(format!("File {} does not exist", tsv_genomes).into());
        }
    }
    if let Some(ref fasta_genomes) = args.fasta_genomes {
        if !Path::new(fasta_genomes).exists() {
            error!("The file {} does not exist", fasta_genomes);
            return Err(format!("File {} does not exist", fasta_genomes).into());
        }
    }

    // Read the reference genome sequence.
    let ref_seq = get_ref(&ref_fasta)?;
    // Read marker positions and lineages from the TSV file.
    let (reference_positions, markers_lineage) = get_positions(&tsv_pos)?;
    // Set kmer length.
    let k: usize = 21;
    // Generate marker kmers from the reference genome.
    let markers_kmers = generate_markerkmer(&reference_positions, &ref_seq, &markers_lineage, k);

    // Open the output file and write the header.
    let mut outfile = File::create(&output_file)?;
    writeln!(outfile, "genome\tk-mer\tk-merPOS\tSNPgenome\tSNPreference\tlineage")?;

    let mut results = Vec::new();

    // Process genomes based on the input type.
    if let Some(tsv_genomes) = args.tsv_genomes {
        // Process genomes using the TSV file with genome paths.
        let genome_paths = get_genomepaths(&tsv_genomes)?;
        let res: Vec<String> = genome_paths
            .par_iter()
            .flat_map_iter(|(genome_name, fasta_path)| {
                analyze_genome(genome_name, fasta_path, &markers_kmers, k)
            })
            .collect();
        results.extend(res);
    } else if let Some(fasta_genomes) = args.fasta_genomes {
        // Process genomes directly from a FASTA file (single or multi-FASTA).
        let genomes = get_genomes_from_fasta(&fasta_genomes)?;
        let res: Vec<String> = genomes
            .par_iter()
            .flat_map_iter(|(genome_name, genome_seq)| {
                analyze_genome_seq(genome_name, genome_seq, &markers_kmers, k)
            })
            .collect();
        results.extend(res);
    }

    // Write the results to the output file.
    for line in results {
        write!(outfile, "{}", line)?;
    }

    Ok(())
}
