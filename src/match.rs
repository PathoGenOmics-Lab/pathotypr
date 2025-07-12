//! src/match.rs
//!
//! Implements the `match` subcommand.
//!
//! This module takes input FASTQ reads and compares them against a collection
//! of reference genomes provided in a single multi-FASTA file. It uses an
//! efficient, k-mer-based weighted containment score to determine the best
//! matching reference genome for the sample.

use crate::errors::{AppError, AppResult};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs;
use std::fs::File;
use std::io::Write;

/// Command-line arguments for the `match` subcommand.
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

    /// Path to a single multi-FASTA file containing all reference genomes.
    #[arg(short = 'r', long = "references", required = true)]
    pub references: String,

    /// Path for the output TSV report. Prints to console if not provided.
    #[arg(short = 'o', long)]
    pub output: Option<String>,

    /// k-mer size for comparison.
    #[arg(short = 'k', long, default_value_t = 31)]
    pub kmer_size: u8,

    /// Number of CPU threads to use. Defaults to all available.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
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
        // The format can be sample_name\tpath1\tpath2...
        // We collect all paths from the second column onwards.
        let paths: Vec<String> = line.split('\t').skip(1).map(String::from).collect();
        if paths.is_empty() {
            warn!("Skipping malformed line in TSV (no paths found): {}", line);
        }
        fastq_paths.extend(paths);
    }
    if fastq_paths.is_empty() {
        return Err(AppError::NotEnoughData(format!("No FASTQ files found in the TSV list: {}", path)));
    }
    Ok(fastq_paths)
}


/// Counts k-mer occurrences from a set of FASTQ files in parallel.
fn count_kmers_from_fastqs(paths: &[String], k: u8) -> AppResult<FxHashMap<u64, u32>> {
    info!("🧬 Counting k-mers from input FASTQ files...");
    let pb = ProgressBar::new(paths.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("  {spinner:.blue} Reading FASTQs: [{bar:40.cyan/blue}] {pos}/{len} {msg}")?
            .progress_chars("=>-"),
    );

    let kmer_counts = paths
        .par_iter()
        .map(|path| {
            let mut local_counts = FxHashMap::default();
            let mut reader = parse_fastx_file(path)
                .unwrap_or_else(|_| panic!("Failed to open or parse FASTQ file: {}", path));

            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid record in FASTQ");
                // Use bit_kmers for canonical k-mers, which is robust to strandness
                for (_, kmer, _) in seqrec.sequence().bit_kmers(k, true) {
                    *local_counts.entry(kmer.0).or_insert(0) += 1;
                }
            }
            pb.inc(1);
            local_counts
        })
        .reduce(FxHashMap::default, |mut a, b| {
            // Combine the k-mer counts from all parallel jobs
            for (kmer, count) in b {
                *a.entry(kmer).or_insert(0) += count;
            }
            a
        });
    
    pb.finish_with_message("FASTQ k-mer counting complete.");
    Ok(kmer_counts)
}

/// Extracts all reference sequences and their headers from a multi-FASTA file.
fn read_references_from_multifasta(path: &str) -> AppResult<Vec<(String, Vec<u8>)>> {
    info!("📖 Loading reference genomes from multi-FASTA file...");
    let mut refs = Vec::new();
    let mut reader = parse_fastx_file(path)
        .map_err(|_| AppError::Generic(format!("Cannot open reference file: {}", path)))?;

    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|_| AppError::Parsing("Invalid record in reference FASTA.".to_string()))?;
        let header = seqrec.id();
        let seq = seqrec.seq().to_vec();
        refs.push((String::from_utf8_lossy(header).to_string(), seq));
    }
    info!("  Found {} reference sequences to compare against.", refs.len());
    Ok(refs)
}

/// Main execution logic for the `match` subcommand.
pub fn run(args: Args) -> AppResult<()> {
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| AppError::Generic(format!("Failed to build thread pool: {}", e)))?;
    }

    // Determine the list of FASTQ files to process from either --input or --input-list
    let fastq_files_to_process = if let Some(list_path) = &args.input_list {
        read_fastq_list_from_tsv(list_path)?
    } else if !args.fastqs.is_empty() {
        args.fastqs
    } else {
        return Err(AppError::Generic(
            "An input source is required: either --input or --input-list must be provided.".to_string(),
        ));
    };

    // Log which files are being processed
    info!("Analyzing FASTQ files: {}", fastq_files_to_process.join(", "));

    // Step 1: Count k-mers from all input FASTQ files.
    let query_kmer_map = count_kmers_from_fastqs(&fastq_files_to_process, args.kmer_size)?;
    if query_kmer_map.is_empty() {
        return Err(AppError::NotEnoughData("No k-mers could be extracted from the input FASTQ files.".to_string()));
    }
    let total_query_kmers = query_kmer_map.values().sum::<u32>() as f64;
    debug!("Total k-mers counted in query: {}", total_query_kmers);

    // Step 2: Load all reference genomes from the multi-FASTA file.
    let references = read_references_from_multifasta(&args.references)?;

    // Step 3: Compare the query against each reference in parallel.
    info!("🔬 Comparing sample against references...");
    let pb = ProgressBar::new(references.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("  {spinner:.green} Matching: [{bar:40.magenta/purple}] {pos}/{len} ({eta})")?
            .progress_chars("=>-"),
    );

    let mut scores: Vec<(String, f64)> = references
        .par_iter()
        .map(|(header, seq)| {
            // Generate the set of unique k-mers for this reference.
            let ref_kmers: FxHashSet<u64> = seq.bit_kmers(args.kmer_size, true).map(|(_, kmer, _)| kmer.0).collect();

            // Calculate the weighted containment score.
            let shared_kmer_count: u32 = query_kmer_map
                .iter()
                .filter(|(kmer, _)| ref_kmers.contains(kmer))
                .map(|(_, count)| count)
                .sum();
            
            let score = if total_query_kmers > 0.0 {
                shared_kmer_count as f64 / total_query_kmers
            } else {
                0.0
            };
            pb.inc(1);
            (header.clone(), score)
        })
        .collect();

    pb.finish_with_message("Matching complete.");

    // Step 4: Sort and report only the best result.
    scores.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    
    info!("✅ Comparison finished. Writing report...");
    let mut writer: Box<dyn Write> = if let Some(path) = &args.output {
        Box::new(File::create(path)?)
    } else {
        Box::new(std::io::stdout())
    };

    writeln!(writer, "Query_Files\tBest_Match_Reference\tShared_Kmer_Fraction")?;
    if let Some((best_header, best_score)) = scores.first() {
        let query_files_str = fastq_files_to_process.join(",");
        writeln!(writer, "{}\t{}\t{:.4}", query_files_str, best_header, best_score)?;
    } else {
        info!("No matching references found.");
    }
    
    Ok(())
}
