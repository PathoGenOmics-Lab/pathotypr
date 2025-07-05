//! This module handles the `classify` subcommand.
//!
//! It performs lineage typing by identifying known genetic markers (k-mers)
//! within genome sequences. The process involves:
//! 1. Reading a list of positions and their corresponding alternate alleles and lineages.
//! 2. Generating "marker" k-mers from a reference sequence by substituting the
//!    alternate allele at each specified position.
//! 3. Scanning input genomes (from a list of FASTA files or a single multifasta)
//!    for the presence of these marker k-mers in parallel.
//! 4. Reporting all found markers for each genome.
//! 5. Generating a summary file that lists the most likely lineage(s) for each genome
//!    based on marker counts.

use anyhow::{anyhow, Result};
use bio::io::fasta;
use clap::Parser;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{error, info};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Parser, Debug)]
pub struct Args {
    /// Path to a TSV file defining markers.
    /// Format: position(1-based)\talt_base\tlineage
    #[arg(short = 'm', long)]
    pub markers: String,

    /// Path to the reference FASTA file.
    #[arg(short = 'r', long)]
    pub reference: String,

    /// Path to a TSV file listing genomes to analyze. Format: name\tpath/to/fasta.
    /// Use this or --genome-fasta.
    #[arg(long, group = "input_method")]
    pub genome_list: Option<String>,

    /// Path to a multifasta file containing all genomes to analyze.
    /// Use this or --genome-list.
    #[arg(long, group = "input_method")]
    pub genome_fasta: Option<String>,

    /// Path for the main output file (detailed marker report).
    #[arg(short = 'o', long)]
    pub output: String,

    /// The size of the k-mers to search for.
    #[arg(short = 'k', long, default_value_t = 21)]
    pub kmer_size: usize,

    /// Number of CPU threads to use for parallel processing.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
}

/// Generates all k-mers from a sequence and stores them in a HashMap for quick lookup.
fn generate_kmers(seq: &str, k: usize) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    if seq.len() < k {
        return m;
    }
    for i in 0..=seq.len() - k {
        m.insert(seq[i..i + k].to_string(), i);
    }
    m
}

/// Finds which of the provided markers are present in a given genome sequence.
pub fn find_markers(
    genome_seq: &str,
    markers: &HashMap<String, (usize, String)>,
    k: usize,
) -> HashMap<String, (usize, usize, String)> {
    let genome_kmers = generate_kmers(genome_seq, k);
    let mut found_markers = HashMap::new();
    for (marker_kmer, (ref_pos, lineage)) in markers {
        if let Some(&genome_pos) = genome_kmers.get(marker_kmer) {
            found_markers.insert(marker_kmer.clone(), (genome_pos, *ref_pos, lineage.clone()));
        }
    }
    found_markers
}

/// Reads a TSV file of marker positions.
pub fn get_positions(tsv: &str) -> Result<(HashMap<usize, String>, HashMap<usize, String>)> {
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(tsv)?;
    let mut ref_pos = HashMap::new();
    let mut lin_map = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        if rec.len() < 3 {
            continue;
        }
        let pos: usize = rec[0].parse()?;
        ref_pos.insert(pos, rec[1].to_string());
        lin_map.insert(pos, rec[2].to_string());
    }
    Ok((ref_pos, lin_map))
}

/// Generates the full set of marker k-mers from a reference sequence.
pub fn generate_marker_kmers(
    ref_pos: &HashMap<usize, String>,
    ref_seq: &str,
    lineage: &HashMap<usize, String>,
    k: usize,
) -> HashMap<String, (usize, String)> {
    let mut marker_map = HashMap::new();
    let half_k = k / 2;
    for (&pos, alt_base) in ref_pos {
        if pos <= half_k || pos + half_k >= ref_seq.len() {
            continue;
        }
        let start = pos - half_k - 1;
        let end = pos + half_k;
        let mut kmer_bytes = ref_seq[start..end].as_bytes().to_vec();
        kmer_bytes[half_k] = alt_base.as_bytes()[0];
        let kmer_str = String::from_utf8(kmer_bytes).unwrap();

        if let Some(lin) = lineage.get(&pos) {
            marker_map.insert(kmer_str, (pos, lin.clone()));
        }
    }
    marker_map
}

/// Reads the first sequence from a FASTA file.
fn get_ref(fa: &str) -> Result<String> {
    let rdr = fasta::Reader::from_file(fa)?;
    let mut recs = rdr.records();
    if let Some(r) = recs.next() {
        let r = r?;
        return Ok(String::from_utf8(r.seq().to_vec())?.to_uppercase());
    }
    Err(anyhow!("Reference FASTA is empty or invalid."))
}

/// Reads a TSV file mapping genome names to their FASTA file paths.
fn get_genome_paths(tsv: &str) -> Result<HashMap<String, String>> {
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(tsv)?;
    let mut m = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        if rec.len() < 2 { continue; }
        let name = rec[0].to_string();
        let path = rec[1].to_string();
        if !Path::new(&path).exists() {
            error!("Missing FASTA for genome {}: {}", name, path);
            continue;
        }
        m.insert(name, path);
    }
    Ok(m)
}

/// Reads all sequences from a multifasta file into a map.
fn get_genomes_from_fasta(path: &str) -> Result<HashMap<String, String>> {
    let rdr = fasta::Reader::from_file(path)?;
    let mut m = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let id = rec.id().to_string();
        let seq = String::from_utf8(rec.seq().to_vec())?.to_uppercase();
        m.insert(id, seq);
    }
    Ok(m)
}

/// Analyzes a single genome from a file path.
fn analyze_path(name: &str, path: &str, markers: &HashMap<String, (usize, String)>, k: usize) -> Vec<String> {
    let seq = match get_ref(path) {
        Ok(s) => s,
        Err(e) => { error!("Could not read {}: {}", path, e); return Vec::new(); }
    };
    analyze_seq(name, &seq, markers, k)
}

/// Analyzes a single genome sequence for markers and formats the output lines.
fn analyze_seq(name: &str, seq: &str, markers: &HashMap<String, (usize, String)>, k: usize) -> Vec<String> {
    let mut v = Vec::new();
    let found_markers = find_markers(seq, markers, k);
    if found_markers.is_empty() {
        v.push(format!("{}\t\t\t\t\t\n", name));
    } else {
        for (kmer, (genome_pos, ref_pos, lineage)) in found_markers {
            let snp_pos_in_genome = genome_pos + k / 2;
            v.push(format!("{}\t{}\t{}\t{}\t{}\t{}\n", name, kmer, genome_pos, snp_pos_in_genome, ref_pos, lineage));
        }
    }
    v
}

/// Summarizes marker counts per lineage for each genome.
fn lineage_summary(lines: &[String]) -> HashMap<String, HashMap<String, usize>> {
    let mut summary_map = HashMap::new();
    for line in lines {
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 6 { continue; }
        let genome_name = fields[0].to_string();
        let lineage = fields[5].to_string();
        summary_map.entry(genome_name).or_insert_with(HashMap::new).entry(lineage).and_modify(|count| *count += 1).or_insert(1);
    }
    summary_map
}

/// Writes a single summary line for a genome.
fn write_summary_line(w: &mut BufWriter<File>, genome: &str, mut lineage_counts: Vec<(String, usize)>) -> Result<()> {
    lineage_counts.sort_by(|a, b| b.1.cmp(&a.1));
    let list = lineage_counts.iter().map(|(l, c)| format!("{}:{}", l, c)).collect::<Vec<_>>().join(",");
    let major_lineage = if lineage_counts.is_empty() {
        "".to_string()
    } else if lineage_counts.len() == 1 || lineage_counts[0].1 > lineage_counts[1].1 {
        lineage_counts[0].0.clone()
    } else {
        let top_count = lineage_counts[0].1;
        lineage_counts.iter().filter(|(_, count)| *count == top_count).map(|(l, _)| l.clone()).collect::<Vec<_>>().join(",")
    };
    writeln!(w, "{}\t{}\t{}", genome, list, major_lineage)?;
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    if let Some(n) = args.threads { rayon::ThreadPoolBuilder::new().num_threads(n).build_global()?; }
    if args.genome_list.is_none() && args.genome_fasta.is_none() {
        return Err(anyhow!("You must provide an input source: either --genome-list or --genome-fasta."));
    }

    info!("▶ Preparing marker k-mers...");
    let k = args.kmer_size;
    let ref_seq = get_ref(&args.reference)?;
    let (pos_map, lin_map) = get_positions(&args.markers)?;
    let marker_kmer_map = generate_marker_kmers(&pos_map, &ref_seq, &lin_map, k);
    info!("  Loaded {} unique marker k-mers.", marker_kmer_map.len());

    info!("▶ Analyzing genomes...");
    let mut lines = Vec::<String>::new();
    let pb_style = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")?
        .progress_chars("#>-");

    if let Some(tsv) = &args.genome_list {
        let paths = get_genome_paths(tsv)?;
        let pb = ProgressBar::new(paths.len() as u64).with_style(pb_style);
        let results: Vec<String> = paths.par_iter().map(|(name, path)| {
            let res = analyze_path(name, path, &marker_kmer_map, k);
            pb.inc(1);
            res
        }).flatten().collect();
        pb.finish_with_message("Genome analysis complete");
        lines.extend(results);
    } else if let Some(mfasta) = &args.genome_fasta {
        let genomes = get_genomes_from_fasta(mfasta)?;
        let pb = ProgressBar::new(genomes.len() as u64).with_style(pb_style);
        let results: Vec<String> = genomes.par_iter().map(|(name, seq)| {
            let res = analyze_seq(name, seq, &marker_kmer_map, k);
            pb.inc(1);
            res
        }).flatten().collect();
        pb.finish_with_message("Genome analysis complete");
        lines.extend(results);
    }

    info!("▶ Writing output files...");
    let mut writer = BufWriter::new(File::create(&args.output)?);
    writeln!(writer, "genome\tk-mer\tk-mer_pos_genome\tsnp_pos_genome\tsnp_pos_reference\tlineage")?;
    for line in &lines {
        write!(writer, "{}", line)?;
    }
    info!("  Detailed marker report written to {}", &args.output);

    let summary_path = format!("{}_summary.tsv", args.output);
    let mut summary_writer = BufWriter::new(File::create(&summary_path)?);
    writeln!(summary_writer, "genome\tlineage_counts\tmajor_lineage")?;
    for (genome, map) in lineage_summary(&lines) {
        write_summary_line(&mut summary_writer, &genome, map.into_iter().collect())?;
    }
    info!("  Lineage summary written to {}", &summary_path);
    
    info!("✅ Process completed.");
    Ok(())
}
