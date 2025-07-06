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

use anyhow::{anyhow, Result};
use bio::io::fasta;
use clap::Parser;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{error, info, warn};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

#[derive(Parser, Debug)]
pub struct Args {
    /// Path to a TSV file defining markers.
    /// Format: position\tREF\tALT\tmarker...
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

    /// The size of the diagnostic k-mers to generate and search for.
    #[arg(short = 'k', long, default_value_t = 31)]
    pub kmer_size: usize,

    /// Number of CPU threads to use for parallel processing.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
}

#[derive(Debug, Clone)]
struct MarkerInfo {
    pos: usize,
    lineage: String,
    annotations: Vec<String>,
}

/// Dynamically builds a map of diagnostic k-mers for alternate alleles.
fn build_diagnostic_kmers(marker_path: &str, ref_seq: &str, k: usize) -> Result<HashMap<String, MarkerInfo>> {
    info!("  Parsing marker file and generating diagnostic k-mers...");
    let mut marker_map = HashMap::new();
    let ref_bytes = ref_seq.as_bytes();
    let ref_len = ref_bytes.len();

    let file = File::open(marker_path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line_str = line?;
        if line_str.trim().is_empty() || line_str.starts_with('#') { continue; }
        let fields: Vec<&str> = line_str.split('\t').collect();
        if fields.len() < 4 { continue; }

        let pos0: usize = match fields[0].parse::<usize>() {
            Ok(p) => p.saturating_sub(1),
            Err(_) => { warn!("Skipping marker with invalid position: {}", fields[0]); continue; }
        };
        let ref_allele = fields[1];
        let alt_allele = fields[2];
        let lineage = fields[3].to_string();
        let annotations = if fields.len() > 4 { fields[4..].iter().map(|s| s.to_string()).collect() } else { Vec::new() };
        
        let alt_len = alt_allele.len();
        let ref_len_allele = ref_allele.len();
        
        if ref_len_allele.max(alt_len) < k {
             // --- Logic for Small Variants ---
            let mut alt_kmer = Vec::with_capacity(k);
            let flank_len = (k - alt_len) / 2;
            let kmer_start_pos = pos0.saturating_sub(flank_len);
            
            if kmer_start_pos + k > ref_len { continue; }

            alt_kmer.extend_from_slice(&ref_bytes[kmer_start_pos..pos0]);
            alt_kmer.extend_from_slice(alt_allele.as_bytes());
            
            let right_flank_start = pos0 + ref_len_allele;
            let needed_right_flank = k.saturating_sub(alt_kmer.len());

            if right_flank_start + needed_right_flank > ref_len { continue; }

            alt_kmer.extend_from_slice(&ref_bytes[right_flank_start .. right_flank_start + needed_right_flank]);

            if let Ok(kmer_str) = String::from_utf8(alt_kmer) {
                marker_map.insert(kmer_str, MarkerInfo { pos: pos0 + 1, lineage, annotations });
            }
        } else {
            // --- NEW LOGIC for Large Variants (SVs) ---
            warn!("Handling large variant at pos {} as a structural variant.", pos0 + 1);
            let k_half = k / 2;
            if ref_len_allele > alt_len { // Large Deletion or complex delins
                let mut junction_kmer = Vec::with_capacity(k);
                if pos0 < k_half { continue; }
                let right_flank_start = pos0 + ref_len_allele;
                if right_flank_start + k_half > ref_len { continue; }
                
                junction_kmer.extend_from_slice(&ref_bytes[pos0 - k_half .. pos0]);
                junction_kmer.extend_from_slice(alt_allele.as_bytes());
                let needed_right_len = k.saturating_sub(junction_kmer.len());

                if right_flank_start + needed_right_len > ref_len { continue; }
                junction_kmer.extend_from_slice(&ref_bytes[right_flank_start .. right_flank_start + needed_right_len]);
                
                if let Ok(kmer_str) = String::from_utf8(junction_kmer) {
                    marker_map.insert(kmer_str, MarkerInfo { pos: pos0 + 1, lineage, annotations });
                }
            } else { // Large Insertion
                if pos0 >= k_half { // Left Junction
                    let mut kmer = Vec::with_capacity(k);
                    kmer.extend_from_slice(&ref_bytes[pos0 - k_half .. pos0]);
                    kmer.extend_from_slice(&alt_allele.as_bytes()[0..k.saturating_sub(k_half)]);
                    if let Ok(kmer_str) = String::from_utf8(kmer) {
                        marker_map.insert(kmer_str, MarkerInfo { pos: pos0 + 1, lineage: lineage.clone(), annotations: annotations.clone() });
                    }
                }
                let right_flank_start = pos0 + ref_len_allele;
                if right_flank_start + k_half <= ref_len { // Right Junction
                    let mut kmer = Vec::with_capacity(k);
                    kmer.extend_from_slice(&alt_allele.as_bytes()[alt_len - k_half ..]);
                    kmer.extend_from_slice(&ref_bytes[right_flank_start .. right_flank_start + k.saturating_sub(k_half)]);
                    if let Ok(kmer_str) = String::from_utf8(kmer) {
                        marker_map.insert(kmer_str, MarkerInfo { pos: pos0 + 1, lineage, annotations });
                    }
                }
            }
        }
    }
    Ok(marker_map)
}

/// Generates all k-mers from a sequence and stores them in a HashMap for quick lookup.
fn generate_kmers(seq: &str, k: usize) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    if seq.len() < k { return m; }
    for i in 0..=seq.len() - k {
        m.insert(seq[i..i + k].to_string(), i);
    }
    m
}

/// Finds which of the provided markers are present in a given genome sequence.
fn find_markers<'a>(genome_seq: &str, markers: &'a HashMap<String, MarkerInfo>, k: usize) -> HashMap<String, (usize, &'a MarkerInfo)> {
    let genome_kmers = generate_kmers(genome_seq, k);
    let mut found_markers = HashMap::new();
    for (marker_kmer, marker_info) in markers {
        if let Some(&genome_pos) = genome_kmers.get(marker_kmer) {
            found_markers.insert(marker_kmer.clone(), (genome_pos, marker_info));
        }
    }
    found_markers
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
        let rec = rec?; if rec.len() < 2 { continue; }
        let name = rec[0].to_string(); let path = rec[1].to_string();
        if !Path::new(&path).exists() { error!("Missing FASTA for genome {}: {}", name, path); continue; }
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
fn analyze_path(name: &str, path: &str, markers: &HashMap<String, MarkerInfo>, k: usize) -> Vec<String> {
    let seq = match get_ref(path) {
        Ok(s) => s,
        Err(e) => { error!("Could not read {}: {}", path, e); return Vec::new(); }
    };
    analyze_seq(name, &seq, markers, k)
}

/// Analyzes a single genome sequence for markers and formats the output lines.
fn analyze_seq(name: &str, seq: &str, markers: &HashMap<String, MarkerInfo>, k: usize) -> Vec<String> {
    let mut v = Vec::new();
    let found_markers = find_markers(seq, markers, k);
    if found_markers.is_empty() {
        v.push(format!("{}\n", name));
    } else {
        for (kmer, (genome_pos, marker_info)) in found_markers {
            let snp_pos_in_genome = genome_pos + k / 2;
            let mut line = format!("{}\t{}\t{}\t{}\t{}\t{}", name, kmer, genome_pos, snp_pos_in_genome, marker_info.pos, marker_info.lineage);
            if !marker_info.annotations.is_empty() {
                line.push('\t');
                line.push_str(&marker_info.annotations.join("\t"));
            }
            v.push(line);
        }
    }
    v
}

/// Summarizes marker counts per lineage for each genome.
fn lineage_summary(lines: &[String]) -> HashMap<String, HashMap<String, usize>> {
    let mut summary_map = HashMap::new();
    for line in lines {
        let fields: Vec<&str> = line.trim().split('\t').collect();
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
    let major_lineage = if lineage_counts.is_empty() { "".to_string() }
    else if lineage_counts.len() == 1 || lineage_counts[0].1 > lineage_counts[1].1 { lineage_counts[0].0.clone() }
    else {
        let top_count = lineage_counts[0].1;
        lineage_counts.iter().filter(|(_, count)| *count == top_count).map(|(l, _)| l.clone()).collect::<Vec<_>>().join(",")
    };
    writeln!(w, "{}\t{}\t{}", genome, list, major_lineage)?;
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    rayon::ThreadPoolBuilder::new().num_threads(args.threads.unwrap_or(0)).build_global()?;
    if args.genome_list.is_none() && args.genome_fasta.is_none() {
        return Err(anyhow!("You must provide an input source: either --genome-list or --genome-fasta."));
    }

    info!("▶ Preparing marker k-mers...");
    let k = args.kmer_size;
    let ref_seq = get_ref(&args.reference)?;
    let marker_kmer_map = build_diagnostic_kmers(&args.markers, &ref_seq, k)?;
    info!("  Loaded {} unique dynamic marker k-mers.", marker_kmer_map.len());

    info!("▶ Analyzing genomes...");
    let mut lines = Vec::<String>::new();
    let pb_style = ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")?.progress_chars("#>-");

    if let Some(tsv) = &args.genome_list {
        let paths = get_genome_paths(tsv)?;
        let pb = ProgressBar::new(paths.len() as u64).with_style(pb_style);
        let results: Vec<String> = paths.par_iter().map(|(name, path)| {
            let res = analyze_path(name, path, &marker_kmer_map, k); pb.inc(1); res
        }).flatten().collect();
        pb.finish_with_message("Genome analysis complete");
        lines.extend(results);
    } else if let Some(mfasta) = &args.genome_fasta {
        let genomes = get_genomes_from_fasta(mfasta)?;
        let pb = ProgressBar::new(genomes.len() as u64).with_style(pb_style);
        let results: Vec<String> = genomes.par_iter().map(|(name, seq)| {
            let res = analyze_seq(name, seq, &marker_kmer_map, k); pb.inc(1); res
        }).flatten().collect();
        pb.finish_with_message("Genome analysis complete");
        lines.extend(results);
    }

    info!("▶ Writing output files...");
    let mut writer = BufWriter::new(File::create(&args.output)?);
    let header = "genome\tk-mer\tk-mer_pos_genome\tsnp_pos_genome\tsnp_pos_reference\tmarker\textra_annotations...";
    writeln!(writer, "{}", header)?;
    for line in &lines { writeln!(writer, "{}", line)?; }
    info!("  Detailed marker report written to {}", &args.output);

    let summary_path = format!("{}_summary.tsv", Path::new(&args.output).file_stem().unwrap().to_str().unwrap());
    let mut summary_writer = BufWriter::new(File::create(&summary_path)?);
    writeln!(summary_writer, "genome\tlineage_counts\tmajor_lineage")?;
    for (genome, map) in lineage_summary(&lines) {
        write_summary_line(&mut summary_writer, &genome, map.into_iter().collect())?;
    }
    info!("  Lineage summary written to {}", &summary_path);
    
    info!("✅ Process completed.");
    Ok(())
}
