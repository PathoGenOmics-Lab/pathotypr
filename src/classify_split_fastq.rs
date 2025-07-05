//! Implements the `SplitFastq` subcommand.
//!
//! This module now supports batch processing of multiple samples, either from
//! direct command-line input or from a TSV file. It automatically determines
//! sample names and generates separate reports for each, plus a final summary.

use crate::split_kmer;
use anyhow::{anyhow, Context, Result};
use clap::Parser;
use log::{info, warn};
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

#[derive(Parser, Debug)]
pub struct SplitFastqArgs {
    #[arg(short = 'i', long = "input", group = "input_method")] pub input: Vec<String>,
    #[arg(short = 'l', long = "input-list", group = "input_method")] pub input_list: Option<String>,
    #[arg(long)] pub paired: bool,
    #[arg(short = 'r', long, required = true)] pub reference: String,
    #[arg(short = 'm', long, required = true)] pub markers: String,
    #[arg(short = 't', long)] pub threads: Option<usize>,
    #[arg(short = 'o', long, default_value = "split")] pub output_prefix: String,
    #[arg(long, default_value_t = 10)] pub min_depth: u32,
    #[arg(long, default_value_t = 95)] pub min_alt_percent: u32,
}

fn derive_sample_name(path_str: &str) -> String {
    let path = Path::new(path_str);
    let file_stem = path.file_stem().unwrap_or_default().to_str().unwrap_or_default();
    file_stem.trim_end_matches(".R1").trim_end_matches(".R2").trim_end_matches("_R1").trim_end_matches("_R2").trim_end_matches("_1").trim_end_matches("_2").to_string()
}

fn read_sample_list(path: &str) -> Result<HashMap<String, Vec<String>>> {
    let mut samples = HashMap::new();
    let reader = fs::read_to_string(path)?;
    for line in reader.lines() {
        if line.trim().is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 { warn!("Skipping malformed line in sample list: {}", line); continue; }
        let sample_name = fields[0].to_string();
        let fastq_paths = fields[1..].iter().map(|s| s.to_string()).collect();
        samples.insert(sample_name, fastq_paths);
    }
    Ok(samples)
}

pub fn run(args: SplitFastqArgs) -> Result<()> {
    if args.input.is_empty() && args.input_list.is_none() {
        return Err(anyhow!("You must provide an input source: either --input or --input-list."));
    }

    rayon::ThreadPoolBuilder::new().num_threads(args.threads.unwrap_or(0)).build_global()?;

    info!("▶ Building split-k-mer table");
    let markers = split_kmer::build_markers(&args.reference, &args.markers)?;
    info!("  Loaded {} markers", markers.len());

    let mut samples_to_process: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(list_path) = &args.input_list {
        samples_to_process = read_sample_list(list_path)?;
    } else {
        for fq_path in &args.input {
            let sample_name = derive_sample_name(fq_path);
            samples_to_process.entry(sample_name).or_default().push(fq_path.clone());
        }
    }

    let mut all_summary_lines = Vec::new();
    for (sample_name, fastq_paths) in &samples_to_process {
        info!("▶ Processing sample: {}", sample_name);
        let counts = split_kmer::scan_fastq(fastq_paths, &markers)?;
        info!("  Finished scan for {}. Analyzing results...", sample_name);

        let detailed_output_path = format!("{}_{}_mutations.tsv", args.output_prefix, sample_name);
        let mut detailed_writer = fs::File::create(&detailed_output_path)?;
        writeln!(detailed_writer, "pos\tref\talt\tA\tC\tG\tT\tlineage")?;

        let mut lineage_counts: HashMap<String, usize> = HashMap::new();

        for (marker_id, base_counts) in counts.into_iter().enumerate() {
            let marker = &markers[marker_id];
            let coverage: u32 = base_counts.iter().sum();
            if coverage < args.min_depth { continue; }
            let alt_base_idx = match marker.alt_base {
                b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => continue,
            };
            let alt_count = base_counts[alt_base_idx];
            if alt_count * 100 < coverage * args.min_alt_percent { continue; }
            writeln!(detailed_writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", marker.pos + 1, marker.ref_base as char, marker.alt_base as char, base_counts[0], base_counts[1], base_counts[2], base_counts[3], marker.lineage)?;
            *lineage_counts.entry(marker.lineage.clone()).or_default() += 1;
        }
        info!("  Detailed report for {} written to {}", sample_name, detailed_output_path);

        let summary_line = if lineage_counts.is_empty() {
            format!("{}\t\t", sample_name)
        } else {
            let mut sorted_lineages: Vec<(String, usize)> = lineage_counts.into_iter().collect();
            sorted_lineages.sort_by(|a, b| b.1.cmp(&a.1));
            let list_str = sorted_lineages.iter().map(|(l, c)| format!("{}:{}", l, c)).collect::<Vec<_>>().join(",");
            let major_lineage = if sorted_lineages.len() == 1 || sorted_lineages[0].1 > sorted_lineages[1].1 {
                sorted_lineages[0].0.clone()
            } else {
                let top_count = sorted_lineages[0].1;
                sorted_lineages.iter().filter(|(_, c)| *c == top_count).map(|(l, _)| l.clone()).collect::<Vec<_>>().join(",")
            };
            format!("{}\t{}\t{}", sample_name, list_str, major_lineage)
        };
        all_summary_lines.push(summary_line);
    }

    let summary_path = format!("{}_summary.tsv", args.output_prefix);
    info!("▶ Writing final summary for all samples to {}", summary_path);
    let mut summary_writer = fs::File::create(&summary_path)?;
    writeln!(summary_writer, "genome\tlineage:count\tmajor_lineage")?;
    for line in all_summary_lines { writeln!(summary_writer, "{}", line)?; }
    
    info!("✅ Process completed.");
    Ok(())
}
