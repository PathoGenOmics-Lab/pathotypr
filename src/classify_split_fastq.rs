//! Implements the `SplitFastq` subcommand.
//!
//! This module now supports batch processing of multiple samples, either from
//! direct command-line input or from a TSV file. It automatically determines
//! sample names and generates separate reports for each, plus a final summary.
//! It uses a dynamic k-mer engine to detect SNPs, MNVs, and Indels.

use crate::split_kmer;
use anyhow::{anyhow, Result};
use clap::Parser;
use log::{info, warn};
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

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

    /// If specified, treats input files as paired-end reads, grouped in pairs.
    #[arg(long)]
    pub paired: bool,

    /// Reference FASTA file used to define the markers.
    #[arg(short = 'r', long, required = true)]
    pub reference: String,

    /// Path to a TSV file defining markers.
    /// Format: position\tREF\tALT\tmarker...
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
}

/// Derives a clean sample name from a file path by taking the part before the first delimiter.
fn derive_sample_name(path_str: &str) -> String {
    let file_stem = Path::new(path_str)
        .file_stem()
        .unwrap_or_default()
        .to_str()
        .unwrap_or_default();
    file_stem
        .split(|c| c == '_' || c == '.')
        .next()
        .unwrap_or(file_stem)
        .to_string()
}

/// Reads a TSV file to get a map of sample names to their FASTQ file paths.
fn read_sample_list(path: &str) -> Result<HashMap<String, Vec<String>>> {
    let mut samples = HashMap::new();
    let reader = fs::read_to_string(path)?;
    for line in reader.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            warn!("Skipping malformed line in sample list: {}", line);
            continue;
        }
        let sample_name = fields[0].to_string();
        let fastq_paths = fields[1..].iter().map(|s| s.to_string()).collect();
        samples.insert(sample_name, fastq_paths);
    }
    Ok(samples)
}

pub fn run(args: SplitFastqArgs) -> Result<()> {
    if args.input.is_empty() && args.input_list.is_none() {
        return Err(anyhow!(
            "You must provide an input source: either --input or --input-list."
        ));
    }

    ThreadPoolBuilder::new()
        .num_threads(args.threads.unwrap_or(0))
        .build_global()?;

    info!("▶ Building dynamic marker database...");
    let markers = split_kmer::build_markers(&args.reference, &args.markers)?;
    info!("  Successfully generated {} dynamic markers.", markers.len());

    let mut samples_to_process: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(list_path) = &args.input_list {
        samples_to_process = read_sample_list(list_path)?;
    } else if args.paired {
        if args.input.len() % 2 != 0 {
            return Err(anyhow!(
                "--paired requires an even number of input files. Found {}.",
                args.input.len()
            ));
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
            samples_to_process.insert(sample_name, vec![r1_path.clone(), r2_path.clone()]);
        }
    } else {
        for fq_path in &args.input {
            let sample_name = derive_sample_name(fq_path);
            samples_to_process
                .entry(sample_name)
                .or_default()
                .push(fq_path.clone());
        }
    }

    info!("Found {} sample(s) to process.", samples_to_process.len());
    let mut all_summary_lines = Vec::new();

    for (sample_name, fastq_paths) in &samples_to_process {
        info!("▶ Processing sample: {}", sample_name);
        let counts = split_kmer::scan_fastq(fastq_paths, &markers)?;
        info!("  Finished scan for {}. Analyzing results...", sample_name);

        let detailed_output_path =
            format!("{}_{}_mutations.tsv", args.output_prefix, sample_name);
        let mut detailed_writer = fs::File::create(&detailed_output_path)?;
        
        let header = "pos\tref_allele\talt_allele\tref_count\talt_count\talt_fraction\tmarker\textra_annotations...";
        writeln!(detailed_writer, "{}", header)?;

        let mut lineage_counts: HashMap<String, usize> = HashMap::new();
        for (marker_id, &[ref_count, alt_count]) in counts.iter().enumerate() {
            let coverage = ref_count + alt_count;
            if coverage < args.min_depth {
                continue;
            }
            let alt_fraction = if coverage > 0 {
                (alt_count as f32 / coverage as f32) * 100.0
            } else {
                0.0
            };
            if alt_fraction < args.min_alt_percent as f32 {
                continue;
            }
            
            let marker = &markers[marker_id];
            
            let mut output_line = format!(
                "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
                marker.pos + 1,
                marker.ref_allele,
                marker.alt_allele,
                ref_count,
                alt_count,
                alt_fraction,
                marker.lineage
            );

            if !marker.annotations.is_empty() {
                output_line.push('\t');
                output_line.push_str(&marker.annotations.join("\t"));
            }
            writeln!(detailed_writer, "{}", output_line)?;

            *lineage_counts
                .entry(marker.lineage.clone())
                .or_default() += 1;
        }
        info!(
            "  Detailed report for {} written to {}",
            sample_name, detailed_output_path
        );

        let summary_line = if lineage_counts.is_empty() {
            format!("{}\t\t", sample_name)
        } else {
            let mut sorted_lineages: Vec<(String, usize)> = lineage_counts.into_iter().collect();
            sorted_lineages.sort_by(|a, b| b.1.cmp(&a.1));
            let list_str = sorted_lineages
                .iter()
                .map(|(l, c)| format!("{}:{}", l, c))
                .collect::<Vec<_>>()
                .join(",");
            let major_lineage =
                if sorted_lineages.len() == 1 || sorted_lineages[0].1 > sorted_lineages[1].1 {
                    sorted_lineages[0].0.clone()
                } else {
                    let top_count = sorted_lineages[0].1;
                    sorted_lineages
                        .iter()
                        .filter(|(_, c)| *c == top_count)
                        .map(|(l, _)| l.clone())
                        .collect::<Vec<_>>()
                        .join(",")
                };
            format!("{}\t{}\t{}", sample_name, list_str, major_lineage)
        };
        all_summary_lines.push(summary_line);
    }

    let summary_path = format!("{}_summary.tsv", args.output_prefix);
    info!(
        "▶ Writing final summary for all samples to {}",
        summary_path
    );
    let mut summary_writer = fs::File::create(&summary_path)?;
    writeln!(summary_writer, "genome\tlineage:count\tmajor_lineage")?;
    for line in all_summary_lines {
        writeln!(summary_writer, "{}", line)?;
    }

    info!("✅ Process completed.");
    Ok(())
}
