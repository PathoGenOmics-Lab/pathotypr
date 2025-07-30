//! Implements the `SplitFastq` subcommand.
//!
//! This module now supports batch processing of multiple samples, either from
//! direct command-line input or from a TSV file. It automatically determines
//! sample names and generates separate reports for each, plus a final summary.
//! It uses a dynamic k-mer engine to detect SNPs, MNVs, and Indels.
//! MODIFIED: Includes strictly corrected logic for nested lineage classification.

use crate::errors::{AppError, AppResult};
use crate::split_kmer;
use clap::Parser;
use log::{debug, info, trace, warn};
use rayon::ThreadPoolBuilder;
use std::collections::{HashMap, HashSet};
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
    /// Format: position\tREF\tALT\tlevel1\tlevel2...
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

    /// [NEW] Enable nested lineage classification logic.
    /// This requires markers to be defined with multiple columns for each level.
    #[arg(long)]
    pub nested_classification: bool,
}

/// Derives a clean sample name from a file path by taking the part before the first delimiter.
fn derive_sample_name(path_str: &str) -> String {
    let file_stem = Path::new(path_str)
        .file_stem()
        .unwrap_or_default()
        .to_str()
        .unwrap_or_default();
    file_stem
        .split(|c| c == '_' || c == ';')
        .next()
        .unwrap_or(file_stem)
        .to_string()
}

/// Reads a TSV file to get a map of sample names to their FASTQ file paths.
fn read_sample_list(path: &str) -> AppResult<HashMap<String, Vec<String>>> {
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
    debug!("Read {} samples from list {}", samples.len(), path);
    Ok(samples)
}

/// [REWRITTEN] Final logic for nested classification.
/// 1. Finds the deepest valid hierarchical path.
/// 2. Finds the most abundant lineage overall.
/// 3. If the most abundant is from a different branch and has more support, it wins.
/// 4. Otherwise, the deepest path wins.
fn get_final_lineage_call(lineage_counts: &HashMap<String, usize>) -> String {
    if lineage_counts.is_empty() {
        return "Unclassified".to_string();
    }

    let supported_lineages: HashSet<String> = lineage_counts.keys().cloned().collect();
    let mut valid_candidates = Vec::new();

    // 1. Identify all lineages with a valid, fully supported path.
    for candidate in supported_lineages.iter() {
        let mut is_path_valid = true;
        if candidate.contains(';') {
            let components: Vec<&str> = candidate.split(';').collect();
            for i in 1..components.len() {
                let parent_path = components[0..i].join(";");
                if !supported_lineages.contains(&parent_path) {
                    is_path_valid = false;
                    break;
                }
            }
        }
        if is_path_valid {
            valid_candidates.push(candidate.clone());
        }
    }

    // If no valid paths exist, fallback to the most abundant lineage.
    if valid_candidates.is_empty() {
        return lineage_counts
            .iter()
            .max_by_key(|&(_, count)| count)
            .map(|(lineage, _)| lineage.clone())
            .unwrap_or_else(|| "Unclassified".to_string());
    }

    // 2. From the valid candidates, find the deepest one.
    // If there's a tie in depth, the higher SNP count for that specific deep lineage wins.
    let best_deep_lineage = valid_candidates
        .iter()
        .max_by(|a, b| {
            let depth_a = a.split(';').count();
            let depth_b = b.split(';').count();
            depth_a.cmp(&depth_b)
                   .then_with(|| lineage_counts[*a].cmp(&lineage_counts[*b]))
        })
        .unwrap(); // Safe because valid_candidates is not empty.

    // 3. Find the most abundant lineage overall.
    let most_abundant_lineage = lineage_counts
        .iter()
        .max_by_key(|&(_, count)| count)
        .map(|(lineage, _)| lineage)
        .unwrap(); // Safe because lineage_counts is not empty.

    // 4. The final decision logic.
    let best_deep_count = lineage_counts[best_deep_lineage];
    let most_abundant_count = lineage_counts[most_abundant_lineage];

    // Check if the most abundant lineage is from a different branch than the deepest one.
    // A shared branch means one starts with the other.
    let shares_branch = best_deep_lineage.starts_with(most_abundant_lineage) ||
                        most_abundant_lineage.starts_with(best_deep_lineage);

    if !shares_branch && most_abundant_count > best_deep_count {
        // The most abundant lineage is from a different branch and has more support, so it wins.
        most_abundant_lineage.clone()
    } else {
        // Otherwise, prioritize the deepest valid path.
        best_deep_lineage.clone()
    }
}


pub fn run(args: SplitFastqArgs) -> AppResult<()> {
    if args.input.is_empty() && args.input_list.is_none() {
        return Err(AppError::Generic(
            "You must provide an input source: either --input or --input-list.".to_string(),
        ));
    }

    if let Some(n) = args.threads {
        debug!("Setting number of threads to: {}", n);
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| AppError::Generic(format!("Failed to build thread pool: {}", e)))?;
    }

    info!("▶ Building dynamic marker database...");
    let markers = split_kmer::build_markers(&args.reference, &args.markers)
        .map_err(|e| AppError::Generic(format!("Failed to build markers: {}", e)))?;
    info!("  Successfully generated {} dynamic markers.", markers.len());
    trace!("First marker generated: {:?}", markers.first());

    let mut samples_to_process: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(list_path) = &args.input_list {
        samples_to_process = read_sample_list(list_path)?;
    } else if args.paired {
        if args.input.len() % 2 != 0 {
            return Err(AppError::Generic(format!(
                "--paired requires an even number of input files. Found {}.",
                args.input.len()
            )));
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
    debug!("Samples to process: {:?}", samples_to_process.keys());
    let mut all_summary_lines = Vec::new();

    for (sample_name, fastq_paths) in &samples_to_process {
        info!("▶ Processing sample: {}", sample_name);
        let counts = split_kmer::scan_fastq(fastq_paths, &markers)
            .map_err(|e| AppError::Generic(format!("Failed to scan FASTQ for {}: {}", sample_name, e)))?;
        info!("  Finished scan for {}. Analyzing results...", sample_name);
        trace!("Raw counts for sample {}: {:?}", sample_name, counts);

        let detailed_output_path =
            format!("{}_{}_mutations.tsv", args.output_prefix, sample_name);
        let mut detailed_writer = fs::File::create(&detailed_output_path)?;
        
        let header = "pos\tref_allele\talt_allele\tref_count\talt_count\talt_fraction\tlineage_path\textra_annotations...";
        writeln!(detailed_writer, "{}", header)?;

        // This map will store counts for ALL levels defined by each found marker.
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
            
            if alt_fraction >= args.min_alt_percent as f32 {
                let marker = &markers[marker_id];
                
                let lineage_path_string = marker.lineages.join(";");

                // Write detailed report line
                let mut output_line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
                    marker.pos + 1,
                    marker.ref_allele,
                    marker.alt_allele,
                    ref_count,
                    alt_count,
                    alt_fraction,
                    lineage_path_string
                );

                if !marker.annotations.is_empty() {
                    output_line.push('\t');
                    output_line.push_str(&marker.annotations.join("\t"));
                }
                writeln!(detailed_writer, "{}", output_line)?;

                // Add a count for each level in this marker's hierarchy to ensure parents are counted.
                let mut current_path_part = String::new();
                for (i, component) in marker.lineages.iter().enumerate() {
                     current_path_part = if i == 0 {
                        component.to_string()
                    } else {
                        format!("{};{}", current_path_part, component)
                    };
                    *lineage_counts.entry(current_path_part.clone()).or_default() += 1;
                }
            }
        }
        
        info!(
            "  Detailed report for {} written to {}",
            sample_name, detailed_output_path
        );

        // Generate summary line
        let mut sorted_lineages: Vec<(String, usize)> = lineage_counts.iter().map(|(k,v)| (k.clone(), *v)).collect();
        sorted_lineages.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        
        let list_str = sorted_lineages
            .iter()
            .map(|(l, c)| format!("{}:{}", l, c))
            .collect::<Vec<_>>()
            .join(" ");

        let major_lineage = if args.nested_classification {
            // Use the new nested logic if the flag is active
            get_final_lineage_call(&lineage_counts)
        } else {
            // Original logic for majority lineage call
            if sorted_lineages.is_empty() {
                "Unclassified".to_string()
            } else {
                sorted_lineages[0].0.clone()
            }
        };

        let summary_line = format!("{}\t{}\t{}", sample_name, list_str, major_lineage);
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
