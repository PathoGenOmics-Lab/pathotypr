//! This module handles the `predict` subcommand.
//!
//! Its primary responsibilities are:
//! 1. Loading a model bundle containing an ensemble of decision trees.
//! 2. Predicting the class for new sequences by aggregating votes from all trees.
//! 3. Calculating confidence metrics and reporting votes for top contenders.

use crate::common::{kmerize, ModelBundle};
use crate::errors::{AppError, AppResult};
use clap::Parser;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, trace};
use smartcore::linalg::basic::matrix::DenseMatrix;
// Add traits needed for matrix operations
use smartcore::linalg::basic::arrays::{Array, Array2};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Parser, Debug)]
pub struct PredictArgs {
    /// Path to the input FASTA file containing sequences to classify.
    #[arg(short = 'i', long)]
    pub input: String,
    /// Path to the unified model file created by the train command.
    #[arg(short = 'm', long)]
    pub model: String,
    /// Path for the output file where predictions will be written in TSV format.
    #[arg(short = 'o', long)]
    pub output: String,
    /// Number of CPU threads to use.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
}

fn read_fasta_for_prediction(path: &str) -> AppResult<Vec<(String, String)>> {
    let file = File::open(path)?;
    let pb = ProgressBar::new(file.metadata()?.len());
    pb.set_style(ProgressStyle::default_bar().template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")?.progress_chars("#>-"));
    let reader = BufReader::new(pb.wrap_read(file));
    let mut records = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                records.push((current_header.clone(), current_seq.clone()));
            }
            current_header = line.trim_start_matches('>').to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_header.is_empty() {
        records.push((current_header, current_seq));
    }
    pb.finish_with_message("FASTA reading complete");
    Ok(records)
}

pub fn run(args: PredictArgs) -> AppResult<()> {
    if let Some(n) = args.threads {
        debug!("Setting number of threads to: {}", n);
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().map_err(|e| AppError::Generic(format!("Failed to build thread pool: {}", e)))?;
    }

    info!("? Loading model bundle from {}", args.model);
    let model_file = File::open(&args.model)?;
    let mut model_decoder = GzDecoder::new(model_file);
    let bundle: ModelBundle = bincode::deserialize_from(&mut model_decoder)?;
    info!("  Model bundle loaded. Using k-mer size: {} and {} trees.", bundle.config.kmer_size, bundle.config.n_trees);

    info!("? Reading input FASTA file: {}", args.input);
    let records = read_fasta_for_prediction(&args.input)?;
    if records.is_empty() { return Err(AppError::NotEnoughData("No records found in the input FASTA file.".to_string())); }

    info!("? Generating k-mers and transforming features...");
    let texts: Vec<String> = records.iter().map(|(_, seq)| kmerize(seq, bundle.config.kmer_size)).collect();
    let x_data = bundle.vectorizer.transform(&texts);
    let x_matrix = DenseMatrix::from_2d_vec(&x_data)?;
    debug!("Transformed {} sequences into a feature matrix.", records.len());

    info!("? Predicting lineages by aggregating tree votes...");
    // The vector will now store: (label, confidence, margin, other_votes_string)
    let mut predictions_with_metrics = Vec::new();
    let total_trees = bundle.trees.len() as f64;
    let default_prediction = "Unknown".to_string();

    for i in 0..x_matrix.shape().0 {
        let row_vec: Vec<f64> = x_matrix.get_row(i).iterator(0).copied().collect();
        let x_row = DenseMatrix::from_2d_vec(&vec![row_vec])?;
        let mut votes: HashMap<usize, u16> = HashMap::new();

        for tree in &bundle.trees {
            let prediction = tree.predict(&x_row)?[0];
            *votes.entry(prediction).or_insert(0) += 1;
        }

        let mut sorted_votes: Vec<(&usize, &u16)> = votes.iter().collect();
        sorted_votes.sort_by(|a, b| b.1.cmp(a.1));

        if let Some(winner) = sorted_votes.get(0) {
            let best_class_idx = *winner.0;
            let winner_vote_count = *winner.1;
            
            let second_vote_count = sorted_votes.get(1).map_or(0, |&(_, count)| *count);

            let confidence = winner_vote_count as f64 / total_trees;
            let margin = (winner_vote_count - second_vote_count) as f64 / total_trees;

            let predicted_label = bundle.label_encoder.int_to_label
                .get(best_class_idx)
                .cloned()
                .unwrap_or_else(|| default_prediction.clone());
            
            // --- NEW LOGIC: Format other votes string ---
            let other_votes_str = sorted_votes
                .iter()
                .skip(1) // Skip the winner
                .take(3) // Take the next top 3 contenders
                .map(|(class_idx, vote_count)| {
                    let label = bundle.label_encoder.int_to_label.get(**class_idx).unwrap_or(&default_prediction);
                    let percentage = **vote_count as f64 / total_trees;
                    format!("{}:{:.2}", label, percentage)
                })
                .collect::<Vec<String>>()
                .join(",");
            
            predictions_with_metrics.push((predicted_label, confidence, margin, other_votes_str));
        } else {
            predictions_with_metrics.push((default_prediction.clone(), 0.0, 0.0, String::new()));
        }
    }
    trace!("First prediction with metrics: {:?}", predictions_with_metrics.first());

    let mut output_file = File::create(&args.output)?;
    // Add the new column to the header
    writeln!(output_file, "Header\tPredicted_Lineage\tConfidence\tConfidence_Margin\tOther_Votes")?;
    for (i, (header, _)) in records.iter().enumerate() {
        if let Some((pred, conf, margin, others)) = predictions_with_metrics.get(i) {
            // Write the new `others` string to the output
            writeln!(output_file, "{}\t{}\t{:.4}\t{:.4}\t{}", header, pred, conf, margin, others)?;
        } else {
            writeln!(output_file, "{}\t{}\t{:.4}\t{:.4}\t", header, default_prediction, 0.0, 0.0)?;
        }
    }
    
    info!("? Predictions written to {}", args.output);
    Ok(())
}
