//! This module handles the `predict` subcommand.
//!
//! Its primary responsibilities are:
//! 1. Loading a unified model bundle from a single file.
//! 2. Reading new genome sequences from a FASTA file.
//! 3. Transforming these sequences into feature vectors in parallel.
//! 4. Using the loaded model to predict the class for each sequence.
//! 5. Writing the predictions to an output file.

use crate::common::{kmerize, ModelBundle};
use crate::errors::{AppError, AppResult};
use clap::Parser;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, trace};
use smartcore::linalg::basic::arrays::Array; // <-- FIX: Import the Array trait
use smartcore::linalg::basic::matrix::DenseMatrix;
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
    let file_size = file.metadata()?.len();
    let pb = ProgressBar::new(file_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})",
            )?
            .progress_chars("#>-"),
    );

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
    debug!("Read {} records from {}", records.len(), path);
    Ok(records)
}

pub fn run(args: PredictArgs) -> AppResult<()> {
    if let Some(n) = args.threads {
        debug!("Setting number of threads to: {}", n);
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| AppError::Generic(format!("Failed to build thread pool: {}", e)))?;
    }

    info!("▶ Loading model bundle from {}", args.model);
    let model_file = File::open(&args.model)?;
    let mut model_decoder = GzDecoder::new(model_file);
    let bundle: ModelBundle = bincode::deserialize_from(&mut model_decoder)?;
    info!(
        "  Model bundle loaded. Using k-mer size: {}",
        bundle.config.kmer_size
    );
    debug!("Model contains {} features and {} classes.", bundle.vectorizer.feature_names.len(), bundle.label_encoder.int_to_label.len());


    info!("▶ Reading input FASTA file: {}", args.input);
    let records = read_fasta_for_prediction(&args.input)?;
    if records.is_empty() {
        return Err(AppError::NotEnoughData(
            "No records found in the input FASTA file.".to_string(),
        ));
    }

    info!("▶ Generating k-mers and transforming features...");
    let texts: Vec<String> = records
        .iter()
        .map(|(_, seq)| kmerize(seq, bundle.config.kmer_size))
        .collect();
    trace!("First k-merized text: '{}'", texts.first().unwrap_or(&"".to_string()));
    let x_data = bundle.vectorizer.transform(&texts);
    let x_matrix = DenseMatrix::from_2d_vec(&x_data)?;
    debug!("Transformed {} sequences into a feature matrix of size {}x{}", records.len(), x_matrix.shape().0, x_matrix.shape().1);


    info!("▶ Predicting lineages...");
    let y_pred = bundle.model.predict(&x_matrix)?;
    debug!("Generated {} predictions.", y_pred.len());

    let default_prediction = "Unknown".to_string();
    let predictions: Vec<String> = y_pred
        .iter()
        .map(|&class| {
            bundle
                .label_encoder
                .int_to_label
                .get(class)
                .cloned()
                .unwrap_or_else(|| default_prediction.clone())
        })
        .collect();
    trace!("First prediction: {}", predictions.first().unwrap_or(&"N/A".to_string()));


    let mut output_file = File::create(&args.output)?;
    writeln!(output_file, "Header\tPredicted_Lineage")?;
    for (i, (header, _)) in records.iter().enumerate() {
        let pred = predictions.get(i).unwrap_or(&default_prediction);
        writeln!(output_file, "{}\t{}", header, pred)?;
    }

    info!("✅ Predictions written to {}", args.output);
    Ok(())
}
