//! This module handles the `predict` subcommand.
//!
//! Its primary responsibilities are:
//! 1. Loading a unified model bundle from a single file.
//! 2. Reading new genome sequences from a FASTA file.
//! 3. Transforming these sequences into feature vectors in parallel.
//! 4. Using the loaded model to predict the class for each sequence.
//! 5. Writing the predictions to an output file.

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use smartcore::ensemble::random_forest_classifier::RandomForestClassifier;
use smartcore::linalg::basic::matrix::DenseMatrix;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[derive(Serialize, Deserialize, Debug)]
pub struct ModelConfig { pub kmer_size: usize, }
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelBundle {
    pub config: ModelConfig,
    pub vectorizer: CountVectorizer,
    pub label_encoder: LabelEncoder,
    pub model: RandomForestClassifier<f64, usize, DenseMatrix<f64>, Vec<usize>>,
}

#[derive(Parser, Debug)]
pub struct PredictArgs {
    #[arg(short = 'i', long)] pub input: String,
    #[arg(short = 'm', long)] pub model: String,
    #[arg(short = 'o', long)] pub output: String,
    #[arg(short = 't', long)] pub threads: Option<usize>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct CountVectorizer { pub vocabulary: HashMap<String, usize>, #[allow(dead_code)] pub feature_names: Vec<String> }
impl CountVectorizer {
    pub fn transform<T: AsRef<str> + Sync>(&self, texts: &[T]) -> Vec<Vec<f64>> {
        texts.par_iter().map(|text| {
            let mut counts = vec![0.0; self.vocabulary.len()];
            for token in text.as_ref().split_whitespace() { if let Some(&idx) = self.vocabulary.get(token) { counts[idx] += 1.0; } }
            counts
        }).collect()
    }
}
#[derive(Serialize, Deserialize, Debug)]
pub struct LabelEncoder { #[allow(dead_code)] pub label_to_int: HashMap<String, usize>, pub int_to_label: Vec<String> }

fn kmerize(sequence: &str, k: usize) -> String {
    if sequence.len() < k { return String::new(); }
    (0..=sequence.len() - k).map(|i| &sequence[i..i + k]).collect::<Vec<&str>>().join(" ")
}

fn read_fasta_for_prediction(path: &str) -> Result<Vec<(String, String)>> {
    let file = File::open(path)?;
    let file_size = file.metadata()?.len();
    let pb = ProgressBar::new(file_size);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")?
        .progress_chars("#>-"));

    let reader = BufReader::new(pb.wrap_read(file));
    let mut records = Vec::new();
    let mut current_header = String::new(); let mut current_seq = String::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_header.is_empty() { records.push((current_header.clone(), current_seq.clone())); }
            current_header = line.trim_start_matches('>').to_string(); current_seq.clear();
        } else { current_seq.push_str(line.trim()); }
    }
    if !current_header.is_empty() { records.push((current_header, current_seq)); }
    
    pb.finish_with_message("FASTA reading complete");
    Ok(records)
}

pub fn run(args: PredictArgs) -> Result<()> {
    if let Some(n) = args.threads { rayon::ThreadPoolBuilder::new().num_threads(n).build_global()?; }

    info!("▶ Loading model bundle from {}", args.model);
    let model_file = File::open(&args.model)?;
    let mut model_decoder = GzDecoder::new(model_file);
    let bundle: ModelBundle = bincode::deserialize_from(&mut model_decoder)
        .context("Failed to deserialize the model bundle.")?;
    info!("  Model bundle loaded. Using k-mer size: {}", bundle.config.kmer_size);

    info!("▶ Reading input FASTA file: {}", args.input);
    let records = read_fasta_for_prediction(&args.input)?;
    if records.is_empty() { return Err(anyhow!("No records found in the input FASTA file.")); }
    
    info!("▶ Generating k-mers and transforming features...");
    let texts: Vec<String> = records.par_iter().map(|(_, seq)| kmerize(seq, bundle.config.kmer_size)).collect();
    let x_data = bundle.vectorizer.transform(&texts);
    let x_matrix = DenseMatrix::from_2d_vec(&x_data)?;

    info!("▶ Predicting lineages...");
    let y_pred = bundle.model.predict(&x_matrix)?;
    
    let default_prediction = "Unknown".to_string();
    let predictions: Vec<String> = y_pred.iter().map(|&class| {
        bundle.label_encoder.int_to_label.get(class).cloned().unwrap_or(default_prediction.clone())
    }).collect();

    let mut output_file = File::create(&args.output)?;
    writeln!(output_file, "Header\tPredicted_Lineage")?;
    for (i, (header, _)) in records.iter().enumerate() {
        let pred = predictions.get(i).unwrap_or(&default_prediction);
        writeln!(output_file, "{}\t{}", header, pred)?;
    }
    
    info!("✅ Predictions written to {}", args.output);
    Ok(())
}
