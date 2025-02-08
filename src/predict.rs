// src/predict.rs

use clap::Parser;
use chrono::Local;
use flate2::read::GzDecoder;
use indicatif::ProgressBar;
use serde::{Deserialize, Serialize};
use smartcore::ensemble::random_forest_classifier::RandomForestClassifier;
use smartcore::linalg::basic::matrix::DenseMatrix;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

const DEFAULT_KMER_SIZE: usize = 6;

/// Converts a genomic sequence into overlapping k-mers.
fn kmerize(sequence: &str, k: usize) -> String {
    if sequence.len() < k {
        return String::new();
    }
    (0..=sequence.len() - k)
        .map(|i| &sequence[i..i + k])
        .collect::<Vec<&str>>()
        .join(" ")
}

/// A simple count vectorizer that splits texts on whitespace.
#[derive(Serialize, Deserialize, Debug)]
pub struct CountVectorizer {
    pub vocabulary: HashMap<String, usize>,
    pub feature_names: Vec<String>,
}

impl CountVectorizer {
    pub fn new() -> Self {
        Self {
            vocabulary: HashMap::new(),
            feature_names: Vec::new(),
        }
    }
    pub fn fit<T: AsRef<str>>(&mut self, texts: &[T]) {
        let mut freq: HashMap<String, usize> = HashMap::new();
        for text in texts {
            for token in text.as_ref().split_whitespace() {
                *freq.entry(token.to_string()).or_insert(0) += 1;
            }
        }
        let mut freq_vec: Vec<(String, usize)> = freq.into_iter().collect();
        freq_vec.sort_by(|a, b| b.1.cmp(&a.1));
        self.vocabulary = freq_vec
            .iter()
            .enumerate()
            .map(|(i, (token, _))| (token.clone(), i))
            .collect();
        self.feature_names = freq_vec.into_iter().map(|(token, _)| token).collect();
    }
    pub fn transform<T: AsRef<str> + Sync>(&self, texts: &[T]) -> Vec<Vec<f64>> {
        texts
            .iter()
            .map(|text| {
                let n_features = self.vocabulary.len();
                let mut counts = vec![0.0; n_features];
                for token in text.as_ref().split_whitespace() {
                    if let Some(&idx) = self.vocabulary.get(token) {
                        counts[idx] += 1.0;
                    }
                }
                counts
            })
            .collect()
    }
}

/// Label encoder that maps labels (strings) to numeric values.
#[derive(Serialize, Deserialize, Debug)]
pub struct LabelEncoder {
    pub label_to_int: HashMap<String, usize>,
    pub int_to_label: Vec<String>,
}

impl LabelEncoder {
    pub fn new() -> Self {
        Self {
            label_to_int: HashMap::new(),
            int_to_label: Vec::new(),
        }
    }
    pub fn fit<T: AsRef<str>>(&mut self, labels: &[T]) {
        for label in labels {
            let label_str = label.as_ref();
            if !self.label_to_int.contains_key(label_str) {
                let index = self.int_to_label.len();
                self.label_to_int.insert(label_str.to_string(), index);
                self.int_to_label.push(label_str.to_string());
            }
        }
    }
    pub fn transform<T: AsRef<str>>(&self, labels: &[T]) -> Vec<usize> {
        labels
            .iter()
            .map(|label| *self.label_to_int.get(label.as_ref()).unwrap())
            .collect()
    }
}

/// Command-line arguments for the predict subcommand.
#[derive(Parser, Debug)]
pub struct PredictArgs {
    /// Input FASTA file (multi-FASTA; header in format "Lineage_sequenceID")
    #[arg(long)]
    pub fasta: String,
    /// Base name of the saved model (expects files: <model_base>_rf_model.bin.gz, etc.)
    #[arg(long)]
    pub model_base: String,
    /// Output file where predictions will be written.
    #[arg(long)]
    pub output: String,
    /// k-mer size (default is 6)
    #[arg(long, default_value_t = DEFAULT_KMER_SIZE)]
    pub kmer_size: usize,
}

fn read_fasta_for_prediction(path: &str) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let pb = ProgressBar::new_spinner();
    pb.set_message("Processing FASTA records for prediction...");
    let mut records = Vec::new();
    let mut current_header = String::new();
    let mut current_seq = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                records.push((current_header.clone(), current_seq.clone()));
                pb.inc(1);
            }
            current_header = line.trim_start_matches('>').to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_header.is_empty() {
        records.push((current_header.clone(), current_seq.clone()));
        pb.inc(1);
    }
    pb.finish_with_message("Finished processing FASTA records.");
    Ok(records)
}

/// Main function for the predict subcommand.
pub fn run(args: PredictArgs) -> Result<(), Box<dyn Error>> {
    println!("INFO: System start time: {}", Local::now().format("%Y-%m-%d %H:%M:%S"));
    let vectorizer_path = format!("{}_vectorizer.bin.gz", args.model_base);
    let label_encoder_path = format!("{}_label_encoder.bin.gz", args.model_base);
    let model_path = format!("{}_rf_model.bin.gz", args.model_base);

    // Load artifacts.
    let vec_file = File::open(&vectorizer_path)?;
    let mut vec_decoder = GzDecoder::new(vec_file);
    let vectorizer: CountVectorizer = bincode::deserialize_from(&mut vec_decoder)?;
    println!("INFO: Loaded vectorizer from {}", vectorizer_path);

    let label_file = File::open(&label_encoder_path)?;
    let mut label_decoder = GzDecoder::new(label_file);
    let label_encoder: LabelEncoder = bincode::deserialize_from(&mut label_decoder)?;
    println!("INFO: Loaded label encoder from {}", label_encoder_path);

    let model_file = File::open(&model_path)?;
    let mut model_decoder = GzDecoder::new(model_file);
    let model: RandomForestClassifier<f64, usize, DenseMatrix<f64>, Vec<usize>> =
        bincode::deserialize_from(&mut model_decoder)?;
    println!("INFO: Loaded model from {}", model_path);

    println!("INFO: Reading input FASTA file: {}", args.fasta);
    let records = read_fasta_for_prediction(&args.fasta)?;
    if records.is_empty() {
        return Err("No records found in the input FASTA file.".into());
    }
    println!("INFO: Read {} records.", records.len());

    let texts: Vec<String> = records.iter().map(|(_, seq)| kmerize(seq, args.kmer_size)).collect();
    let x_data = vectorizer.transform(&texts);
    let x_matrix = DenseMatrix::from_2d_vec(&x_data)
        .map_err(|_| "Failed to create feature matrix")?;

    let predict_start = Instant::now();
    let y_pred = model.predict(&x_matrix)
        .map_err(|e| format!("Error during prediction: {:?}", e))?;
    println!("INFO: Prediction completed in {:.2} seconds.", predict_start.elapsed().as_secs_f32());

    let default_prediction = String::from("Unknown");
    let predictions: Vec<String> = y_pred.iter().map(|&class| {
        label_encoder.int_to_label.get(class).cloned().unwrap_or_else(|| default_prediction.clone())
    }).collect();

    let mut output_file = File::create(&args.output)?;
    writeln!(output_file, "Header\tPredicted_Lineage")?;
    for (i, (header, _)) in records.iter().enumerate() {
        let pred = predictions.get(i).unwrap_or(&default_prediction);
        writeln!(output_file, "{}\t{}", header, pred)?;
    }
    println!("INFO: Predictions written to {}", args.output);
    println!("INFO: System finish time: {}", Local::now().format("%Y-%m-%d %H:%M:%S"));
    Ok(())
}
