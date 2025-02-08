use clap::Parser;
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
use chrono::Local;

// Default k-mer size if not provided.
const DEFAULT_KMER_SIZE: usize = 4;

/// Converts a genomic sequence into overlapping k-mers separated by spaces.
/// For example, "ATGCAT" with k=3 becomes "ATG TGC GCA CAT".
fn kmerize(sequence: &str, k: usize) -> String {
    if sequence.len() < k {
        return String::new();
    }
    (0..=sequence.len() - k)
        .map(|i| &sequence[i..i + k])
        .collect::<Vec<&str>>()
        .join(" ")
}

/// Domain type representing a genome sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct GenomeSequence(String);

impl GenomeSequence {
    /// Creates a new GenomeSequence after validating that it is not empty.
    pub fn new<S: Into<String>>(s: S) -> Result<Self, String> {
        let s = s.into();
        if s.trim().is_empty() {
            Err("Genome sequence cannot be empty".into())
        } else {
            Ok(Self(s))
        }
    }
    /// Returns the inner string slice.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::ops::Deref for GenomeSequence {
    type Target = str;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Domain type representing a lineage.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct Lineage(String);

impl Lineage {
    /// Creates a new Lineage after validating that it is not empty.
    pub fn new<S: Into<String>>(s: S) -> Result<Self, String> {
        let s = s.into();
        if s.trim().is_empty() {
            Err("Lineage cannot be empty".into())
        } else {
            Ok(Self(s))
        }
    }
    /// Returns the inner string slice.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::ops::Deref for Lineage {
    type Target = str;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Command-line arguments for prediction.
/// The input is a FASTA file, and the model base name is used to infer the artifact filenames.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct PredictArgs {
    /// Input FASTA file (multifasta; header is expected to be in the format "Lineage_sequenceID")
    #[arg(long)]
    fasta: String,
    /// Base name for the saved model artifacts.
    /// For example, if you provide "A2_mode", then the program will load:
    ///   - "A2_mode_rf_model.bin.gz"
    ///   - "A2_mode_vectorizer.bin.gz"
    ///   - "A2_mode_label_encoder.bin.gz"
    #[arg(long)]
    model_base: String,
    /// Output file where predictions will be written.
    #[arg(long)]
    output: String,
    /// k-mer size (default is 21)
    #[arg(long, default_value_t = DEFAULT_KMER_SIZE)]
    kmer_size: usize,
}

/// A simple count vectorizer that splits texts on whitespace.
#[derive(Serialize, Deserialize, Debug)]
struct CountVectorizer {
    vocabulary: HashMap<String, usize>,
    feature_names: Vec<String>,
}

impl CountVectorizer {
    pub fn new() -> Self {
        Self {
            vocabulary: HashMap::new(),
            feature_names: Vec::new(),
        }
    }
    /// Fits the vectorizer on a collection of texts.
    /// T can be any type that can be converted to a &str.
    pub fn fit<T: AsRef<str>>(&mut self, texts: &[T]) {
        let mut freq: HashMap<String, usize> = HashMap::new();
        for text in texts {
            for token in text.as_ref().split_whitespace() {
                *freq.entry(token.to_string()).or_insert(0) += 1;
            }
        }
        // Sort tokens by frequency in descending order.
        let mut freq_vec: Vec<(String, usize)> = freq.into_iter().collect();
        freq_vec.sort_by(|a, b| b.1.cmp(&a.1));
        self.vocabulary = freq_vec
            .iter()
            .enumerate()
            .map(|(i, (token, _))| (token.clone(), i))
            .collect();
        self.feature_names = freq_vec.into_iter().map(|(token, _)| token).collect();
    }
    /// Transforms a collection of texts into a 2D vector (one row per text).
    /// Sequential processing is used for memory efficiency.
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
struct LabelEncoder {
    label_to_int: HashMap<String, usize>,
    int_to_label: Vec<String>,
}

impl LabelEncoder {
    pub fn new() -> Self {
        Self {
            label_to_int: HashMap::new(),
            int_to_label: Vec::new(),
        }
    }
    /// Fits the encoder on a collection of labels.
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
    /// Transforms a collection of labels into numeric values.
    pub fn transform<T: AsRef<str>>(&self, labels: &[T]) -> Vec<usize> {
        labels
            .iter()
            .map(|label| *self.label_to_int.get(label.as_ref()).unwrap())
            .collect()
    }
}

/// Reads a FASTA file and returns a vector of (header, sequence) tuples.
/// The header is expected to be in the format "Lineage_sequenceID".
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
    pb.finish_with_message("Finished processing FASTA records for prediction.");
    Ok(records)
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line arguments.
    let args = PredictArgs::parse();
    let kmer_size = args.kmer_size;

    // Log system start time.
    println!("INFO: System start time: {}", Local::now().format("%Y-%m-%d %H:%M:%S"));

    // Build file paths from the provided model base name.
    let model_path = format!("{}_rf_model.bin.gz", args.model_base);
    let vectorizer_path = format!("{}_vectorizer.bin.gz", args.model_base);
    let label_encoder_path = format!("{}_label_encoder.bin.gz", args.model_base);

    // Load saved vectorizer.
    let vec_file = File::open(&vectorizer_path)?;
    let mut vec_decoder = GzDecoder::new(vec_file);
    let vectorizer: CountVectorizer = bincode::deserialize_from(&mut vec_decoder)?;
    println!("INFO: Loaded vectorizer from {}", vectorizer_path);

    // Load saved label encoder.
    let label_file = File::open(&label_encoder_path)?;
    let mut label_decoder = GzDecoder::new(label_file);
    let label_encoder: LabelEncoder = bincode::deserialize_from(&mut label_decoder)?;
    println!("INFO: Loaded label encoder from {}", label_encoder_path);

    // Load saved model.
    let model_file = File::open(&model_path)?;
    let mut model_decoder = GzDecoder::new(model_file);
    let model: RandomForestClassifier<f64, usize, DenseMatrix<f64>, Vec<usize>> =
        bincode::deserialize_from(&mut model_decoder)?;
    println!("INFO: Loaded model from {}", model_path);

    // Read FASTA input.
    println!("INFO: Reading input FASTA file: {}", args.fasta);
    let records = read_fasta_for_prediction(&args.fasta)?;
    if records.is_empty() {
        return Err("No records found in the input FASTA file.".into());
    }
    println!("INFO: Read {} records.", records.len());

    // Convert each sequence to overlapping k-mers.
    let texts: Vec<String> = records
        .iter()
        .map(|(_, seq)| kmerize(seq, kmer_size))
        .collect();

    // Transform sequences using the loaded vectorizer.
    let x_data = vectorizer.transform(&texts);
    let x_matrix = DenseMatrix::from_2d_vec(&x_data)
        .expect("Failed to create feature matrix");

    // Predict using the loaded model.
    let predict_start = Instant::now();
    let y_pred = model
        .predict(&x_matrix)
        .map_err(|e| format!("Error during prediction: {:?}", e))?;
    println!(
        "INFO: Prediction completed in {:.2} seconds.",
        predict_start.elapsed().as_secs_f32()
    );

    // Map numeric predictions back to labels.
    let default_prediction = String::from("Unknown");
    let predictions: Vec<String> = y_pred
        .iter()
        .map(|&class| {
            label_encoder
                .int_to_label
                .get(class)
                .cloned()
                .unwrap_or_else(|| default_prediction.clone())
        })
        .collect();

    // Write predictions to the output file.
    let mut output_file = File::create(&args.output)?;
    writeln!(output_file, "Header\tPredicted_Lineage")?;
    // Bind a default string to avoid temporary drop issues.
    let default_str = String::from("Unknown");
    for (i, (header, _)) in records.iter().enumerate() {
        // Use a local binding for the default so that the reference lives long enough.
        let pred = predictions.get(i).unwrap_or(&default_str);
        writeln!(output_file, "{}\t{}", header, pred)?;
    }
    println!("INFO: Predictions written to {}", args.output);

    // Log system finish time.
    println!("INFO: System finish time: {}", Local::now().format("%Y-%m-%d %H:%M:%S"));

    Ok(())
}
