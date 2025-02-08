use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::ProgressBar;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};
use smartcore::ensemble::random_forest_classifier::{
    RandomForestClassifier, RandomForestClassifierParameters,
};
use smartcore::linalg::basic::matrix::DenseMatrix;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use chrono;
use smartcore::tree::decision_tree_classifier::SplitCriterion;
// Default k-mer size if not provided.
const DEFAULT_KMER_SIZE: usize = 6;

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

/// Command-line arguments for training. Only FASTA input is accepted.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTA file (multifasta; header is expected to be in the format "Lineage_sequenceID")
    #[arg(long)]
    fasta: String,
    /// Base name for the output model artifacts.
    #[arg(short, long)]
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
        // Sort tokens by frequency (descending).
        let mut freq_vec: Vec<(String, usize)> = freq.into_iter().collect();
        freq_vec.sort_by(|a, b| b.1.cmp(&a.1));
        // (Optional) You can limit the vocabulary size here if needed.
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

/// Reads a FASTA file and returns a tuple: (vector of GenomeSequence, vector of Lineage).
/// The header is expected to be in the format "Lineage_sequenceID"; the lineage is taken as the part before the underscore.
/// A progress bar is used to indicate processing.
fn read_fasta(path: &str, _k: usize) -> Result<(Vec<GenomeSequence>, Vec<Lineage>), Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let pb = ProgressBar::new_spinner();
    pb.set_message("Processing FASTA records...");
    let mut sequences = Vec::new();
    let mut lineages = Vec::new();
    let mut current_seq = String::new();
    let mut current_lineage = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_lineage.is_empty() {
                sequences.push(GenomeSequence::new(current_seq.clone())?);
                lineages.push(Lineage::new(current_lineage.clone())?);
                pb.inc(1);
            }
            // Extract lineage from header (portion before the underscore).
            let header = line.trim_start_matches('>');
            let lineage_part = header.split('_').next().unwrap_or(header);
            current_lineage = lineage_part.to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_lineage.is_empty() {
        sequences.push(GenomeSequence::new(current_seq)?);
        lineages.push(Lineage::new(current_lineage)?);
        pb.inc(1);
    }
    pb.finish_with_message("Finished processing FASTA file.");
    Ok((sequences, lineages))
}

/// Loads input data from a FASTA file, returning vectors of domain types.
fn load_input_data(args: &Args, k: usize) -> Result<(Vec<GenomeSequence>, Vec<Lineage>), Box<dyn Error>> {
    println!("INFO: Reading input FASTA file: {}", args.fasta);
    read_fasta(&args.fasta, k)
}

/// Prepares the feature matrix and label vector by converting sequences into k-mer strings,
/// then fitting the vectorizer and label encoder.
fn prepare_data(texts: &[String], labels: &[String], kmer_size: usize) -> Result<(CountVectorizer, LabelEncoder, Vec<Vec<f64>>, Vec<usize>), Box<dyn Error>> {
    // Convert each genome sequence into overlapping k-mers.
    let kmer_texts: Vec<String> = texts.iter().map(|s| kmerize(s, kmer_size)).collect();
    let mut vectorizer = CountVectorizer::new();
    vectorizer.fit(&kmer_texts);
    let x_data = vectorizer.transform(&kmer_texts);
    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(labels);
    let y = label_encoder.transform(labels);
    if label_encoder.int_to_label.len() < 2 {
        return Err("Training data must contain at least two distinct classes.".into());
    }
    Ok((vectorizer, label_encoder, x_data, y))
}

/// Splits the data into training and testing sets.
fn split_train_test(x_data: &[Vec<f64>], y: &[usize], test_ratio: f64) -> Result<(DenseMatrix<f64>, Vec<usize>, DenseMatrix<f64>, Vec<usize>), Box<dyn Error>> {
    let n_samples = x_data.len();
    let raw_test_size = ((n_samples as f64) * test_ratio).round() as usize;
    let test_size = if raw_test_size == 0 && n_samples > 1 { 1 } else { raw_test_size };
    if n_samples - test_size == 0 {
        return Err("Not enough samples for training after splitting.".into());
    }
    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    indices.shuffle(&mut rng);
    let test_indices = &indices[..test_size];
    let train_indices = &indices[test_size..];
    let train_data: Vec<Vec<f64>> = train_indices.iter().map(|&i| x_data[i].clone()).collect();
    let test_data: Vec<Vec<f64>> = test_indices.iter().map(|&i| x_data[i].clone()).collect();
    let y_test: Vec<usize> = test_indices.iter().map(|&i| y[i]).collect();
    let y_train: Vec<usize> = train_indices.iter().map(|&i| y[i]).collect();
    let x_train = DenseMatrix::from_2d_vec(&train_data).expect("Failed to create training matrix");
    let x_test = DenseMatrix::from_2d_vec(&test_data).expect("Failed to create test matrix");
    Ok((x_train, y_train, x_test, y_test))
}

/// Saves the model artifacts (vectorizer, label encoder, and the trained model).
fn save_artifacts<M>(model: &M, vectorizer: &CountVectorizer, label_encoder: &LabelEncoder, output_base: &str) -> Result<(), Box<dyn Error>>
where
    M: Serialize,
{
    let vectorizer_filename = format!("{}_vectorizer.bin.gz", output_base);
    let vec_file = File::create(&vectorizer_filename)?;
    let mut vec_encoder = GzEncoder::new(vec_file, Compression::default());
    bincode::serialize_into(&mut vec_encoder, vectorizer)?;
    vec_encoder.finish()?;
    println!("INFO: Saved the vectorizer to {}", vectorizer_filename);

    let label_filename = format!("{}_label_encoder.bin.gz", output_base);
    let label_file = File::create(&label_filename)?;
    let mut label_encoder_gz = GzEncoder::new(label_file, Compression::default());
    bincode::serialize_into(&mut label_encoder_gz, label_encoder)?;
    label_encoder_gz.finish()?;
    println!("INFO: Saved the label encoder to {}", label_filename);

    let model_filename = format!("{}_rf_model.bin.gz", output_base);
    let model_file = File::create(&model_filename)?;
    let mut model_encoder = GzEncoder::new(model_file, Compression::default());
    bincode::serialize_into(&mut model_encoder, model)?;
    model_encoder.finish()?;
    println!("INFO: Saved the model to {}", model_filename);

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line arguments.
    let args = Args::parse();
    let kmer_size = args.kmer_size;

    // Log the current system start time.
    println!("INFO: System start time: {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"));

    let overall_start = Instant::now();

    // Load FASTA input.
    let input_start = Instant::now();
    let (genome_vec, lineage_vec) = load_input_data(&args, kmer_size)?;
    println!(
        "INFO: Finished reading input FASTA file in {:.2} seconds.",
        input_start.elapsed().as_secs_f32()
    );

    let texts: Vec<String> = genome_vec.iter().map(|g| g.as_str().to_string()).collect();
    let labels: Vec<String> = lineage_vec.iter().map(|l| l.as_str().to_string()).collect();

    // Prepare data: convert sequences into k-mer strings, then fit vectorizer and label encoder.
    let prep_start = Instant::now();
    let (vectorizer, label_encoder, x_data, y) = prepare_data(&texts, &labels, kmer_size)?;
    println!(
        "INFO: Data preparation completed in {:.2} seconds.",
        prep_start.elapsed().as_secs_f32()
    );

    // Split data.
    let (x_train, y_train, x_test, y_test) = split_train_test(&x_data, &y, 0.2)?;

    println!("INFO: Starting to train the model: {}", args.output);
    let train_start = Instant::now();

    let rf_params = RandomForestClassifierParameters {
        max_depth: None,
        min_samples_leaf: 1,
        min_samples_split: 2,
        n_trees: 100,
        m: Some((vectorizer.vocabulary.len() as f64).sqrt().floor() as usize),
        seed: 42,
        criterion: SplitCriterion::Gini,
        keep_samples: true,
    };

    let clf = RandomForestClassifier::fit(&x_train, &y_train, rf_params)
        .map_err(|e| format!("Error training model: {:?}", e))?;

    println!(
        "INFO: Model training completed in {:.2} seconds.",
        train_start.elapsed().as_secs_f32()
    );

    let pred_start = Instant::now();
    let y_pred = clf
        .predict(&x_test)
        .map_err(|e| format!("Error during prediction: {:?}", e))?;
    println!(
        "INFO: Prediction on test set completed in {:.2} seconds.",
        pred_start.elapsed().as_secs_f32()
    );

    let correct = y_pred
        .iter()
        .zip(y_test.iter())
        .filter(|(&pred, &true_val)| pred == true_val)
        .count();
    let accuracy = correct as f64 / y_test.len() as f64;
    println!("INFO: Model's accuracy on the test set is: {:.2}%", accuracy * 100.0);
    println!("INFO: Finished training the model: {}", args.output);

    let save_start = Instant::now();
    save_artifacts(&clf, &vectorizer, &label_encoder, &args.output)?;
    println!(
        "INFO: Saving artifacts completed in {:.2} seconds.",
        save_start.elapsed().as_secs_f32()
    );

    println!(
        "INFO: Overall process completed in {:.2} seconds.",
        overall_start.elapsed().as_secs_f32()
    );

    // Log the current system finish time.
    println!("INFO: System finish time: {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"));
    
    Ok(())
}
