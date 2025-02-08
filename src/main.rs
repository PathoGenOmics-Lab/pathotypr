use clap::Parser;
use csv::StringRecord;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::ProgressBar;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use smartcore::ensemble::random_forest_classifier::{
    RandomForestClassifier, RandomForestClassifierParameters,
};
use nalgebra::DenseMatrix; // Changed from smartcore::linalg::dense_matrix
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader}; // Removed unused Write import

/// Command line arguments.
/// Either --tsv or --fasta must be provided for the input.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input TSV file (tab-separated with headers "genome_sequence" and "lineage")
    #[arg(long, conflicts_with = "fasta", required_unless_present = "fasta")]
    tsv: Option<String>,

    /// Input FASTA file (multifasta file; header is used as label)
    #[arg(long, conflicts_with = "tsv", required_unless_present = "tsv")]
    fasta: Option<String>,

    /// Base name for the output model files.
    #[arg(short, long)]
    output: String,
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

    /// Fit the vectorizer on an array of texts.
    pub fn fit(&mut self, texts: &[String]) {
        for text in texts {
            for token in text.split_whitespace() {
                if !self.vocabulary.contains_key(token) {
                    let index = self.vocabulary.len();
                    self.vocabulary.insert(token.to_string(), index);
                }
            }
        }
        // Ensure that the feature names are in order of insertion.
        let mut vocab_vec: Vec<(String, usize)> = self
            .vocabulary
            .iter()
            .map(|(token, &index)| (token.clone(), index))
            .collect();
        vocab_vec.sort_by_key(|&(_, index)| index);
        self.feature_names = vocab_vec.into_iter().map(|(token, _)| token).collect();
    }

    /// Transform an array of texts into a DenseMatrix of counts.
    pub fn transform(&self, texts: &[String]) -> DenseMatrix<f64> {
        let n_samples = texts.len();
        let n_features = self.vocabulary.len();
        // Process each text in parallel.
        let row_data: Vec<Vec<f64>> = texts
            .par_iter()
            .map(|text| {
                let mut counts = vec![0.0; n_features];
                for token in text.split_whitespace() {
                    if let Some(&idx) = self.vocabulary.get(token) {
                        counts[idx] += 1.0;
                    }
                }
                counts
            })
            .collect();

        // Flatten the 2D vector into a 1D vector.
        let mut data = Vec::with_capacity(n_samples * n_features);
        for row in row_data {
            data.extend_from_slice(&row);
        }
        DenseMatrix::from_array(n_samples, n_features, &data)
    }

    /// Get the feature names.
    pub fn get_feature_names(&self) -> &Vec<String> {
        &self.feature_names
    }
}

/// Label encoder for mapping categorical labels (strings) to numeric values.
#[derive(Serialize, Deserialize, Debug)]
struct LabelEncoder {
    label_to_int: HashMap<String, usize>,
    int_to_label: Vec<String>, // index -> label
}

impl LabelEncoder {
    pub fn new() -> Self {
        Self {
            label_to_int: HashMap::new(),
            int_to_label: Vec::new(),
        }
    }

    pub fn fit(&mut self, labels: &[String]) {
        for label in labels {
            if !self.label_to_int.contains_key(label) {
                let index = self.int_to_label.len();
                self.label_to_int.insert(label.clone(), index);
                self.int_to_label.push(label.clone());
            }
        }
    }

    pub fn transform(&self, labels: &[String]) -> Vec<usize> {
        labels
            .iter()
            .map(|label| *self.label_to_int.get(label).unwrap())
            .collect()
    }
}

/// Splits a sequence into k-mers.
fn split_kmers(sequence: &str, k: usize) -> Vec<String> {
    let seq_len = sequence.len();
    if seq_len < k {
        return Vec::new();
    }
    (0..=seq_len - k)
        .map(|i| sequence[i..i + k].to_string())
        .collect()
}

/// Reads a multifasta file and returns a tuple: (vector of sequences, vector of labels).
fn read_fasta(path: &str) -> Result<(Vec<String>, Vec<String>), Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut labels = Vec::new();
    let mut current_seq = String::new();
    let mut current_label = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Save previous record if exists.
            if !current_label.is_empty() {
                sequences.push(current_seq.clone());
                labels.push(current_label.clone());
            }
            current_label = line.trim_start_matches('>').to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    // Save the last record.
    if !current_label.is_empty() {
        sequences.push(current_seq);
        labels.push(current_label);
    }
    Ok((sequences, labels))
}

fn main() -> Result<(), Box<dyn Error>> {
    // Parse command line arguments.
    let args = Args::parse();

    // Parameters
    let chunksize = 10;
    let k = 21;
    let mut texts: Vec<String> = Vec::new();
    let mut labels: Vec<String> = Vec::new();

    // Process input based on whether --tsv or --fasta was provided.
    if let Some(tsv_file) = args.tsv {
        println!("INFO: Reading input TSV file: {}", tsv_file);
        // Open the CSV file with tab delimiter.
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&tsv_file)?;

        // The TSV file is expected to have headers.
        let headers = rdr.headers()?.clone();
        // Determine the indices of the needed columns.
        let seq_idx = headers
            .iter()
            .position(|h| h == "genome_sequence")
            .ok_or("Missing 'genome_sequence' column in header")?;
        let lineage_idx = headers
            .iter()
            .position(|h| h == "lineage")
            .ok_or("Missing 'lineage' column in header")?;

        // Initialize a progress bar.
        let pb = ProgressBar::new_spinner();
        pb.set_message("Processing TSV records...");

        let mut record_iter = rdr.records();
        // Process TSV file in chunks.
        loop {
            let mut chunk_texts = Vec::with_capacity(chunksize);
            let mut chunk_labels = Vec::with_capacity(chunksize);
            for _ in 0..chunksize {
                if let Some(result) = record_iter.next() {
                    let record: StringRecord = result?;
                    let genome_sequence = record
                        .get(seq_idx)
                        .ok_or("Missing genome_sequence field")?;
                    let lineage = record
                        .get(lineage_idx)
                        .ok_or("Missing lineage field")?;
                    // Compute kmers and join them with whitespace.
                    let kmers = split_kmers(genome_sequence, k);
                    let joined = kmers.join(" ");
                    chunk_texts.push(joined);
                    chunk_labels.push(lineage.to_string());
                    pb.inc(1);
                } else {
                    break;
                }
            }
            if chunk_texts.is_empty() {
                break;
            }
            texts.extend(chunk_texts);
            labels.extend(chunk_labels);
        }
        pb.finish_with_message("Finished processing TSV file.");
        println!("INFO: Finished processing the input TSV file: {}", tsv_file);
    } else if let Some(fasta_file) = args.fasta {
        println!("INFO: Reading input FASTA file: {}", fasta_file);
        let (sequences, fasta_labels) = read_fasta(&fasta_file)?;
        // Process each sequence in parallel to compute kmers.
        texts = sequences
            .par_iter()
            .map(|seq| {
                let kmers = split_kmers(seq, k);
                kmers.join(" ")
            })
            .collect();
        labels = fasta_labels;
        println!("INFO: Finished processing the input FASTA file: {}", fasta_file);
    } else {
        return Err("Either --tsv or --fasta must be provided as input.".into());
    }

    // Fit the vectorizer on all texts.
    let mut vectorizer = CountVectorizer::new();
    vectorizer.fit(&texts);
    let X = vectorizer.transform(&texts);

    // Encode the labels into integers.
    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(&labels);
    let y: Vec<usize> = label_encoder.transform(&labels);

    // Perform an 80/20 train/test split.
    let n_samples = X.shape().0;
    let test_size = ((n_samples as f64) * 0.2).round() as usize;
    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    indices.shuffle(&mut rng);
    let test_indices = &indices[..test_size];
    let train_indices = &indices[test_size..];

    // Create train and test matrices.
    let n_features = X.shape().1;
    let mut x_train_data = Vec::with_capacity(train_indices.len() * n_features);
    let mut y_train = Vec::with_capacity(train_indices.len());
    let mut x_test_data = Vec::with_capacity(test_indices.len() * n_features);
    let mut y_test = Vec::with_capacity(test_indices.len());

    for &i in train_indices {
        let row = X.get_row(i);
        x_train_data.extend_from_slice(row);
        y_train.push(y[i]);
    }
    for &i in test_indices {
        let row = X.get_row(i);
        x_test_data.extend_from_slice(row);
        y_test.push(y[i]);
    }
    let x_train = DenseMatrix::from_array(train_indices.len(), n_features, &x_train_data);
    let x_test = DenseMatrix::from_array(test_indices.len(), n_features, &x_test_data);

    println!("INFO: Starting to train the model: {}", args.output);

    // Train the random forest classifier.
    let rf_params = RandomForestClassifierParameters {
        n_trees: 100,
        criterion: Default::default(), // Default criterion
        keep_samples: true, // Whether to keep samples for OOB error estimation
        seed: 42, // Changed from Some(42) to 42
    };

    let clf = RandomForestClassifier::fit(&x_train, &y_train, rf_params)
        .map_err(|e| format!("Error training model: {:?}", e))?;

    // Predict on the test set.
    let y_pred = clf
        .predict(&x_test)
        .map_err(|e| format!("Error during prediction: {:?}", e))?;

    // Compute accuracy.
    let correct = y_pred
        .iter()
        .zip(y_test.iter())
        .filter(|(&pred, &true_val)| pred == true_val)
        .count();
    let accuracy = correct as f64 / y_test.len() as f64;
    println!("INFO: Model's accuracy on the test set is: {}", accuracy);
    println!("INFO: Finished training the model: {}", args.output);

    // Save the model, vectorizer, and label encoder to gzipped files.
    {
        let model_filename = format!("{}_rfm.bin.gz", args.output);
        let model_file = File::create(&model_filename)?;
        let mut encoder = GzEncoder::new(model_file, Compression::default());
        bincode::serialize_into(&mut encoder, &clf)?;
        encoder.finish()?;
        println!("INFO: Saved the random forest model to {}", model_filename);
    }
    {
        let vectorizer_filename = format!("{}_vectorizer.bin.gz", args.output);
        let vec_file = File::create(&vectorizer_filename)?;
        let mut vec_encoder = GzEncoder::new(vec_file, Compression::default());
        bincode::serialize_into(&mut vec_encoder, &vectorizer)?;
        vec_encoder.finish()?;
        println!("INFO: Saved the vectorizer to {}", vectorizer_filename);
    }
    {
        let label_filename = format!("{}_label_encoder.bin.gz", args.output);
        let label_file = File::create(&label_filename)?;
        let mut label_encoder_gz = GzEncoder::new(label_file, Compression::default());
        bincode::serialize_into(&mut label_encoder_gz, &label_encoder)?;
        label_encoder_gz.finish()?;
        println!("INFO: Saved the label encoder to {}", label_filename);
    }

    Ok(())
}