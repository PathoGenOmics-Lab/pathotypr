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
use smartcore::linalg::basic::matrix::DenseMatrix;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use smartcore::tree::decision_tree_classifier::SplitCriterion;

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

    /// Transform an array of texts into a 2D vector (one row per text).
    pub fn transform(&self, texts: &[String]) -> Vec<Vec<f64>> {
        texts
            .par_iter()
            .map(|text| {
                let n_features = self.vocabulary.len();
                let mut counts = vec![0.0; n_features];
                for token in text.split_whitespace() {
                    if let Some(&idx) = self.vocabulary.get(token) {
                        counts[idx] += 1.0;
                    }
                }
                counts
            })
            .collect()
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
            // Save previous record if it exists.
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

    let k = 21;
    let mut texts: Vec<String> = Vec::new();
    let mut labels: Vec<String> = Vec::new();

    // Process input based on whether --tsv or --fasta was provided.
    if let Some(tsv_file) = args.tsv {
        println!("INFO: Reading input TSV file: {}", tsv_file);
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&tsv_file)?;
        let headers = rdr.headers()?.clone();
        let seq_idx = headers
            .iter()
            .position(|h| h == "genome_sequence")
            .ok_or("Missing 'genome_sequence' column in header")?;
        let lineage_idx = headers
            .iter()
            .position(|h| h == "lineage")
            .ok_or("Missing 'lineage' column in header")?;
        let pb = ProgressBar::new_spinner();
        pb.set_message("Processing TSV records...");

        for result in rdr.records() {
            let record: StringRecord = result?;
            let genome_sequence = record
                .get(seq_idx)
                .ok_or("Missing genome_sequence field")?;
            let lineage = record
                .get(lineage_idx)
                .ok_or("Missing lineage field")?;
            let kmers = split_kmers(genome_sequence, k);
            let joined = kmers.join(" ");
            texts.push(joined);
            labels.push(lineage.to_string());
            pb.inc(1);
        }
        pb.finish_with_message("Finished processing TSV file.");
        println!("INFO: Finished processing the input TSV file: {}", tsv_file);
    } else if let Some(fasta_file) = args.fasta {
        println!("INFO: Reading input FASTA file: {}", fasta_file);
        let (seqs, fasta_labels) = read_fasta(&fasta_file)?;
        texts = seqs
            .into_par_iter()
            .map(|seq| {
                let kmers = split_kmers(&seq, k);
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
    // Instead of returning a DenseMatrix directly, we return a 2D vector.
    let x_data = vectorizer.transform(&texts);
    // Build a DenseMatrix from the 2D vector (if needed for training).
    let _x = DenseMatrix::from_2d_vec(&x_data).expect("Failed to create matrix");

    // Encode the labels.
    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(&labels);
    let y = label_encoder.transform(&labels);

    // Use the number of samples and features.
    let n_samples = x_data.len();
    let n_features = vectorizer.vocabulary.len();

    // Perform an 80/20 train/test split.
    let test_size = ((n_samples as f64) * 0.2).round() as usize;
    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    indices.shuffle(&mut rng);
    let test_indices = &indices[..test_size];
    let train_indices = &indices[test_size..];

    // Partition the 2D data and labels.
    let train_data: Vec<Vec<f64>> = train_indices.iter().map(|&i| x_data[i].clone()).collect();
    let test_data: Vec<Vec<f64>> = test_indices.iter().map(|&i| x_data[i].clone()).collect();
    let y_test: Vec<usize> = test_indices.iter().map(|&i| y[i]).collect();

    let x_train = DenseMatrix::from_2d_vec(&train_data).expect("Failed to create training matrix");
    let x_test = DenseMatrix::from_2d_vec(&test_data).expect("Failed to create test matrix");

    println!("INFO: Starting to train the model: {}", args.output);

    // Build the RandomForest parameters.
    let rf_params = RandomForestClassifierParameters {
        max_depth: None,
        min_samples_leaf: 1,
        min_samples_split: 2,
        n_trees: 100,
        m: Some((n_features as f64).sqrt().floor() as usize),
        seed: 42,
        criterion: SplitCriterion::Gini, // Using the default criterion.
        keep_samples: true,
    };

    let clf = RandomForestClassifier::fit(&x_train, &y, rf_params)
        .map_err(|e| format!("Error training model: {:?}", e))?;

    let y_pred = clf
        .predict(&x_test)
        .map_err(|e| format!("Error during prediction: {:?}", e))?;

    let correct = y_pred
        .iter()
        .zip(y_test.iter())
        .filter(|(&pred, &true_val)| pred == true_val)
        .count();
    let accuracy = correct as f64 / y_test.len() as f64;
    println!("INFO: Model's accuracy on the test set is: {}", accuracy);
    println!("INFO: Finished training the model: {}", args.output);

    // Serialize the vectorizer.
    {
        let vectorizer_filename = format!("{}_vectorizer.bin.gz", args.output);
        let vec_file = File::create(&vectorizer_filename)?;
        let mut vec_encoder = GzEncoder::new(vec_file, Compression::default());
        bincode::serialize_into(&mut vec_encoder, &vectorizer)?;
        vec_encoder.finish()?;
        println!("INFO: Saved the vectorizer to {}", vectorizer_filename);
    }
    // Serialize the label encoder.
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
