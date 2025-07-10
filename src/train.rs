//! This module handles the `train` subcommand.
//!
//! It includes functionality for:
//! 1. Reading and parsing FASTA files, and warning about imbalanced classes.
//! 2. Converting genome sequences into k-mer features in parallel.
//! 3. Training a Random Forest classifier using the `smartcore` library.
//! 4. Bundling the model, vectorizer, label encoder, and configuration into a single file.

use crate::common::{kmerize, CountVectorizer, LabelEncoder, ModelBundle, ModelConfig};
use crate::errors::{AppError, AppResult};
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, trace, warn};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use smartcore::ensemble::random_forest_classifier::{
    RandomForestClassifier, RandomForestClassifierParameters,
};
use smartcore::linalg::basic::matrix::DenseMatrix;
use smartcore::tree::decision_tree_classifier::SplitCriterion;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Deref;
use std::time::Instant;

const DEFAULT_KMER_SIZE: usize = 6;
const DEFAULT_TEST_SPLIT: f64 = 0.2;

#[derive(Parser, Debug)]
pub struct Args {
    /// Path to the input multifasta file. Headers must be in Lineage_sequenceID format.
    #[arg(short = 'i', long)]
    pub input: String,
    /// Path for the unified output model file (e.g., my_model.pathotypr.gz).
    #[arg(short = 'o', long)]
    pub output: String,
    /// The size of the k-mers to generate from sequences.
    #[arg(short = 'k', long, default_value_t = DEFAULT_KMER_SIZE)]
    pub kmer_size: usize,
    /// Proportion of the data to use for the test set.
    #[arg(short = 's', long, default_value_t = DEFAULT_TEST_SPLIT)]
    pub test_split: f64,
    /// Number of CPU threads to use.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct GenomeSequence(String);
impl GenomeSequence {
    pub fn new<S: Into<String>>(s: S) -> Result<Self, String> {
        let s = s.into();
        if s.trim().is_empty() {
            Err("Genome sequence cannot be empty".into())
        } else {
            Ok(Self(s))
        }
    }
    pub fn inner(&self) -> &str {
        &self.0
    }
}
impl Deref for GenomeSequence {
    type Target = str;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Lineage(String);
impl Lineage {
    pub fn new<S: Into<String>>(s: S) -> Result<Self, String> {
        let s = s.into();
        if s.trim().is_empty() {
            Err("Lineage cannot be empty".into())
        } else {
            Ok(Self(s))
        }
    }
    pub fn inner(&self) -> &str {
        &self.0
    }
}
impl Deref for Lineage {
    type Target = str;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

fn check_class_balance(labels: &[Lineage]) {
    if labels.is_empty() {
        return;
    }
    let mut counts = HashMap::new();
    for label in labels {
        *counts.entry(label.inner()).or_insert(0) += 1;
    }
    trace!("Class counts: {:?}", counts);
    let min_count = counts.values().min().unwrap_or(&0);
    let max_count = counts.values().max().unwrap_or(&0);
    if *min_count < 5 {
        warn!(
            "Class imbalance: A class has only {} samples. Model may be unreliable.",
            min_count
        );
    }
    if *min_count > 0 && *max_count as f64 / *min_count as f64 > 10.0 {
        warn!(
            "Class imbalance: Ratio of largest ({}) to smallest ({}) class is > 10:1.",
            max_count, min_count
        );
    }
}

fn read_fasta(path: &str) -> AppResult<(Vec<GenomeSequence>, Vec<Lineage>)> {
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
    let mut sequences = Vec::new();
    let mut lineages = Vec::new();
    let mut current_seq = String::new();
    let mut current_lineage = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_lineage.is_empty() {
                sequences.push(GenomeSequence::new(current_seq.clone()).map_err(AppError::Parsing)?);
                lineages.push(Lineage::new(current_lineage.clone()).map_err(AppError::Parsing)?);
            }
            let header = line.trim_start_matches('>');
            let lineage_part = header.split('_').next().unwrap_or(header);
            current_lineage = lineage_part.to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_lineage.is_empty() {
        sequences.push(GenomeSequence::new(current_seq).map_err(AppError::Parsing)?);
        lineages.push(Lineage::new(current_lineage).map_err(AppError::Parsing)?);
    }
    pb.finish_with_message("FASTA reading complete");
    debug!("Read {} sequences and lineages from {}", sequences.len(), path);
    Ok((sequences, lineages))
}

fn load_input_data(args: &Args) -> AppResult<(Vec<GenomeSequence>, Vec<Lineage>)> {
    info!("▶ Reading input FASTA file: {}", args.input);
    let (sequences, lineages) = read_fasta(&args.input)?;
    check_class_balance(&lineages);
    Ok((sequences, lineages))
}

fn prepare_data(
    texts: &[String],
    labels: &[String],
    kmer_size: usize,
) -> AppResult<(CountVectorizer, LabelEncoder, Vec<Vec<f64>>, Vec<usize>)> {
    info!("▶ Generating k-mers and vectorizing features...");
    let kmer_texts: Vec<String> = texts.iter().map(|s| kmerize(s, kmer_size)).collect();
    let mut vectorizer = CountVectorizer::new();
    vectorizer.fit(&kmer_texts);
    debug!("Vectorizer fitted with a vocabulary of {} k-mers.", vectorizer.vocabulary.len());
    let x_data = vectorizer.transform(&kmer_texts);
    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(labels);
    debug!("Label encoder fitted with {} classes.", label_encoder.int_to_label.len());
    trace!("Class mapping: {:?}", label_encoder.label_to_int);
    let y = label_encoder.transform(labels);
    if label_encoder.int_to_label.len() < 2 {
        return Err(AppError::NotEnoughData(
            "Training data must contain at least two distinct classes.".to_string(),
        ));
    }
    Ok((vectorizer, label_encoder, x_data, y))
}

fn split_train_test(
    x_data: &[Vec<f64>],
    y: &[usize],
    test_ratio: f64,
) -> AppResult<(DenseMatrix<f64>, Vec<usize>, DenseMatrix<f64>, Vec<usize>)> {
    let n_samples = x_data.len();
    let raw_test_size = (n_samples as f64 * test_ratio).round() as usize;
    let test_size = if raw_test_size == 0 && n_samples > 1 {
        1
    } else {
        raw_test_size
    };
    let train_size = n_samples - test_size;
    debug!("Splitting data: {} training samples, {} test samples.", train_size, test_size);

    if train_size == 0 {
        return Err(AppError::NotEnoughData(
            "Not enough samples for training after splitting.".to_string(),
        ));
    }
    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    indices.shuffle(&mut rng);
    let test_indices = &indices[..test_size];
    let train_indices = &indices[test_size..];
    let train_data: Vec<Vec<f64>> = train_indices.iter().map(|&i| x_data[i].clone()).collect();
    let test_data: Vec<Vec<f64>> = test_indices.iter().map(|&i| x_data[i].clone()).collect();
    let y_train: Vec<usize> = train_indices.iter().map(|&i| y[i]).collect();
    let y_test: Vec<usize> = test_indices.iter().map(|&i| y[i]).collect();
    let x_train = DenseMatrix::from_2d_vec(&train_data)?;
    let x_test = DenseMatrix::from_2d_vec(&test_data)?;
    Ok((x_train, y_train, x_test, y_test))
}

fn save_model_bundle(bundle: &ModelBundle, output_path: &str) -> AppResult<()> {
    let file = File::create(output_path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    bincode::serialize_into(&mut encoder, bundle)?;
    encoder.finish()?;
    Ok(())
}

pub fn run(args: Args) -> AppResult<()> {
    if !(0.0..=1.0).contains(&args.test_split) {
        return Err(AppError::Generic(
            "--test-split must be between 0.0 and 1.0.".to_string(),
        ));
    }
    if let Some(n) = args.threads {
        debug!("Setting number of threads to: {}", n);
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .map_err(|e| AppError::Generic(format!("Failed to build thread pool: {}", e)))?;
    }

    let overall_start = Instant::now();
    let (genome_vec, lineage_vec) = load_input_data(&args)?;
    let texts: Vec<String> = genome_vec.iter().map(|g| g.inner().to_string()).collect();
    let labels: Vec<String> = lineage_vec.iter().map(|l| l.inner().to_string()).collect();
    let (vectorizer, label_encoder, x_data, y) = prepare_data(&texts, &labels, args.kmer_size)?;
    let (x_train, y_train, x_test, y_test) = split_train_test(&x_data, &y, args.test_split)?;

    info!("▶ Starting to train the model -> {}", args.output);
    let n_features = (vectorizer.vocabulary.len() as f64).sqrt().floor() as usize;
    let rf_params = RandomForestClassifierParameters {
        max_depth: None,
        min_samples_leaf: 1,
        min_samples_split: 2,
        n_trees: 100,
        m: Some(n_features),
        seed: 42,
        criterion: SplitCriterion::Gini,
        keep_samples: true,
    };
    debug!("Random Forest parameters: {:?}", rf_params);

    let train_start = Instant::now();
    let clf = RandomForestClassifier::fit(&x_train, &y_train, rf_params)?;
    info!(
        "  Model training completed in {:.2} seconds.",
        train_start.elapsed().as_secs_f32()
    );

    let y_pred = clf.predict(&x_test)?;
    let correct = y_pred
        .iter()
        .zip(y_test.iter())
        .filter(|(&pred, &true_val)| pred == true_val)
        .count();
    let accuracy = if y_test.is_empty() { 0.0 } else { correct as f64 / y_test.len() as f64 };
    info!(
        "  Model's accuracy on the test set is: {:.2}%",
        accuracy * 100.0
    );

    let config = ModelConfig {
        kmer_size: args.kmer_size,
    };
    let bundle = ModelBundle {
        config,
        vectorizer,
        label_encoder,
        model: clf,
    };
    save_model_bundle(&bundle, &args.output)?;
    info!("  Saved the complete model bundle to {}", args.output);

    info!(
        "✅ Process completed in {:.2} seconds.",
        overall_start.elapsed().as_secs_f32()
    );
    Ok(())
}
