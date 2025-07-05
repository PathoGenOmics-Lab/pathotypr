//! This module handles the `train` subcommand.
//!
//! It includes functionality for:
//! 1. Reading and parsing FASTA files, and warning about imbalanced classes.
//! 2. Converting genome sequences into k-mer features in parallel.
//! 3. Training a Random Forest classifier using the `smartcore` library.
//! 4. Bundling the model, vectorizer, label encoder, and configuration into a single file.

use anyhow::{anyhow, Result};
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
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

#[derive(Serialize, Deserialize, Debug)]
pub struct ModelConfig { pub kmer_size: usize }
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelBundle {
    pub config: ModelConfig,
    pub vectorizer: CountVectorizer,
    pub label_encoder: LabelEncoder,
    pub model: RandomForestClassifier<f64, usize, DenseMatrix<f64>, Vec<usize>>,
}

#[derive(Parser, Debug)]
pub struct Args {
    #[arg(short = 'i', long)] pub input: String,
    #[arg(short = 'o', long)] pub output: String,
    #[arg(short = 'k', long, default_value_t = DEFAULT_KMER_SIZE)] pub kmer_size: usize,
    #[arg(short = 's', long, default_value_t = DEFAULT_TEST_SPLIT)] pub test_split: f64,
    #[arg(short = 't', long)] pub threads: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct GenomeSequence(String);
impl GenomeSequence {
    pub fn new<S: Into<String>>(s: S) -> std::result::Result<Self, String> {
        let s = s.into(); if s.trim().is_empty() { Err("Genome sequence cannot be empty".into()) } else { Ok(Self(s)) }
    }
    pub fn inner(&self) -> &str { &self.0 }
}
impl Deref for GenomeSequence { type Target = str; fn deref(&self) -> &Self::Target { &self.0 } }

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
struct Lineage(String);
impl Lineage {
    pub fn new<S: Into<String>>(s: S) -> std::result::Result<Self, String> {
        let s = s.into(); if s.trim().is_empty() { Err("Lineage cannot be empty".into()) } else { Ok(Self(s)) }
    }
    pub fn inner(&self) -> &str { &self.0 }
}
impl Deref for Lineage { type Target = str; fn deref(&self) -> &Self::Target { &self.0 } }

#[derive(Serialize, Deserialize, Debug)]
pub struct CountVectorizer { vocabulary: HashMap<String, usize>, feature_names: Vec<String> }
impl CountVectorizer {
    pub fn new() -> Self { Self { vocabulary: HashMap::new(), feature_names: Vec::new() } }
    pub fn fit<T: AsRef<str>>(&mut self, texts: &[T]) {
        let mut freq: HashMap<String, usize> = HashMap::new();
        for text in texts { for token in text.as_ref().split_whitespace() { *freq.entry(token.to_string()).or_insert(0) += 1; } }
        let mut freq_vec: Vec<(String, usize)> = freq.into_iter().collect();
        freq_vec.sort_by(|a, b| b.1.cmp(&a.1));
        self.vocabulary = freq_vec.iter().enumerate().map(|(i, (token, _))| (token.clone(), i)).collect();
        self.feature_names = freq_vec.into_iter().map(|(token, _)| token).collect();
    }
    pub fn transform<T: AsRef<str> + Sync>(&self, texts: &[T]) -> Vec<Vec<f64>> {
        texts.par_iter().map(|text| {
            let mut counts = vec![0.0; self.vocabulary.len()];
            for token in text.as_ref().split_whitespace() { if let Some(&idx) = self.vocabulary.get(token) { counts[idx] += 1.0; } }
            counts
        }).collect()
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LabelEncoder { label_to_int: HashMap<String, usize>, int_to_label: Vec<String> }
impl LabelEncoder {
    pub fn new() -> Self { Self { label_to_int: HashMap::new(), int_to_label: Vec::new() } }
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
    pub fn transform<T: AsRef<str>>(&self, labels: &[T]) -> Vec<usize> { labels.iter().map(|label| *self.label_to_int.get(label.as_ref()).unwrap()).collect() }
}

fn check_class_balance(labels: &[Lineage]) {
    if labels.is_empty() { return; }
    let mut counts = HashMap::new();
    for label in labels { *counts.entry(label.inner()).or_insert(0) += 1; }
    let min_count = counts.values().min().unwrap_or(&0);
    let max_count = counts.values().max().unwrap_or(&0);
    if *min_count < 5 { warn!("Class imbalance: A class has only {} samples. Model may be unreliable.", min_count); }
    if *min_count > 0 && *max_count as f64 / *min_count as f64 > 10.0 { warn!("Class imbalance: Ratio of largest ({}) to smallest ({}) class is > 10:1.", max_count, min_count); }
}

fn kmerize(sequence: &str, k: usize) -> String {
    if sequence.len() < k { return String::new(); }
    (0..=sequence.len() - k).map(|i| &sequence[i..i + k]).collect::<Vec<&str>>().join(" ")
}

fn read_fasta(path: &str) -> Result<(Vec<GenomeSequence>, Vec<Lineage>)> {
    let file = File::open(path)?;
    let file_size = file.metadata()?.len();
    let pb = ProgressBar::new(file_size);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")?
        .progress_chars("#>-"));
    
    let reader = BufReader::new(pb.wrap_read(file));
    let mut sequences = Vec::new(); let mut lineages = Vec::new();
    let mut current_seq = String::new(); let mut current_lineage = String::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_lineage.is_empty() {
                sequences.push(GenomeSequence::new(current_seq.clone()).map_err(|e| anyhow!(e))?);
                lineages.push(Lineage::new(current_lineage.clone()).map_err(|e| anyhow!(e))?);
            }
            let header = line.trim_start_matches('>');
            let lineage_part = header.split('_').next().unwrap_or(header);
            current_lineage = lineage_part.to_string(); current_seq.clear();
        } else { current_seq.push_str(line.trim()); }
    }
    if !current_lineage.is_empty() {
        sequences.push(GenomeSequence::new(current_seq).map_err(|e| anyhow!(e))?);
        lineages.push(Lineage::new(current_lineage).map_err(|e| anyhow!(e))?);
    }
    pb.finish_with_message("FASTA reading complete");
    Ok((sequences, lineages))
}

fn load_input_data(args: &Args) -> Result<(Vec<GenomeSequence>, Vec<Lineage>)> {
    info!("▶ Reading input FASTA file: {}", args.input);
    let (sequences, lineages) = read_fasta(&args.input)?;
    check_class_balance(&lineages);
    Ok((sequences, lineages))
}

fn prepare_data(texts: &[String], labels: &[String], kmer_size: usize) -> Result<(CountVectorizer, LabelEncoder, Vec<Vec<f64>>, Vec<usize>)> {
    info!("▶ Generating k-mers and vectorizing features...");
    let kmer_texts: Vec<String> = texts.par_iter().map(|s| kmerize(s, kmer_size)).collect();
    let mut vectorizer = CountVectorizer::new();
    vectorizer.fit(&kmer_texts);
    let x_data = vectorizer.transform(&kmer_texts);
    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(labels);
    let y = label_encoder.transform(labels);
    if label_encoder.int_to_label.len() < 2 { return Err(anyhow!("Training data must contain at least two distinct classes.")); }
    Ok((vectorizer, label_encoder, x_data, y))
}

fn split_train_test(x_data: &[Vec<f64>], y: &[usize], test_ratio: f64) -> Result<(DenseMatrix<f64>, Vec<usize>, DenseMatrix<f64>, Vec<usize>)> {
    let n_samples = x_data.len();
    let raw_test_size = (n_samples as f64 * test_ratio).round() as usize;
    let test_size = if raw_test_size == 0 && n_samples > 1 { 1 } else { raw_test_size };
    if n_samples - test_size == 0 { return Err(anyhow!("Not enough samples for training after splitting.")); }
    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    indices.shuffle(&mut rng);
    let test_indices = &indices[..test_size]; let train_indices = &indices[test_size..];
    let train_data: Vec<Vec<f64>> = train_indices.iter().map(|&i| x_data[i].clone()).collect();
    let test_data: Vec<Vec<f64>> = test_indices.iter().map(|&i| x_data[i].clone()).collect();
    let y_train: Vec<usize> = train_indices.iter().map(|&i| y[i]).collect();
    let y_test: Vec<usize> = test_indices.iter().map(|&i| y[i]).collect();
    let x_train = DenseMatrix::from_2d_vec(&train_data)?;
    let x_test = DenseMatrix::from_2d_vec(&test_data)?;
    Ok((x_train, y_train, x_test, y_test))
}

fn save_model_bundle(bundle: &ModelBundle, output_path: &str) -> Result<()> {
    let file = File::create(output_path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    bincode::serialize_into(&mut encoder, bundle)?;
    encoder.finish()?;
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    if !(0.0..=1.0).contains(&args.test_split) { return Err(anyhow!("--test-split must be between 0.0 and 1.0.")); }
    if let Some(n) = args.threads { rayon::ThreadPoolBuilder::new().num_threads(n).build_global()?; }

    let overall_start = Instant::now();
    let (genome_vec, lineage_vec) = load_input_data(&args)?;
    let texts: Vec<String> = genome_vec.iter().map(|g| g.inner().to_string()).collect();
    let labels: Vec<String> = lineage_vec.iter().map(|l| l.inner().to_string()).collect();
    let (vectorizer, label_encoder, x_data, y) = prepare_data(&texts, &labels, args.kmer_size)?;
    let (x_train, y_train, x_test, y_test) = split_train_test(&x_data, &y, args.test_split)?;

    info!("▶ Starting to train the model -> {}", args.output);
    let train_start = Instant::now();
    let rf_params = RandomForestClassifierParameters {
        max_depth: None, min_samples_leaf: 1, min_samples_split: 2,
        n_trees: 100, m: Some((vectorizer.vocabulary.len() as f64).sqrt().floor() as usize),
        seed: 42, criterion: SplitCriterion::Gini, keep_samples: true,
    };
    let clf = RandomForestClassifier::fit(&x_train, &y_train, rf_params)?;
    info!("  Model training completed in {:.2} seconds.", train_start.elapsed().as_secs_f32());

    let y_pred = clf.predict(&x_test)?;
    let correct = y_pred.iter().zip(y_test.iter()).filter(|(&pred, &true_val)| pred == true_val).count();
    let accuracy = correct as f64 / y_test.len() as f64;
    info!("  Model's accuracy on the test set is: {:.2}%", accuracy * 100.0);

    let config = ModelConfig { kmer_size: args.kmer_size };
    let bundle = ModelBundle { config, vectorizer, label_encoder, model: clf };
    save_model_bundle(&bundle, &args.output)?;
    info!("  Saved the complete model bundle to {}", args.output);

    info!("✅ Process completed in {:.2} seconds.", overall_start.elapsed().as_secs_f32());
    Ok(())
}
