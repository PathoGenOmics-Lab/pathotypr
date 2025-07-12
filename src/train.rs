//! This module handles the `train` subcommand.
//!
//! It has been refactored for high performance and low memory usage by:
//! 1. Using `needletail` for efficient, allocation-free k-merization.
//! 2. Representing the main feature matrix in a sparse format to avoid OOM errors.
//! 3. Using `rayon` for parallel k-mer counting and `FxHashMap` for speed.
//! 4. Compressing the final model bundle with `zstd`.

// --- Standard and External Crates ---
use crate::common::{LabelEncoder, ModelBundle, ModelConfig, CountVectorizer};
use crate::errors::{AppError, AppResult};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use needletail::Sequence;
use rand::seq::SliceRandom;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
// MODIFIED: All matrix operations now use f32.
use smartcore::linalg::basic::matrix::DenseMatrix;
use smartcore::tree::decision_tree_classifier::{
    DecisionTreeClassifier, DecisionTreeClassifierParameters,
};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

// --- Module Constants ---
const DEFAULT_KMER_SIZE: usize = 21;
const DEFAULT_TEST_SPLIT: f64 = 0.2;
const N_TREES: u16 = 100;
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Command-line arguments for the `train` subcommand.
#[derive(Parser, Debug)]
pub struct Args {
    #[arg(short = 'i', long)]
    pub input: String,
    #[arg(short = 'o', long)]
    pub output: String,
    #[arg(short = 'k', long, default_value_t = DEFAULT_KMER_SIZE)]
    pub kmer_size: usize,
    #[arg(short = 's', long, default_value_t = DEFAULT_TEST_SPLIT)]
    pub test_split: f64,
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
}

/// Reads sequences and their corresponding lineages from a FASTA file.
fn read_fasta(path: &str) -> AppResult<(Vec<String>, Vec<String>)> {
    let file = File::open(path)?;
    let pb = ProgressBar::new(file.metadata()?.len());
    let style = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")?
        .progress_chars("#>-");
    pb.set_style(style);

    let reader = BufReader::new(pb.wrap_read(file));
    let mut sequences = Vec::new();
    let mut lineages = Vec::new();
    let mut current_seq = String::new();
    let mut current_lineage = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_lineage.is_empty() {
                sequences.push(current_seq.clone());
                lineages.push(current_lineage.clone());
            }
            let header = line.trim_start_matches('>');
            current_lineage = header.split('_').next().unwrap_or(header).to_string();
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_lineage.is_empty() {
        sequences.push(current_seq);
        lineages.push(current_lineage);
    }
    pb.finish_with_message("FASTA reading complete");
    Ok((sequences, lineages))
}

/// Performs parallel k-mer counting on all sequences using a map-reduce strategy.
fn parallel_kmer_counting(sequences: &[String], k: usize) -> FxHashMap<u64, u32> {
    info!("🧬 Counting k-mers in parallel...");
    let pb = ProgressBar::new(sequences.len() as u64);
    let style = ProgressStyle::default_bar().template("  {spinner:.blue} Counting k-mers: [{bar:40.magenta/purple}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("##-");
    pb.set_style(style);

    let final_counts = sequences
        .par_iter()
        .map(|seq| {
            let mut local_counts = FxHashMap::default();
            for (_, bitkmer_tuple, _) in seq.as_bytes().bit_kmers(k as u8, true) {
                *local_counts.entry(bitkmer_tuple.0).or_insert(0) += 1;
            }
            pb.inc(1);
            local_counts
        })
        .reduce(FxHashMap::default, |mut a, b| {
            for (kmer, count) in b {
                *a.entry(kmer).or_insert(0) += count;
            }
            a
        });
        
    pb.finish_with_message("K-mer counting finished.");
    info!("  Total unique k-mers found: {}", final_counts.len());
    final_counts
}

/// MODIFIED: Prepares data, returning a sparse representation `Vec<Vec<(usize, f32)>>`.
fn prepare_data(
    sequences: &[String],
    labels: &[String],
    kmer_size: usize,
) -> AppResult<(CountVectorizer, LabelEncoder, Vec<Vec<(usize, f32)>>, Vec<usize>)> {
    info!("▶ Vectorizing features and encoding labels...");

    let all_kmer_counts = parallel_kmer_counting(sequences, kmer_size);
    if all_kmer_counts.is_empty() {
        return Err(AppError::NotEnoughData("No k-mers could be generated from the input sequences.".to_string()));
    }

    let mut vectorizer = CountVectorizer::new();
    vectorizer.fit(&all_kmer_counts);
    debug!("  Vectorizer fitted with a vocabulary of {} features.", vectorizer.num_features);

    // MODIFIED: Use the new sparse transform method.
    let x_data_sparse = vectorizer.transform_sparse(sequences, kmer_size);
    debug!("  Transformed {} sequences into sparse feature vectors.", sequences.len());

    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(labels);
    let y_data = label_encoder.transform(labels);
    if label_encoder.int_to_label.len() < 2 {
        return Err(AppError::NotEnoughData(
            "Training data must contain at least two distinct classes.".to_string(),
        ));
    }
    debug!("  Labels encoded into {} classes.", label_encoder.int_to_label.len());
    
    Ok((vectorizer, label_encoder, x_data_sparse, y_data))
}

/// MODIFIED: This function now takes sparse f32 data as input.
fn split_train_test<'a>(
    x_data: &'a [Vec<(usize, f32)>],
    y: &'a [usize],
    test_ratio: f64,
    seed: u64,
) -> (
    Vec<&'a Vec<(usize, f32)>>,
    Vec<usize>,
    Vec<&'a Vec<(usize, f32)>>,
    Vec<usize>,
) {
    let n_samples = x_data.len();
    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = SmallRng::seed_from_u64(seed);
    indices.shuffle(&mut rng);

    let test_size = (n_samples as f64 * test_ratio).round() as usize;
    
    let test_indices = &indices[..test_size];
    let train_indices = &indices[test_size..];

    let x_train: Vec<&Vec<(usize, f32)>> = train_indices.iter().map(|&i| &x_data[i]).collect();
    let y_train: Vec<usize> = train_indices.iter().map(|&i| y[i]).collect();
    
    let x_test: Vec<&Vec<(usize, f32)>> = test_indices.iter().map(|&i| &x_data[i]).collect();
    let y_test: Vec<usize> = test_indices.iter().map(|&i| y[i]).collect();

    info!("  Data split into {} training samples and {} test samples.", x_train.len(), x_test.len());
    (x_train, y_train, x_test, y_test)
}

/// Saves the complete model bundle using `bincode` for serialization and `zstd` for compression.
fn save_model_bundle(bundle: &ModelBundle, output_path: &str) -> AppResult<()> {
    info!("  Compressing and saving the model bundle to {}...", output_path);
    let serialized = bincode::serialize(bundle).map_err(|e| AppError::Serialization(e))?;
    let compressed = zstd::encode_all(&serialized[..], 3).map_err(|e| AppError::Io(e))?;
    
    let mut file = File::create(output_path)?;
    file.write_all(&compressed)?;
    info!("  Model bundle saved successfully ({:.2} MB).", compressed.len() as f64 / 1_048_576.0);
    Ok(())
}

/// Main execution logic for the `train` subcommand.
pub fn run(args: Args) -> AppResult<()> {
    if !(0.0..1.0).contains(&args.test_split) {
        return Err(AppError::Generic("--test-split must be between 0.0 and 1.0.".to_string()));
    }
    if args.kmer_size == 0 || args.kmer_size > 31 {
         return Err(AppError::Generic("--kmer-size must be between 1 and 31 to fit in a u64.".to_string()));
    }
    if let Some(n) = args.threads {
        info!("▶ Setting number of threads to: {}", n);
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global()
            .map_err(|e| AppError::Generic(format!("Failed to build thread pool: {}", e)))?;
    }

    let overall_start = Instant::now();

    info!("▶ Loading and preparing data...");
    let (sequences, lineages) = read_fasta(&args.input)?;
    if sequences.len() < 10 {
        warn!("Input contains very few sequences ({}), the trained model may not be reliable.", sequences.len());
    }
    
    let (vectorizer, label_encoder, x_data_sparse, y_data) =
        prepare_data(&sequences, &lineages, args.kmer_size)?;
    
    let (x_train_sparse, y_train, x_test_sparse, y_test) =
        split_train_test(&x_data_sparse, &y_data, args.test_split, 42);

    info!("▶ Training an ensemble of {} decision trees...", N_TREES);
    let train_start = Instant::now();
    let mut trees = Vec::with_capacity(N_TREES as usize);
    let mut rng = SmallRng::seed_from_u64(42);
    let n_samples = x_train_sparse.len();
    let n_features = vectorizer.num_features;

    let pb = ProgressBar::new(N_TREES as u64);
    let style = ProgressStyle::default_bar().template("  Training trees: [{bar:40.cyan/blue}] {pos}/{len} ({eta})")?.progress_chars("=>-");
    pb.set_style(style);

    for _ in 0..N_TREES {
        // Bootstrap aggregating (bagging): sample with replacement.
        let sample_indices: Vec<usize> = (0..n_samples).map(|_| rng.gen_range(0..n_samples)).collect();
        
        // MODIFIED: Create a small DENSE matrix for this sample only, using f32.
        let mut x_sample_dense_vec: Vec<Vec<f32>> = Vec::with_capacity(n_samples);
        let mut y_sample: Vec<usize> = Vec::with_capacity(n_samples);

        for &idx in &sample_indices {
            let sparse_row = x_train_sparse[idx];
            let mut dense_row = vec![0.0_f32; n_features];
            for &(feature_idx, val) in sparse_row {
                dense_row[feature_idx] = val;
            }
            x_sample_dense_vec.push(dense_row);
            y_sample.push(y_train[idx]);
        }
        
        let x_sample = DenseMatrix::from_2d_vec(&x_sample_dense_vec)?;
        // MODIFIED: Explicitly drop the large vector to free memory immediately.
        drop(x_sample_dense_vec);

        let tree = DecisionTreeClassifier::fit(&x_sample, &y_sample, DecisionTreeClassifierParameters::default())?;
        trees.push(tree);
        pb.inc(1);
    }
    pb.finish_with_message("Tree training complete.");
    info!("  Model training finished in {:.2} seconds.", train_start.elapsed().as_secs_f32());

    // MODIFIED: Evaluation on the sparse test set.
    if !y_test.is_empty() {
        info!("▶ Evaluating model on the test set...");
        let mut y_pred_votes = Vec::with_capacity(y_test.len());

        for sparse_row in x_test_sparse {
            // Create a single dense row for prediction.
            let mut dense_row = vec![0.0_f32; n_features];
            for &(feature_idx, val) in sparse_row {
                dense_row[feature_idx] = val;
            }
            let x_row = DenseMatrix::from_2d_vec(&vec![dense_row])?;

            let mut votes = FxHashMap::default();
            for tree in &trees {
                let prediction = tree.predict(&x_row)?[0];
                *votes.entry(prediction).or_insert(0) += 1;
            }
            let best_class = votes.into_iter().max_by_key(|&(_, count)| count).map(|(class, _)| class).unwrap_or(0);
            y_pred_votes.push(best_class);
        }

        let correct = y_pred_votes.iter().zip(y_test.iter()).filter(|(&pred, &true_val)| pred == true_val).count();
        let accuracy = correct as f64 / y_test.len() as f64;
        info!("  ✅ Test Set Accuracy: {:.2}%", accuracy * 100.0);
    } else {
        info!("▶ Test set is empty, skipping evaluation.");
    }

    info!("▶ Finalizing and saving the model...");
    let config = ModelConfig {
        pathotypr_version: VERSION.to_string(),
        kmer_size: args.kmer_size,
        n_trees: N_TREES,
    };
    let bundle = ModelBundle {
        config,
        vectorizer,
        label_encoder,
        trees,
    };

    save_model_bundle(&bundle, &args.output)?;
    
    info!("🏁 Process completed in {:.2} seconds.", overall_start.elapsed().as_secs_f32());
    Ok(())
}