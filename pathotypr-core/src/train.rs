//! Train subcommand — Random Forest on k-mer feature-hashed sparse vectors.
//!
//! Pipeline: read labeled FASTA → feature-hash k-mers → train ensemble →
//! evaluate → save model + feature importance with genomic coordinates.

use crate::common::configure_thread_pool;
use crate::errors::{check_cancelled, AppError, AppResult, CancellationToken, ParallelCancellation};
use crate::fasta_io::read_fasta;
use crate::model::{
    LabelEncoder, ModelBundle, ModelConfig, DEFAULT_HASH_BUCKETS, MODEL_FORMAT_VERSION,
};
use crate::vectorizer::FeatureHasher;
use crate::sparse_tree::{ensemble_feature_importance, SparseDecisionTree, TreeParams};

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

// --- Module Constants ---
const DEFAULT_KMER_SIZE: usize = 21;
const DEFAULT_TEST_SPLIT: f64 = 0.2;
const N_TREES: u16 = 100;
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured report returned after a successful training run.
#[derive(Debug, Clone, Serialize)]
pub struct TrainReport {
    pub accuracy_pct: f64,
    pub n_training: usize,
    pub n_test: usize,
    pub n_classes: usize,
    pub n_features: usize,
    pub n_trees: u16,
    pub kmer_size: usize,
    pub model_size_mb: f64,
    pub training_time_secs: f32,
    pub class_names: Vec<String>,
    /// Per-fold accuracies when --cv-folds is used; empty otherwise.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub cv_fold_accuracies: Vec<f64>,
    /// Mean CV accuracy (only set when --cv-folds is used).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cv_mean_accuracy_pct: Option<f64>,
    /// Std dev of CV accuracy (only set when --cv-folds is used).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cv_std_accuracy_pct: Option<f64>,
    /// Out-of-bag accuracy from the final ensemble (always computed).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub oob_accuracy_pct: Option<f64>,
}

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
    /// Number of cross-validation folds for accuracy estimation (e.g., 5 or 10).
    /// When set, stratified k-fold CV is used instead of a single train/test split,
    /// giving more reliable accuracy estimates for small datasets.
    /// The final model is always trained on ALL data regardless of this setting.
    #[arg(long)]
    pub cv_folds: Option<usize>,
    /// Maximum tree depth. Lower values reduce overfitting on small datasets.
    /// Default: 20.
    #[arg(long, default_value_t = 20)]
    pub max_depth: usize,
    /// Minimum number of samples required in a leaf node.
    /// Higher values act as regularization, preventing the tree from memorizing
    /// individual training samples. Default: 5.
    #[arg(long, default_value_t = 5)]
    pub min_samples_leaf: usize,
    /// Cancellation token for stopping the task (GUI only, not CLI).
    #[arg(skip)]
    pub cancel_token: Option<CancellationToken>,
}

/// Prepares data: hashes k-mers into fixed-size buckets and encodes labels.
pub fn prepare_data(
    sequences: &[String],
    labels: &[String],
    kmer_size: usize,
    cancel_token: &Option<CancellationToken>,
) -> AppResult<(
    FeatureHasher,
    LabelEncoder,
    Vec<Vec<(usize, f32)>>,
    Vec<usize>,
)> {
    info!("▶ Vectorizing features with feature hashing...");

    check_cancelled(cancel_token)?;

    let hasher = FeatureHasher::new(DEFAULT_HASH_BUCKETS);
    info!(
        "  Feature hashing: {} buckets (2^{}).",
        hasher.num_buckets,
        (hasher.num_buckets as f64).log2() as u32
    );
    let x_data_sparse = hasher.transform_sparse(sequences, kmer_size);

    debug!(
        "  Transformed {} sequences into sparse feature vectors ({} buckets).",
        sequences.len(),
        hasher.num_features()
    );

    let mut label_encoder = LabelEncoder::new();
    label_encoder.fit(labels);
    let y_data = label_encoder.transform(labels)
        .map_err(|e| AppError::Parsing(e))?;
    if label_encoder.int_to_label.len() < 2 {
        return Err(AppError::NotEnoughData(
            "Training data must contain at least two distinct classes.".to_string(),
        ));
    }
    debug!(
        "  Labels encoded into {} classes.",
        label_encoder.int_to_label.len()
    );

    Ok((hasher, label_encoder, x_data_sparse, y_data))
}

/// Split data into train/test sets.
/// Kept for backward compatibility and unit tests; the main `run` function
/// now computes indices directly for both CV and single-split paths.
#[cfg(test)]
fn split_train_test<'a>(
    x_data: &'a [Vec<(usize, f32)>],
    y: &'a [usize],
    test_ratio: f64,
    seed: u64,
) -> AppResult<(
    Vec<&'a Vec<(usize, f32)>>,
    Vec<usize>,
    Vec<&'a Vec<(usize, f32)>>,
    Vec<usize>,
)> {
    let n_samples = x_data.len();
    if n_samples == 0 || y.is_empty() || n_samples != y.len() {
        return Err(AppError::NotEnoughData(
            "Training requires non-empty, aligned feature and label sets.".to_string(),
        ));
    }

    // Keep both train and test sets non-empty whenever possible.
    let mut test_size = (n_samples as f64 * test_ratio).round() as usize;
    if n_samples > 1 {
        if test_size == 0 && test_ratio > 0.0 {
            test_size = 1;
        }
        test_size = test_size.min(n_samples - 1);
    } else {
        test_size = 0;
    }

    let mut indices: Vec<usize> = (0..n_samples).collect();
    let mut rng = SmallRng::seed_from_u64(seed);
    indices.shuffle(&mut rng);

    let test_indices = &indices[..test_size];
    let train_indices = &indices[test_size..];
    if train_indices.is_empty() {
        return Err(AppError::NotEnoughData(
            "Training split produced an empty training set. Adjust --test-split.".to_string(),
        ));
    }

    let x_train: Vec<&Vec<(usize, f32)>> = train_indices.iter().map(|&i| &x_data[i]).collect();
    let y_train: Vec<usize> = train_indices.iter().map(|&i| y[i]).collect();

    let x_test: Vec<&Vec<(usize, f32)>> = test_indices.iter().map(|&i| &x_data[i]).collect();
    let y_test: Vec<usize> = test_indices.iter().map(|&i| y[i]).collect();

    info!(
        "  Data split into {} training samples and {} test samples.",
        x_train.len(),
        x_test.len()
    );
    Ok((x_train, y_train, x_test, y_test))
}

/// Wrapper to count bytes written through a `Write` impl.
struct CountingWriter<W: Write> {
    inner: W,
    count: usize,
}
impl<W: Write> CountingWriter<W> {
    fn new(inner: W) -> Self { Self { inner, count: 0 } }
    fn bytes_written(&self) -> usize { self.count }
}
impl<W: Write> Write for CountingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.count += n;
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> { self.inner.flush() }
}

/// Saves the model bundle using streaming serialization: bincode → zstd → file.
/// Unlike the previous approach (serialize to Vec then compress to Vec), this
/// never holds the full uncompressed + compressed buffers simultaneously,
/// reducing peak RAM by 50-200 MB for large models.
fn save_model_bundle(bundle: &ModelBundle, output_path: &str) -> AppResult<usize> {
    info!(
        "  Compressing and saving the model bundle to {}...",
        output_path
    );
    let file = BufWriter::new(File::create(output_path)?);
    let counting = CountingWriter::new(file);
    let mut encoder = zstd::Encoder::new(counting, 3).map_err(AppError::Io)?;
    bincode::serialize_into(&mut encoder, bundle).map_err(AppError::Serialization)?;
    let mut counting = encoder.finish().map_err(AppError::Io)?;
    // Explicitly flush the buffered writer so a failure on the final buffered
    // chunk (e.g. ENOSPC) surfaces as an error instead of being swallowed by
    // BufWriter's Drop, which would leave a silently truncated model file.
    counting.flush().map_err(AppError::Io)?;
    let compressed_size = counting.bytes_written();
    info!(
        "  Model bundle saved successfully ({:.2} MB).",
        compressed_size as f64 / 1_048_576.0
    );
    Ok(compressed_size)
}

/// Train an ensemble of decision trees on the given data.
///
/// Shared by both the main training path and each CV fold.
fn train_ensemble(
    x_data: &[Vec<(usize, f32)>],
    y_data: &[usize],
    sample_indices: &[usize],
    n_features: usize,
    num_classes: usize,
    n_trees: u16,
    max_depth: usize,
    min_samples_leaf: usize,
    cancel_token: &Option<CancellationToken>,
    show_progress: bool,
) -> AppResult<(Vec<SparseDecisionTree>, Vec<u64>)> {
    let n_samples = sample_indices.len();
    let max_features = (n_features as f64).sqrt().ceil() as usize;

    let params = TreeParams {
        max_depth,
        min_samples_leaf,
        max_features,
        n_classes: num_classes,
    };

    // Collect training rows for this subset
    let x_subset: Vec<Vec<(usize, f32)>> = sample_indices.iter().map(|&i| x_data[i].clone()).collect();
    let y_subset: Vec<usize> = sample_indices.iter().map(|&i| y_data[i]).collect();

    let mut seed_rng = SmallRng::seed_from_u64(42);
    let tree_seeds: Vec<u64> = (0..n_trees).map(|_| seed_rng.gen()).collect();
    let tree_cancellation = ParallelCancellation::new(cancel_token);

    let pb = if show_progress {
        let pb = ProgressBar::new(n_trees as u64);
        let style = ProgressStyle::default_bar()
            .template("  Training trees: [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap_or_else(|_| ProgressStyle::default_bar());
        pb.set_style(style);
        Some(Arc::new(pb))
    } else {
        None
    };

    let trees: Vec<SparseDecisionTree> = tree_seeds
        .par_iter()
        .map(|&seed| -> AppResult<SparseDecisionTree> {
            if tree_cancellation.is_cancelled() {
                return Err(AppError::Cancelled);
            }
            let mut rng = SmallRng::seed_from_u64(seed);
            let bootstrap: Vec<usize> = (0..n_samples)
                .map(|_| rng.gen_range(0..n_samples))
                .collect();
            let tree = SparseDecisionTree::fit(&x_subset, &y_subset, &bootstrap, &params, &mut rng);
            if let Some(ref pb) = pb {
                pb.inc(1);
            }
            Ok(tree)
        })
        .collect::<AppResult<Vec<_>>>()?;

    if let Some(pb) = pb {
        pb.finish_with_message("Tree training complete.");
    }
    Ok((trees, tree_seeds))
}

/// Compute out-of-bag (OOB) accuracy for a bagged ensemble.
///
/// Each tree was trained on a bootstrap sample (~63% of data). The remaining
/// ~37% are "out-of-bag" for that tree. For each sample, we aggregate votes
/// only from trees where it was OOB, giving a nearly unbiased accuracy
/// estimate without needing a held-out test set.
fn compute_oob_accuracy(
    trees: &[SparseDecisionTree],
    x_data: &[Vec<(usize, f32)>],
    y_data: &[usize],
    sample_indices: &[usize],
    num_classes: usize,
    tree_seeds: &[u64],
) -> Option<f64> {
    let n_samples = sample_indices.len();
    if n_samples == 0 || trees.is_empty() {
        return None;
    }

    // For each sample (local index), accumulate votes from OOB trees
    let mut oob_votes: Vec<Vec<u32>> = vec![vec![0u32; num_classes]; n_samples];
    let mut oob_count: Vec<u32> = vec![0; n_samples];

    for (tree_idx, seed) in tree_seeds.iter().enumerate() {
        // Regenerate the bootstrap for this tree
        let mut rng = SmallRng::seed_from_u64(*seed);
        let mut in_bag = vec![false; n_samples];
        for _ in 0..n_samples {
            let idx = rng.gen_range(0..n_samples);
            in_bag[idx] = true;
        }

        // Vote only for OOB samples
        for local_idx in 0..n_samples {
            if !in_bag[local_idx] {
                let global_idx = sample_indices[local_idx];
                let pred = trees[tree_idx].predict_one(&x_data[global_idx]);
                if pred < num_classes {
                    oob_votes[local_idx][pred] += 1;
                }
                oob_count[local_idx] += 1;
            }
        }
    }

    // Count correct predictions (only for samples that appeared OOB at least once)
    let mut correct = 0usize;
    let mut evaluated = 0usize;
    for local_idx in 0..n_samples {
        if oob_count[local_idx] == 0 {
            continue; // sample was in every bootstrap (very unlikely with 100 trees)
        }
        let predicted = oob_votes[local_idx]
            .iter()
            .enumerate()
            .max_by_key(|(_, &c)| c)
            .map(|(cls, _)| cls)
            .unwrap_or(0);
        if predicted == y_data[sample_indices[local_idx]] {
            correct += 1;
        }
        evaluated += 1;
    }

    if evaluated == 0 {
        return None;
    }
    Some((correct as f64 / evaluated as f64) * 100.0)
}

/// Evaluate an ensemble on a test set, returning accuracy percentage.
fn evaluate_ensemble(
    trees: &[SparseDecisionTree],
    x_data: &[Vec<(usize, f32)>],
    y_data: &[usize],
    test_indices: &[usize],
    num_classes: usize,
) -> f64 {
    if test_indices.is_empty() {
        return 0.0;
    }
    let mut correct = 0usize;
    for &idx in test_indices {
        let row = &x_data[idx];
        let mut vote_counts = vec![0u16; num_classes];
        for tree in trees {
            let pred = tree.predict_one(row);
            if pred < num_classes {
                vote_counts[pred] += 1;
            }
        }
        let best_class = vote_counts
            .iter()
            .enumerate()
            .max_by_key(|(_, &c)| c)
            .map(|(idx, _)| idx)
            .unwrap_or(0);
        if best_class == y_data[idx] {
            correct += 1;
        }
    }
    (correct as f64 / test_indices.len() as f64) * 100.0
}

/// Stratified k-fold cross-validation.
///
/// Splits data into k folds preserving class proportions. For each fold,
/// trains an ensemble on k-1 folds and evaluates on the held-out fold.
/// Returns per-fold accuracies.
fn stratified_cv(
    x_data: &[Vec<(usize, f32)>],
    y_data: &[usize],
    k: usize,
    n_features: usize,
    num_classes: usize,
    n_trees: u16,
    max_depth: usize,
    min_samples_leaf: usize,
    cancel_token: &Option<CancellationToken>,
) -> AppResult<Vec<f64>> {
    let n = x_data.len();

    // Group indices by class
    let mut class_indices: Vec<Vec<usize>> = vec![Vec::new(); num_classes];
    for (i, &label) in y_data.iter().enumerate() {
        if label < num_classes {
            class_indices[label].push(i);
        }
    }

    // Shuffle within each class for randomness, then distribute into folds round-robin
    let mut rng = SmallRng::seed_from_u64(42);
    let mut fold_assignments = vec![0usize; n];
    for class_group in &mut class_indices {
        class_group.shuffle(&mut rng);
        for (i, &idx) in class_group.iter().enumerate() {
            fold_assignments[idx] = i % k;
        }
    }

    let mut fold_accuracies = Vec::with_capacity(k);

    for fold in 0..k {
        check_cancelled(cancel_token)?;
        info!("  📊 CV fold {}/{}", fold + 1, k);

        let train_idx: Vec<usize> = (0..n).filter(|&i| fold_assignments[i] != fold).collect();
        let test_idx: Vec<usize> = (0..n).filter(|&i| fold_assignments[i] == fold).collect();

        if train_idx.is_empty() || test_idx.is_empty() {
            warn!("  Fold {} has empty train or test set, skipping.", fold + 1);
            continue;
        }

        let (trees, _seeds) = train_ensemble(
            x_data, y_data, &train_idx,
            n_features, num_classes, n_trees,
            max_depth, min_samples_leaf,
            cancel_token, false,
        )?;

        let acc = evaluate_ensemble(&trees, x_data, y_data, &test_idx, num_classes);
        info!("  Fold {} accuracy: {:.2}%  (train={}, test={})", fold + 1, acc, train_idx.len(), test_idx.len());
        fold_accuracies.push(acc);
    }

    Ok(fold_accuracies)
}

/// Main execution logic for the `train` subcommand.
pub fn run(args: Args) -> AppResult<TrainReport> {
    if !(0.0..1.0).contains(&args.test_split) {
        return Err(AppError::Generic(
            "--test-split must be between 0.0 and 1.0.".to_string(),
        ));
    }
    crate::common::validate_kmer_size(args.kmer_size)?;
    if let Some(k) = args.cv_folds {
        if k < 2 {
            return Err(AppError::Generic(
                "--cv-folds must be at least 2.".to_string(),
            ));
        }
    }
    configure_thread_pool(args.threads);

    let cancel_token = &args.cancel_token;
    check_cancelled(cancel_token)?;

    let overall_start = Instant::now();

    info!("▶ Loading and preparing data...");
    let (sequences, lineages, headers) = read_fasta(&args.input)?;
    if sequences.len() < 10 {
        warn!(
            "Input contains very few sequences ({}), the trained model may not be reliable.",
            sequences.len()
        );
    }

    check_cancelled(cancel_token)?;

    let (vectorizer, label_encoder, x_data_sparse, y_data) =
        prepare_data(&sequences, &lineages, args.kmer_size, cancel_token)?;

    check_cancelled(cancel_token)?;

    let n_features = vectorizer.num_features();
    let num_classes = label_encoder.int_to_label.len();
    let n_total = x_data_sparse.len();

    // --- Accuracy estimation: CV or single split ---
    let mut accuracy_pct = 0.0_f64;
    let n_training;
    let n_test;
    let mut cv_fold_accuracies: Vec<f64> = Vec::new();
    let mut cv_mean_accuracy_pct: Option<f64> = None;
    let mut cv_std_accuracy_pct: Option<f64> = None;

    if let Some(k) = args.cv_folds {
        // --- Stratified k-fold cross-validation ---
        if n_total < k {
            return Err(AppError::NotEnoughData(format!(
                "Cannot perform {}-fold CV with only {} samples. Reduce --cv-folds or add more data.",
                k, n_total
            )));
        }
        info!("▶ Running {}-fold stratified cross-validation...", k);
        let fold_accs = stratified_cv(
            &x_data_sparse, &y_data, k,
            n_features, num_classes, N_TREES,
            args.max_depth, args.min_samples_leaf,
            cancel_token,
        )?;

        if !fold_accs.is_empty() {
            let mean = fold_accs.iter().sum::<f64>() / fold_accs.len() as f64;
            // Bessel's correction: divide by (n-1) for sample standard deviation
            let n = fold_accs.len() as f64;
            let variance = fold_accs.iter().map(|a| (a - mean).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);
            let std_dev = variance.sqrt();

            info!("  ✅ {}-fold CV accuracy: {:.2}% ± {:.2}%", k, mean, std_dev);
            accuracy_pct = mean;
            cv_mean_accuracy_pct = Some(mean);
            cv_std_accuracy_pct = Some(std_dev);
            cv_fold_accuracies = fold_accs;
        }

        // Final model trained on ALL data (below)
        n_training = n_total;
        n_test = 0usize;
    } else {
        // --- Single train/test split (legacy behavior) ---
        // Reproduce the same shuffle as split_train_test with seed 42
        let mut all_indices: Vec<usize> = (0..n_total).collect();
        let mut rng = SmallRng::seed_from_u64(42);
        all_indices.shuffle(&mut rng);

        let mut test_size = (n_total as f64 * args.test_split).round() as usize;
        if n_total > 1 {
            if test_size == 0 && args.test_split > 0.0 { test_size = 1; }
            test_size = test_size.min(n_total - 1);
        } else {
            test_size = 0;
        }

        let test_indices = &all_indices[..test_size];
        let train_indices = &all_indices[test_size..];
        n_training = train_indices.len();
        n_test = test_size;

        if !test_indices.is_empty() {
            info!("▶ Training evaluation ensemble ({} train, {} test)...", n_training, n_test);
            let (eval_trees, _seeds) = train_ensemble(
                &x_data_sparse, &y_data, &train_indices.to_vec(),
                n_features, num_classes, N_TREES,
                args.max_depth, args.min_samples_leaf,
                cancel_token, false,
            )?;
            accuracy_pct = evaluate_ensemble(&eval_trees, &x_data_sparse, &y_data, &test_indices.to_vec(), num_classes);
            info!("  ✅ Test Set Accuracy: {:.2}%", accuracy_pct);
        } else {
            info!("▶ Test set is empty, skipping evaluation.");
        }
    }

    check_cancelled(cancel_token)?;

    // --- Train FINAL model on ALL data ---
    info!("▶ Training final model on all {} samples ({} trees)...", n_total, N_TREES);
    let train_start = Instant::now();
    let all_indices: Vec<usize> = (0..n_total).collect();
    info!("  Tree params: max_depth={}, min_samples_leaf={}", args.max_depth, args.min_samples_leaf);
    let (trees, tree_seeds) = train_ensemble(
        &x_data_sparse, &y_data, &all_indices,
        n_features, num_classes, N_TREES,
        args.max_depth, args.min_samples_leaf,
        cancel_token, true,
    )?;
    info!(
        "  Final model trained in {:.2} seconds ({} trees on {} samples).",
        train_start.elapsed().as_secs_f32(),
        trees.len(),
        n_total,
    );

    // --- OOB accuracy (free estimate from bootstrap) ---
    // Uses the exact seeds from train_ensemble to reproduce bootstrap samples.
    let oob_accuracy_pct = compute_oob_accuracy(
        &trees, &x_data_sparse, &y_data, &all_indices,
        num_classes, &tree_seeds,
    );
    if let Some(oob) = oob_accuracy_pct {
        info!("  🎯 Out-of-bag accuracy: {:.2}%", oob);
    }

    info!("▶ Finalizing and saving the model...");
    let class_names = label_encoder.int_to_label.clone();
    let config = ModelConfig {
        pathotypr_version: VERSION.to_string(),
        kmer_size: args.kmer_size,
        n_trees: N_TREES,
        format_version: MODEL_FORMAT_VERSION,
    };
    let bundle = ModelBundle {
        config,
        vectorizer,
        label_encoder,
        trees,
    };

    let model_bytes = save_model_bundle(&bundle, &args.output)?;

    // --- Feature Importance Analysis ---
    check_cancelled(cancel_token)?;
    info!("▶ Computing feature importance and resolving discriminant k-mers...");
    let top_n = 500;
    let importance = ensemble_feature_importance(&bundle.trees);
    let top_buckets: FxHashSet<usize> = importance
        .iter()
        .take(top_n)
        .map(|&(feat, _)| feat as usize)
        .collect();

    let bucket_kmers = bundle
        .vectorizer
        .reverse_map_buckets(&sequences, args.kmer_size, &top_buckets);

    // Write importance TSV alongside the model file
    let importance_path = Path::new(&args.output)
        .with_extension("importance.tsv");
    let total_splits: u32 = importance.iter().map(|&(_, c)| c).sum();
    {
        let mut out = BufWriter::new(File::create(&importance_path)?);
        writeln!(out, "rank\tbucket\tsplit_count\timportance_pct\tkmers")?;
        for (rank, &(feat, count)) in importance.iter().take(top_n).enumerate() {
            let pct = if total_splits > 0 {
                count as f64 / total_splits as f64 * 100.0
            } else {
                0.0
            };
            let kmers_str = bucket_kmers
                .get(&(feat as usize))
                .map(|v| v.join(","))
                .unwrap_or_default();
            writeln!(out, "{}\t{}\t{}\t{:.4}\t{}", rank + 1, feat, count, pct, kmers_str)?;
        }
        // Flush explicitly: BufWriter's Drop ignores the result, so a failure on
        // the final buffered chunk would be swallowed and a truncated report
        // reported as written.
        out.flush()?;
    }
    info!(
        "  Feature importance written to {} ({} top features).",
        importance_path.display(),
        top_n.min(importance.len())
    );

    // --- Genomic Coordinates for discriminant k-mers ---
    let coords = bundle
        .vectorizer
        .reverse_map_with_coords(&sequences, args.kmer_size, &top_buckets);

    // Build a lookup: bucket -> (rank, split_count, importance_pct)
    let bucket_rank: FxHashMap<usize, (usize, u32, f64)> = importance
        .iter()
        .take(top_n)
        .enumerate()
        .map(|(rank, &(feat, count))| {
            let pct = if total_splits > 0 {
                count as f64 / total_splits as f64 * 100.0
            } else {
                0.0
            };
            (feat as usize, (rank + 1, count, pct))
        })
        .collect();

    let coords_path = Path::new(&args.output).with_extension("importance.coords.tsv");
    let mut sorted_coords = coords;
    sorted_coords.sort_by(|a, b| {
        let ra = bucket_rank.get(&a.bucket).map(|r| r.0).unwrap_or(usize::MAX);
        let rb = bucket_rank.get(&b.bucket).map(|r| r.0).unwrap_or(usize::MAX);
        ra.cmp(&rb)
            .then(a.seq_index.cmp(&b.seq_index))
            .then(a.position.cmp(&b.position))
    });
    {
        let mut out = BufWriter::new(File::create(&coords_path)?);
        writeln!(out, "rank\tbucket\tsplit_count\timportance_pct\tkmer\tsequence\tlineage\tposition")?;
        for hit in &sorted_coords {
            let (rank, splits, pct) = bucket_rank
                .get(&hit.bucket)
                .copied()
                .unwrap_or((0, 0, 0.0));
            let header = headers.get(hit.seq_index).map(|s| s.as_str()).unwrap_or("?");
            let lineage = lineages.get(hit.seq_index).map(|s| s.as_str()).unwrap_or("?");
            writeln!(
                out,
                "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
                rank, hit.bucket, splits, pct, hit.kmer, header, lineage, hit.position
            )?;
        }
        // See above: an unflushed BufWriter can silently truncate the report.
        out.flush()?;
    }
    info!(
        "  Genomic coordinates written to {} ({} hits).",
        coords_path.display(),
        sorted_coords.len()
    );

    let training_time_secs = overall_start.elapsed().as_secs_f32();

    info!(
        "🏁 Process completed in {:.2} seconds.",
        training_time_secs
    );
    Ok(TrainReport {
        accuracy_pct,
        n_training,
        n_test,
        n_classes: num_classes,
        n_features,
        n_trees: N_TREES,
        kmer_size: args.kmer_size,
        model_size_mb: model_bytes as f64 / 1_048_576.0,
        training_time_secs,
        class_names,
        cv_fold_accuracies,
        cv_mean_accuracy_pct,
        cv_std_accuracy_pct,
        oob_accuracy_pct,
    })
}

#[cfg(test)]
mod tests {
    use super::split_train_test;

    #[test]
    fn split_keeps_training_set_non_empty() {
        let x = vec![vec![(0usize, 1.0f32)], vec![(0usize, 2.0f32)]];
        let y = vec![0usize, 1usize];
        let (x_train, y_train, x_test, y_test) = split_train_test(&x, &y, 0.99, 42).unwrap();
        assert!(!x_train.is_empty());
        assert_eq!(x_train.len(), y_train.len());
        assert_eq!(x_test.len(), y_test.len());
    }

    #[test]
    fn split_rejects_misaligned_inputs() {
        let x = vec![vec![(0usize, 1.0f32)]];
        let y = vec![];
        assert!(split_train_test(&x, &y, 0.2, 42).is_err());
    }
}
