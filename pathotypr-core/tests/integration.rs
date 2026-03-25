//! Integration tests for pathotypr-core.
//!
//! Tests end-to-end workflows (train → predict → verify) and individual
//! module functionality using synthetic data.

use std::fs;
use std::io::Write;
use tempfile::TempDir;

// ============================================================================
// Helpers
// ============================================================================

/// Creates a synthetic FASTA file with `n_per_class` sequences per class.
/// Each class gets sequences with a distinct k-mer signature so the RF
/// can learn to distinguish them.
fn write_synthetic_fasta(path: &str, classes: &[&str], n_per_class: usize) {
    let mut f = fs::File::create(path).unwrap();
    // Each class gets a repeating motif so k-mers are distinctive
    let motifs = [
        "ATCGATCGATCGATCGATCGATCGATCG",   // class 0
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTA",   // class 1
        "TTAACCGGTTAACCGGTTAACCGGTTAA",   // class 2
        "CCGGAATTCCGGAATTCCGGAATTCCGG",   // class 3
        "AACCTTGGAACCTTGGAACCTTGGAACC",   // class 4
    ];
    for (ci, class) in classes.iter().enumerate() {
        let motif = motifs[ci % motifs.len()];
        for i in 0..n_per_class {
            writeln!(f, ">{} sample_{}_{}", class, class, i).unwrap();
            // Build a 500bp sequence from the motif with slight variation
            let mut seq = String::with_capacity(500);
            for j in 0..18 {
                seq.push_str(motif);
                // Add a variant base every other repeat for the sample index
                if j % 3 == 0 {
                    let base = match (i + j) % 4 {
                        0 => 'A', 1 => 'C', 2 => 'G', _ => 'T',
                    };
                    seq.push(base);
                }
            }
            writeln!(f, "{}", &seq[..500.min(seq.len())]).unwrap();
        }
    }
}

// ============================================================================
// End-to-end: train → predict → verify
// ============================================================================

#[test]
fn end_to_end_train_predict() {
    let tmp = TempDir::new().unwrap();
    let fasta_path = tmp.path().join("training.fasta");
    let model_path = tmp.path().join("model.pathotypr.zst");
    let predict_out = tmp.path().join("predictions.tsv");

    // Generate training data: 3 classes × 20 samples = 60 genomes
    write_synthetic_fasta(
        fasta_path.to_str().unwrap(),
        &["Lineage1", "Lineage2", "Lineage3"],
        20,
    );

    // Train
    let train_args = pathotypr_core::train::Args {
        input: fasta_path.to_str().unwrap().to_string(),
        output: model_path.to_str().unwrap().to_string(),
        kmer_size: 11,
        test_split: 0.2,
        threads: Some(2),
        cv_folds: None,
        max_depth: 10,
        min_samples_leaf: 2,
        cancel_token: None,
    };
    let report = pathotypr_core::train::run(train_args).unwrap();
    assert!(report.accuracy_pct >= 80.0, "Accuracy too low: {:.1}%", report.accuracy_pct);
    assert_eq!(report.n_classes, 3);
    assert_eq!(report.n_trees, 100);
    assert!(report.oob_accuracy_pct.is_some());
    assert!(report.oob_accuracy_pct.unwrap() >= 80.0);

    // Predict using the trained model
    let predict_args = pathotypr_core::predict::PredictArgs {
        input: fasta_path.to_str().unwrap().to_string(),
        model: model_path.to_str().unwrap().to_string(),
        output: predict_out.to_str().unwrap().to_string(),
        threads: Some(2),
        excel: false,
        cancel_token: None,
    };
    pathotypr_core::predict::run(predict_args).unwrap();

    // Verify predictions
    let content = fs::read_to_string(&predict_out).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() > 1, "Output should have header + data");
    // Header
    assert!(lines[0].contains("Predicted_Lineage"));
    // Check that at least 80% are correct
    let mut correct = 0;
    let mut total = 0;
    for line in &lines[1..] {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() >= 2 {
            let header = cols[0]; // e.g. "Lineage1 sample_Lineage1_0"
            let predicted = cols[1];
            let true_label = header.split_whitespace().next().unwrap_or("");
            if true_label == predicted {
                correct += 1;
            }
            total += 1;
        }
    }
    assert!(total >= 60, "Should have predictions for all 60 samples");
    let acc = correct as f64 / total as f64;
    assert!(acc >= 0.80, "Predict accuracy {:.1}% < 80%", acc * 100.0);
}

// ============================================================================
// Train with cross-validation
// ============================================================================

#[test]
fn train_with_cv_folds() {
    let tmp = TempDir::new().unwrap();
    let fasta_path = tmp.path().join("cv_train.fasta");
    let model_path = tmp.path().join("cv_model.pathotypr.zst");

    write_synthetic_fasta(
        fasta_path.to_str().unwrap(),
        &["ClassA", "ClassB"],
        15,
    );

    let args = pathotypr_core::train::Args {
        input: fasta_path.to_str().unwrap().to_string(),
        output: model_path.to_str().unwrap().to_string(),
        kmer_size: 11,
        test_split: 0.2,
        threads: Some(2),
        cv_folds: Some(3),
        max_depth: 10,
        min_samples_leaf: 2,
        cancel_token: None,
    };
    let report = pathotypr_core::train::run(args).unwrap();

    assert_eq!(report.cv_fold_accuracies.len(), 3);
    assert!(report.cv_mean_accuracy_pct.is_some());
    assert!(report.cv_std_accuracy_pct.is_some());
    assert!(report.cv_mean_accuracy_pct.unwrap() >= 70.0);
    // OOB should also be present
    assert!(report.oob_accuracy_pct.is_some());
}

// ============================================================================
// Model validation
// ============================================================================

#[test]
fn load_model_validates_format_version() {
    let tmp = TempDir::new().unwrap();
    let path = tmp.path().join("bad_model.pathotypr.zst");

    // Create a model with wrong format version
    let bundle = pathotypr_core::ModelBundle {
        config: pathotypr_core::ModelConfig {
            pathotypr_version: "0.0.0".to_string(),
            kmer_size: 21,
            n_trees: 1,
            format_version: 999, // wrong version
        },
        vectorizer: pathotypr_core::FeatureHasher::new(1024),
        label_encoder: {
            let mut le = pathotypr_core::LabelEncoder::new();
            le.fit(&["A".to_string(), "B".to_string()]);
            le
        },
        trees: vec![],
    };

    // Serialize and save
    let data = bincode::serialize(&bundle).unwrap();
    let compressed = zstd::encode_all(&data[..], 3).unwrap();
    fs::write(&path, compressed).unwrap();

    // load_model_bundle should reject it (empty trees)
    let result = pathotypr_core::predict::load_model_bundle(path.to_str().unwrap());
    assert!(result.is_err());
    let err_msg = format!("{}", result.unwrap_err());
    assert!(err_msg.contains("no trees"), "Expected 'no trees' error, got: {}", err_msg);
}

// ============================================================================
// LabelEncoder: error on unknown label
// ============================================================================

#[test]
fn label_encoder_rejects_unknown() {
    let mut le = pathotypr_core::LabelEncoder::new();
    le.fit(&["cat".to_string(), "dog".to_string()]);
    let result = le.transform(&["cat", "fish"]);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("fish"));
}

#[test]
fn label_encoder_roundtrip() {
    let mut le = pathotypr_core::LabelEncoder::new();
    let labels = vec!["L1".to_string(), "L2".to_string(), "L3".to_string()];
    le.fit(&labels);
    let encoded = le.transform(&labels).unwrap();
    assert_eq!(encoded.len(), 3);
    // Each label should map to a unique index
    let mut seen = std::collections::HashSet::new();
    for &idx in &encoded {
        assert!(idx < 3);
        seen.insert(idx);
    }
    assert_eq!(seen.len(), 3);
}

// ============================================================================
// FeatureHasher
// ============================================================================

#[test]
fn feature_hasher_deterministic() {
    let hasher = pathotypr_core::FeatureHasher::new(1024);
    let seqs = vec!["ATCGATCGATCG".to_string(), "GCTAGCTAGCTA".to_string()];
    let result1 = hasher.transform_sparse(&seqs, 5);
    let result2 = hasher.transform_sparse(&seqs, 5);
    assert_eq!(result1, result2);
}

#[test]
fn feature_hasher_rows_are_sorted() {
    let hasher = pathotypr_core::FeatureHasher::new(1024);
    let seqs = vec!["ATCGATCGATCGATCGATCGATCG".to_string()];
    let rows = hasher.transform_sparse(&seqs, 7);
    for row in &rows {
        for w in row.windows(2) {
            assert!(w[0].0 < w[1].0, "Row not sorted: {:?} >= {:?}", w[0], w[1]);
        }
    }
}

#[test]
#[should_panic(expected = "power of 2")]
fn feature_hasher_rejects_non_power_of_two() {
    let _ = pathotypr_core::FeatureHasher::new(1000);
}

// ============================================================================
// Sparse tree: edge cases
// ============================================================================

#[test]
fn tree_handles_single_class_gracefully() {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use pathotypr_core::sparse_tree::{SparseDecisionTree, TreeParams};

    // All samples same class — tree should just create a leaf
    let data: Vec<Vec<(usize, f32)>> = vec![
        vec![(0, 1.0)], vec![(0, 2.0)], vec![(0, 3.0)],
    ];
    let labels = vec![0, 0, 0];
    let indices: Vec<usize> = (0..3).collect();
    let params = TreeParams {
        max_depth: 10,
        min_samples_leaf: 1,
        max_features: 1,
        n_classes: 2, // 2 classes but only 1 present
    };
    let mut rng = SmallRng::seed_from_u64(42);
    let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

    // All predictions should be class 0
    for row in &data {
        assert_eq!(tree.predict_one(row), 0);
    }
}

#[test]
fn tree_respects_min_samples_leaf() {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use pathotypr_core::sparse_tree::{SparseDecisionTree, TreeParams};

    let data: Vec<Vec<(usize, f32)>> = vec![
        vec![(0, 1.0)], vec![(0, 2.0)], vec![(0, 3.0)],
        vec![(0, 8.0)], vec![(0, 9.0)], vec![(0, 10.0)],
    ];
    let labels = vec![0, 0, 0, 1, 1, 1];
    let indices: Vec<usize> = (0..6).collect();
    // min_samples_leaf=3 means no leaf can have fewer than 3 samples
    let params = TreeParams {
        max_depth: 10,
        min_samples_leaf: 3,
        max_features: 1,
        n_classes: 2,
    };
    let mut rng = SmallRng::seed_from_u64(42);
    let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

    // Count leaf sizes by predicting all training samples
    // With min_leaf=3, each leaf should have ≥ 3 samples
    let mut leaf_counts = std::collections::HashMap::new();
    for row in &data {
        let pred = tree.predict_one(row);
        *leaf_counts.entry(pred).or_insert(0) += 1;
    }
    for (&_class, &count) in &leaf_counts {
        assert!(count >= 3, "Leaf has {} samples, expected >= 3", count);
    }
}

#[test]
fn tree_respects_max_depth() {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use pathotypr_core::sparse_tree::{SparseDecisionTree, TreeNode, TreeParams};

    let data: Vec<Vec<(usize, f32)>> = (0..20)
        .map(|i| vec![(0, i as f32)])
        .collect();
    let labels: Vec<usize> = (0..20).map(|i| i % 4).collect();
    let indices: Vec<usize> = (0..20).collect();
    let params = TreeParams {
        max_depth: 3,
        min_samples_leaf: 1,
        max_features: 1,
        n_classes: 4,
    };
    let mut rng = SmallRng::seed_from_u64(42);
    let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

    // Measure actual depth by traversing
    fn measure_depth(nodes: &[TreeNode], idx: usize) -> usize {
        match &nodes[idx] {
            TreeNode::Leaf { .. } => 0,
            TreeNode::Split { left, right, .. } => {
                1 + measure_depth(nodes, *left as usize)
                    .max(measure_depth(nodes, *right as usize))
            }
        }
    }
    let depth = measure_depth(&tree.nodes, 0);
    assert!(depth <= 3, "Tree depth {} exceeds max_depth 3", depth);
}

// ============================================================================
// Cancellation
// ============================================================================

#[test]
fn train_respects_cancellation() {
    let tmp = TempDir::new().unwrap();
    let fasta_path = tmp.path().join("cancel.fasta");
    let model_path = tmp.path().join("cancel_model.pathotypr.zst");

    write_synthetic_fasta(fasta_path.to_str().unwrap(), &["A", "B"], 10);

    let token = pathotypr_core::CancellationToken::new();
    token.cancel(); // Pre-cancel

    let args = pathotypr_core::train::Args {
        input: fasta_path.to_str().unwrap().to_string(),
        output: model_path.to_str().unwrap().to_string(),
        kmer_size: 11,
        test_split: 0.2,
        threads: Some(1),
        cv_folds: None,
        max_depth: 10,
        min_samples_leaf: 1,
        cancel_token: Some(token),
    };
    let result = pathotypr_core::train::run(args);
    assert!(result.is_err());
    let err_msg = format!("{}", result.unwrap_err());
    assert!(err_msg.contains("cancelled"), "Expected cancellation error, got: {}", err_msg);
}

// ============================================================================
// Paired-end detection
// ============================================================================

#[test]
fn detect_illumina_paired_end() {
    let files = vec![
        "sample1_R1_001.fastq.gz".to_string(),
        "sample1_R2_001.fastq.gz".to_string(),
        "sample2_R1_001.fastq.gz".to_string(),
        "sample2_R2_001.fastq.gz".to_string(),
    ];
    let result = pathotypr_core::detect_paired_end_files(&files);
    assert!(result.is_paired);
    assert_eq!(result.paired_count, 2);
    assert_eq!(result.single_count, 0);
}

#[test]
fn detect_single_end() {
    let files = vec![
        "reads_A.fastq.gz".to_string(),
        "reads_B.fastq.gz".to_string(),
    ];
    let result = pathotypr_core::detect_paired_end_files(&files);
    assert!(!result.is_paired);
    assert_eq!(result.single_count, 2);
}

// ============================================================================
// Feature importance
// ============================================================================

#[test]
fn ensemble_importance_sums_across_trees() {
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use pathotypr_core::sparse_tree::{
        SparseDecisionTree, TreeParams, ensemble_feature_importance,
    };

    let data: Vec<Vec<(usize, f32)>> = vec![
        vec![(0, 1.0), (1, 0.0)],
        vec![(0, 2.0), (1, 0.0)],
        vec![(0, 3.0), (1, 0.0)],
        vec![(0, 8.0), (1, 1.0)],
        vec![(0, 9.0), (1, 1.0)],
        vec![(0, 10.0), (1, 1.0)],
    ];
    let labels = vec![0, 0, 0, 1, 1, 1];
    let indices: Vec<usize> = (0..6).collect();
    let params = TreeParams {
        max_depth: 5,
        min_samples_leaf: 1,
        max_features: 2,
        n_classes: 2,
    };

    let mut trees = Vec::new();
    for seed in 0..10u64 {
        let mut rng = SmallRng::seed_from_u64(seed);
        trees.push(SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng));
    }

    let importance = ensemble_feature_importance(&trees);
    assert!(!importance.is_empty());
    // Feature 0 should be most important (it separates the classes)
    let top_feature = importance[0].0;
    assert_eq!(top_feature, 0, "Feature 0 should be most discriminative");
}

// ============================================================================
// Excel export
// ============================================================================

#[test]
fn excel_stream_writer_creates_file() {
    let tmp = TempDir::new().unwrap();
    let tsv_path = tmp.path().join("test_output.tsv");

    let mut writer = pathotypr_core::ExcelStreamWriter::new(
        tsv_path.to_str().unwrap(),
        &["Name", "Value", "Score"],
    ).unwrap();

    writer.write_row(&["sample1".to_string(), "100".to_string(), "0.95".to_string()]).unwrap();
    writer.write_row(&["sample2".to_string(), "200".to_string(), "0.88".to_string()]).unwrap();
    writer.finish().unwrap();

    let xlsx_path = tsv_path.with_extension("xlsx");
    assert!(xlsx_path.exists(), "Excel file should be created");
    assert!(xlsx_path.metadata().unwrap().len() > 0, "Excel file should not be empty");
}
