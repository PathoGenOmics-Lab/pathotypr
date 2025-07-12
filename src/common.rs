//! src/common.rs
//!
//! This module provides shared, high-performance data structures and utility
//! functions for the `pathotypr` application. It has been optimized to reduce
//! memory allocations and improve processing speed, especially for the `train`
//! and `predict` subcommands.

// --- Crates for performance and serialization ---
use rayon::prelude::*;
use rustc_hash::FxHashMap; // A much faster hasher for integer keys
use serde::{Deserialize, Serialize};
use needletail::Sequence; // Import the Sequence trait to use its methods

// --- SmartCore components for machine learning ---
// MODIFIED: All types now use f32 for memory efficiency.
use smartcore::linalg::basic::matrix::DenseMatrix;
use smartcore::tree::decision_tree_classifier::DecisionTreeClassifier;

// --- Model Bundle Structs ---

/// Defines the configuration of a trained model.
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelConfig {
    pub pathotypr_version: String,
    pub kmer_size: usize,
    pub n_trees: u16,
}

/// A unified, compressed bundle containing everything needed for prediction.
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelBundle {
    pub config: ModelConfig,
    pub vectorizer: CountVectorizer,
    pub label_encoder: LabelEncoder,
    // MODIFIED: The DecisionTreeClassifier now uses f32.
    pub trees: Vec<DecisionTreeClassifier<f32, usize, DenseMatrix<f32>, Vec<usize>>>,
}

// --- Feature Processing Components ---

/// Transforms sequences into k-mer count vectors using `u64` representation.
#[derive(Serialize, Deserialize, Debug)]
pub struct CountVectorizer {
    pub vocabulary: FxHashMap<u64, usize>,
    pub num_features: usize,
}

impl CountVectorizer {
    /// Creates a new, empty `CountVectorizer`.
    pub fn new() -> Self {
        Self {
            vocabulary: FxHashMap::default(),
            num_features: 0,
        }
    }

    /// Builds the vocabulary from a pre-computed map of k-mer counts.
    pub fn fit(&mut self, kmer_counts: &FxHashMap<u64, u32>) {
        let mut vocab_idx = 0;
        for &kmer_hash in kmer_counts.keys() {
            self.vocabulary.insert(kmer_hash, vocab_idx);
            vocab_idx += 1;
        }
        self.num_features = self.vocabulary.len();
    }

    /// MODIFIED: Transforms sequences into a sparse data format using f32.
    /// Instead of a giant Vec<Vec<f64>>, this returns a Vec of sparse vectors.
    /// Each sparse vector is a Vec of (feature_index, count).
    pub fn transform_sparse(&self, sequences: &[String], k: usize) -> Vec<Vec<(usize, f32)>> {
        sequences
            .par_iter()
            .map(|seq| {
                // Use a temporary HashMap to count k-mers for this sequence only.
                // This is memory-efficient as it's local to the sequence.
                let mut sequence_kmer_counts: FxHashMap<usize, f32> = FxHashMap::default();
                for (_, bitkmer_tuple, _) in seq.as_bytes().bit_kmers(k as u8, true) {
                    let kmer_hash = bitkmer_tuple.0;
                    if let Some(&idx) = self.vocabulary.get(&kmer_hash) {
                        *sequence_kmer_counts.entry(idx).or_insert(0.0) += 1.0;
                    }
                }
                // Convert the map to a vector of (index, value) tuples.
                // Sorting is good practice for some sparse matrix formats.
                let mut features: Vec<(usize, f32)> = sequence_kmer_counts.into_iter().collect();
                features.sort_unstable_by_key(|&(idx, _)| idx);
                features
            })
            .collect()
    }
}

/// Encodes string labels into integer representations and vice-versa.
#[derive(Serialize, Deserialize, Debug)]
pub struct LabelEncoder {
    pub label_to_int: FxHashMap<String, usize>,
    pub int_to_label: Vec<String>,
}

impl LabelEncoder {
    /// Creates a new, empty `LabelEncoder`.
    pub fn new() -> Self {
        Self {
            label_to_int: FxHashMap::default(),
            int_to_label: Vec::new(),
        }
    }

    /// Learns the mapping from a slice of string labels.
    pub fn fit<T: AsRef<str> + std::hash::Hash + std::cmp::Eq>(&mut self, labels: &[T]) {
        for label in labels {
            let label_str = label.as_ref().to_string();
            self.label_to_int.entry(label_str.clone()).or_insert_with(|| {
                let index = self.int_to_label.len();
                self.int_to_label.push(label_str);
                index
            });
        }
    }

    /// Transforms a slice of string labels into their integer representations.
    pub fn transform<T: AsRef<str>>(&self, labels: &[T]) -> Vec<usize> {
        labels
            .iter()
            .map(|label| {
                *self.label_to_int.get(label.as_ref())
                    .expect("Label not found in encoder. `fit` must be called first with all possible labels.")
            })
            .collect()
    }
}