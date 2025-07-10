//! src/common.rs
//!
//! Module for shared data structures and utility functions across the application.
//!
//! This module centralizes components like the model bundle, vectorizer,
//! and label encoder to avoid code duplication between the `train` and `predict`
//! modules.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use smartcore::ensemble::random_forest_classifier::RandomForestClassifier;
use smartcore::linalg::basic::matrix::DenseMatrix;
use std::collections::HashMap;

// --- Model Bundle Structs ---

/// Configuration for the model, such as the k-mer size.
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelConfig {
    pub kmer_size: usize,
}

/// A unified bundle containing everything needed for prediction.
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelBundle {
    pub config: ModelConfig,
    pub vectorizer: CountVectorizer,
    pub label_encoder: LabelEncoder,
    pub model: RandomForestClassifier<f64, usize, DenseMatrix<f64>, Vec<usize>>,
}

// --- Feature Processing Components ---

/// Transforms text into k-mer count vectors.
#[derive(Serialize, Deserialize, Debug)]
pub struct CountVectorizer {
    pub vocabulary: HashMap<String, usize>,
    pub feature_names: Vec<String>,
}

impl CountVectorizer {
    pub fn new() -> Self {
        Self {
            vocabulary: HashMap::new(),
            feature_names: Vec::new(),
        }
    }

    /// Builds the vocabulary from a collection of texts.
    pub fn fit<T: AsRef<str>>(&mut self, texts: &[T]) {
        let mut freq: HashMap<String, usize> = HashMap::new();
        for text in texts {
            for token in text.as_ref().split_whitespace() {
                *freq.entry(token.to_string()).or_insert(0) += 1;
            }
        }
        let mut freq_vec: Vec<(String, usize)> = freq.into_iter().collect();
        freq_vec.sort_by(|a, b| b.1.cmp(&a.1));
        self.vocabulary = freq_vec
            .iter()
            .enumerate()
            .map(|(i, (token, _))| (token.clone(), i))
            .collect();
        self.feature_names = freq_vec.into_iter().map(|(token, _)| token).collect();
    }

    /// Transforms a collection of texts into a feature matrix.
    pub fn transform<T: AsRef<str> + Sync>(&self, texts: &[T]) -> Vec<Vec<f64>> {
        texts
            .par_iter()
            .map(|text| {
                let mut counts = vec![0.0; self.vocabulary.len()];
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

/// Encodes class labels (strings) into integers.
#[derive(Serialize, Deserialize, Debug)]
pub struct LabelEncoder {
    pub label_to_int: HashMap<String, usize>,
    pub int_to_label: Vec<String>,
}

impl LabelEncoder {
    pub fn new() -> Self {
        Self {
            label_to_int: HashMap::new(),
            int_to_label: Vec::new(),
        }
    }

    /// Learns the label mapping from a collection of labels.
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

    /// Transforms a collection of labels into their integer representations.
    pub fn transform<T: AsRef<str>>(&self, labels: &[T]) -> Vec<usize> {
        labels
            .iter()
            .map(|label| *self.label_to_int.get(label.as_ref()).unwrap())
            .collect()
    }
}

// --- Utility Functions ---

/// Generates a space-separated string of k-mers from a DNA sequence.
///
/// # Arguments
/// * `sequence` - The input DNA sequence.
/// * `k` - The k-mer size.
///
/// # Returns
/// A `String` of k-mers. Returns an empty string if the sequence is shorter than `k`.
pub fn kmerize(sequence: &str, k: usize) -> String {
    if sequence.len() < k {
        return String::new();
    }
    (0..=sequence.len() - k)
        .map(|i| &sequence[i..i + k])
        .collect::<Vec<&str>>()
        .join(" ")
}