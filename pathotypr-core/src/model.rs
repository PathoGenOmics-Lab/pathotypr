//! Model bundle definitions and label encoding.
//!
//! Contains the serializable structures for trained models:
//! `ModelBundle`, `ModelConfig`, and `LabelEncoder`.

use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

use crate::sparse_tree::SparseDecisionTree;
use crate::vectorizer::FeatureHasher;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Current model format version. Increment when breaking serialization compatibility.
pub const MODEL_FORMAT_VERSION: u32 = 3;

/// Default number of hash buckets: 2^20 = 1,048,576.
pub const DEFAULT_HASH_BUCKETS: usize = 1 << 20;

// ---------------------------------------------------------------------------
// Model bundle
// ---------------------------------------------------------------------------

/// Configuration metadata stored alongside a trained model.
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelConfig {
    pub pathotypr_version: String,
    pub kmer_size: usize,
    pub n_trees: u16,
    /// Model format version for compatibility checks.
    #[serde(default)]
    pub format_version: u32,
}

/// A unified, compressed bundle containing everything needed for prediction.
#[derive(Serialize, Deserialize, Debug)]
pub struct ModelBundle {
    pub config: ModelConfig,
    pub vectorizer: FeatureHasher,
    pub label_encoder: LabelEncoder,
    pub trees: Vec<SparseDecisionTree>,
}

// ---------------------------------------------------------------------------
// Label encoder
// ---------------------------------------------------------------------------

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
            self.label_to_int
                .entry(label_str.clone())
                .or_insert_with(|| {
                    let index = self.int_to_label.len();
                    self.int_to_label.push(label_str);
                    index
                });
        }
    }

    /// Transforms a slice of string labels into their integer representations.
    ///
    /// Returns `Err` if any label was not seen during `fit`.
    pub fn transform<T: AsRef<str>>(&self, labels: &[T]) -> Result<Vec<usize>, String> {
        labels
            .iter()
            .map(|label| {
                self.label_to_int
                    .get(label.as_ref())
                    .copied()
                    .ok_or_else(|| format!("Unknown label '{}' not found in encoder. Call `fit` first with all labels.", label.as_ref()))
            })
            .collect()
    }
}
