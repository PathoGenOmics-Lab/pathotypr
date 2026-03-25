//! Sparse Decision Tree and Random Forest implementation.
//!
//! Custom CART decision tree that operates directly on sparse feature vectors
//! `Vec<(usize, f32)>`, eliminating the need for dense matrix conversion.
//! Key advantages over smartcore's DecisionTreeClassifier:
//! - No dense buffer allocation (O(nnz) vs O(n_features) per sample)
//! - Feature subsampling: evaluates sqrt(n_features) per split (standard RF)
//! - Prediction in O(tree_depth) per sample via binary search on sparse rows

use rand::rngs::SmallRng;
use rand::Rng;
use serde::{Deserialize, Serialize};

/// A single node in the decision tree, stored in an arena (Vec<TreeNode>).
#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum TreeNode {
    /// Internal split node.
    Split {
        feature: u32,
        threshold: f32,
        left: u32,
        right: u32,
    },
    /// Leaf node with class prediction.
    Leaf {
        class: u32,
    },
}

/// A decision tree that operates on sparse feature vectors.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SparseDecisionTree {
    pub nodes: Vec<TreeNode>,
    pub n_classes: u32,
}

/// Hyperparameters for tree construction.
pub struct TreeParams {
    pub max_depth: usize,
    pub min_samples_leaf: usize,
    pub max_features: usize,
    pub n_classes: usize,
}

/// Look up a feature value in a sorted sparse row via binary search.
/// Returns 0.0 if the feature is not present.
#[inline]
fn sparse_lookup(row: &[(usize, f32)], feature: usize) -> f32 {
    match row.binary_search_by_key(&feature, |&(f, _)| f) {
        Ok(idx) => row[idx].1,
        Err(_) => 0.0,
    }
}

/// Compute Gini impurity for a class distribution.
#[inline]
fn gini_impurity(class_counts: &[u32], total: u32) -> f32 {
    if total == 0 {
        return 0.0;
    }
    let inv_total = 1.0 / total as f32;
    let mut sum_sq = 0.0_f32;
    for &count in class_counts {
        let p = count as f32 * inv_total;
        sum_sq += p * p;
    }
    1.0 - sum_sq
}

impl SparseDecisionTree {
    /// Train a decision tree on sparse data with bootstrap sample indices.
    ///
    /// `data` contains the sparse rows (sorted by feature index).
    /// `labels` contains the class label for each sample.
    /// `sample_indices` are the indices into data/labels for this bootstrap sample.
    pub fn fit(
        data: &[Vec<(usize, f32)>],
        labels: &[usize],
        sample_indices: &[usize],
        params: &TreeParams,
        rng: &mut SmallRng,
    ) -> Self {
        let mut nodes = Vec::with_capacity(256);

        // Root: placeholder, will be filled
        nodes.push(TreeNode::Leaf { class: 0 });

        // DFS build via explicit stack: (node_idx, sample_indices, depth)
        let mut queue: Vec<(usize, Vec<usize>, usize)> = Vec::new();
        queue.push((0, sample_indices.to_vec(), 0));

        while let Some((node_idx, indices, depth)) = queue.pop() {
            let n = indices.len();

            // Count classes
            let mut class_counts = vec![0u32; params.n_classes];
            for &idx in &indices {
                class_counts[labels[idx]] += 1;
            }

            // Find majority class (class_counts always has n_classes ≥ 2 elements)
            let (majority_class, &majority_count) = class_counts
                .iter()
                .enumerate()
                .max_by_key(|(_, &c)| c)
                .unwrap_or((0, &0));

            // Leaf conditions: pure node, max depth, min samples
            if majority_count == n as u32
                || depth >= params.max_depth
                || n <= params.min_samples_leaf * 2
            {
                nodes[node_idx] = TreeNode::Leaf {
                    class: majority_class as u32,
                };
                continue;
            }

            // Find best split among random subset of features
            let best = find_best_split(
                data,
                labels,
                &indices,
                &class_counts,
                params,
                rng,
            );

            if let Some((best_feature, best_threshold, left_indices, right_indices)) = best {
                // Allocate child nodes
                let left_idx = nodes.len() as u32;
                nodes.push(TreeNode::Leaf { class: 0 }); // placeholder
                let right_idx = nodes.len() as u32;
                nodes.push(TreeNode::Leaf { class: 0 }); // placeholder

                nodes[node_idx] = TreeNode::Split {
                    feature: best_feature as u32,
                    threshold: best_threshold,
                    left: left_idx,
                    right: right_idx,
                };

                queue.push((left_idx as usize, left_indices, depth + 1));
                queue.push((right_idx as usize, right_indices, depth + 1));
            } else {
                // No valid split found
                nodes[node_idx] = TreeNode::Leaf {
                    class: majority_class as u32,
                };
            }
        }

        SparseDecisionTree {
            nodes,
            n_classes: params.n_classes as u32,
        }
    }

    /// Predict class for a single sparse row.
    #[inline]
    pub fn predict_one(&self, row: &[(usize, f32)]) -> usize {
        let mut node_idx = 0;
        loop {
            match &self.nodes[node_idx] {
                TreeNode::Leaf { class } => return *class as usize,
                TreeNode::Split {
                    feature,
                    threshold,
                    left,
                    right,
                } => {
                    let val = sparse_lookup(row, *feature as usize);
                    node_idx = if val <= *threshold {
                        *left as usize
                    } else {
                        *right as usize
                    };
                }
            }
        }
    }

    /// Predict classes for a batch of sparse rows.
    pub fn predict_batch(&self, rows: &[&[(usize, f32)]]) -> Vec<usize> {
        rows.iter().map(|row| self.predict_one(row)).collect()
    }

    /// Count how many times each feature (bucket) is used as a split in this tree.
    /// Returns a map of feature_index -> split_count.
    pub fn split_feature_counts(&self) -> Vec<(u32, u32)> {
        let mut counts: rustc_hash::FxHashMap<u32, u32> = rustc_hash::FxHashMap::default();
        for node in &self.nodes {
            if let TreeNode::Split { feature, .. } = node {
                *counts.entry(*feature).or_insert(0) += 1;
            }
        }
        counts.into_iter().collect()
    }
}

/// Compute aggregate feature importance across an ensemble of trees.
/// Returns a sorted Vec of (feature_index, total_split_count) in descending order.
pub fn ensemble_feature_importance(trees: &[SparseDecisionTree]) -> Vec<(u32, u32)> {
    let mut total: rustc_hash::FxHashMap<u32, u32> = rustc_hash::FxHashMap::default();
    for tree in trees {
        for (feat, count) in tree.split_feature_counts() {
            *total.entry(feat).or_insert(0) += count;
        }
    }
    let mut sorted: Vec<(u32, u32)> = total.into_iter().collect();
    sorted.sort_unstable_by(|a, b| b.1.cmp(&a.1));
    sorted
}

/// Find the best split among a random subset of features.
fn find_best_split(
    data: &[Vec<(usize, f32)>],
    labels: &[usize],
    indices: &[usize],
    parent_class_counts: &[u32],
    params: &TreeParams,
    rng: &mut SmallRng,
) -> Option<(usize, f32, Vec<usize>, Vec<usize>)> {
    let n = indices.len();
    let n_classes = params.n_classes;
    let parent_total = n as u32;
    let parent_gini = gini_impurity(parent_class_counts, parent_total);

    let mut best_gain = 0.0_f32;
    let mut best_feature = 0usize;
    let mut best_threshold = 0.0_f32;

    // Sample max_features random features to evaluate
    // Use reservoir sampling over the feature space
    let n_features_to_try = params.max_features;

    // Collect unique features present in this node's data for smarter sampling
    let mut present_features: Vec<usize> = Vec::new();
    {
        use rustc_hash::FxHashSet;
        let mut seen = FxHashSet::default();
        for &idx in indices {
            for &(f, _) in &data[idx] {
                if seen.insert(f) {
                    present_features.push(f);
                }
            }
        }
    }

    // If very few features present, evaluate all of them; otherwise sample
    let features_to_eval: Vec<usize> = if present_features.len() <= n_features_to_try {
        present_features
    } else {
        // Fisher-Yates partial shuffle
        let mut feats = present_features;
        let k = n_features_to_try.min(feats.len());
        for i in 0..k {
            let j = rng.gen_range(i..feats.len());
            feats.swap(i, j);
        }
        feats.truncate(k);
        feats
    };

    // Temporary buffer for feature values
    let mut feature_vals: Vec<(f32, usize)> = Vec::with_capacity(n);

    for &feature in &features_to_eval {
        // Collect (value, label) for this feature across node samples
        feature_vals.clear();
        for &idx in indices {
            let val = sparse_lookup(&data[idx], feature);
            feature_vals.push((val, labels[idx]));
        }

        // Sort by feature value
        feature_vals.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Scan sorted values to find best threshold
        let mut left_counts = vec![0u32; n_classes];
        let mut right_counts: Vec<u32> = parent_class_counts.to_vec();
        let mut left_total = 0u32;
        let mut right_total = parent_total;

        for i in 0..n - 1 {
            let (val, label) = feature_vals[i];
            left_counts[label] += 1;
            right_counts[label] -= 1;
            left_total += 1;
            right_total -= 1;

            // Only consider split if next value is different
            let next_val = feature_vals[i + 1].0;
            if (next_val - val).abs() < f32::EPSILON {
                continue;
            }

            // Skip if either side would be too small
            if (left_total as usize) < params.min_samples_leaf
                || (right_total as usize) < params.min_samples_leaf
            {
                continue;
            }

            let left_gini = gini_impurity(&left_counts, left_total);
            let right_gini = gini_impurity(&right_counts, right_total);

            let weighted_gini = (left_total as f32 * left_gini
                + right_total as f32 * right_gini)
                / parent_total as f32;
            let gain = parent_gini - weighted_gini;

            if gain > best_gain {
                best_gain = gain;
                best_feature = feature;
                best_threshold = (val + next_val) * 0.5;
            }
        }
    }

    if best_gain <= 0.0 {
        return None;
    }

    // Partition indices by best split
    let mut left_indices = Vec::new();
    let mut right_indices = Vec::new();
    for &idx in indices {
        let val = sparse_lookup(&data[idx], best_feature);
        if val <= best_threshold {
            left_indices.push(idx);
        } else {
            right_indices.push(idx);
        }
    }

    if left_indices.is_empty() || right_indices.is_empty() {
        return None;
    }

    Some((best_feature, best_threshold, left_indices, right_indices))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn make_params(n_features: usize, n_classes: usize) -> TreeParams {
        TreeParams {
            max_depth: 10,
            min_samples_leaf: 1,
            max_features: (n_features as f64).sqrt().ceil() as usize,
            n_classes,
        }
    }

    #[test]
    fn test_pure_split() {
        // Two linearly separable classes based on feature 0
        let data: Vec<Vec<(usize, f32)>> = vec![
            vec![(0, 1.0)],
            vec![(0, 2.0)],
            vec![(0, 3.0)],
            vec![(0, 8.0)],
            vec![(0, 9.0)],
            vec![(0, 10.0)],
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let indices: Vec<usize> = (0..6).collect();
        let params = make_params(1, 2);
        let mut rng = SmallRng::seed_from_u64(42);

        let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

        // Should perfectly classify
        for (i, row) in data.iter().enumerate() {
            assert_eq!(tree.predict_one(row), labels[i], "Failed on sample {}", i);
        }
    }

    #[test]
    fn test_sparse_lookup_missing() {
        let row = vec![(2, 3.0), (5, 7.0), (10, 1.0)];
        assert_eq!(sparse_lookup(&row, 0), 0.0);
        assert_eq!(sparse_lookup(&row, 2), 3.0);
        assert_eq!(sparse_lookup(&row, 5), 7.0);
        assert_eq!(sparse_lookup(&row, 7), 0.0);
    }

    #[test]
    fn test_multiclass() {
        // 3 classes, each with a distinct feature pattern
        let data: Vec<Vec<(usize, f32)>> = vec![
            vec![(0, 10.0)],
            vec![(0, 11.0)],
            vec![(1, 10.0)],
            vec![(1, 11.0)],
            vec![(2, 10.0)],
            vec![(2, 11.0)],
        ];
        let labels = vec![0, 0, 1, 1, 2, 2];
        let indices: Vec<usize> = (0..6).collect();
        let params = make_params(3, 3);
        let mut rng = SmallRng::seed_from_u64(42);

        let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

        for (i, row) in data.iter().enumerate() {
            assert_eq!(tree.predict_one(row), labels[i], "Failed on sample {}", i);
        }
    }

    #[test]
    fn test_serialization_roundtrip() {
        // Use enough samples per class for reliable splitting
        let data: Vec<Vec<(usize, f32)>> = vec![
            vec![(0, 1.0)], vec![(0, 2.0)], vec![(0, 3.0)],
            vec![(0, 8.0)], vec![(0, 9.0)], vec![(0, 10.0)],
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let indices: Vec<usize> = (0..6).collect();
        let params = make_params(1, 2);
        let mut rng = SmallRng::seed_from_u64(42);

        let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

        // Serialize and deserialize
        let bytes = bincode::serialize(&tree).unwrap();
        let tree2: SparseDecisionTree = bincode::deserialize(&bytes).unwrap();

        for (i, row) in data.iter().enumerate() {
            assert_eq!(tree2.predict_one(row), labels[i], "Failed on sample {} after roundtrip", i);
        }
    }

    #[test]
    fn test_predict_batch() {
        let data: Vec<Vec<(usize, f32)>> = vec![
            vec![(0, 1.0)], vec![(0, 2.0)], vec![(0, 3.0)],
            vec![(0, 8.0)], vec![(0, 9.0)], vec![(0, 10.0)],
        ];
        let labels = vec![0, 0, 0, 1, 1, 1];
        let indices: Vec<usize> = (0..6).collect();
        let params = make_params(1, 2);
        let mut rng = SmallRng::seed_from_u64(42);

        let tree = SparseDecisionTree::fit(&data, &labels, &indices, &params, &mut rng);

        let rows: Vec<&[(usize, f32)]> = data.iter().map(|r| r.as_slice()).collect();
        let preds = tree.predict_batch(&rows);
        assert_eq!(preds, labels);
    }
}
