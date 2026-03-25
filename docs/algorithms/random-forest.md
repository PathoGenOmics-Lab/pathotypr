# Random Forest on Sparse Vectors

## Overview

pathotypr implements a custom Random Forest classifier that operates **directly on sparse feature vectors** — no dense matrix conversion required. This eliminates the O(N × 1M) memory overhead that a library like smartcore would require for N samples with 1M feature buckets.

## Decision Tree (CART)

Each tree is a binary classifier built with the CART algorithm (Classification and Regression Trees):

### Node Types

```rust
enum TreeNode {
    Split { feature: u32, threshold: f32, left: u32, right: u32 },
    Leaf { class: u32 },
}
```

Nodes are stored in a flat arena (`Vec<TreeNode>`) with index-based child pointers. This gives cache-friendly traversal and compact serialization.

### Tree Construction

**Input**: Sparse feature vectors `[(bucket_idx, count)]`, class labels, bootstrap sample indices.

**Algorithm** (iterative DFS via explicit stack):

1. Start with all bootstrap samples at the root
2. For each node:
   - Count class frequencies in the current sample set
   - **Stopping conditions**: pure node (single class), `depth >= max_depth`, or `n_samples <= 2 × min_samples_leaf`
   - If not stopped, find the best split (see below)
   - Partition samples by the split and push children onto the stack

### Finding the Best Split

1. **Feature subsampling**: Collect all features present in the current node's data. If more than `sqrt(n_features)` features are present, randomly sample `sqrt(n_features)` via partial Fisher-Yates shuffle.

2. **For each candidate feature**:
   - Extract `(value, label)` pairs via binary search on each sample's sparse row
   - Sort by value
   - Scan sorted values to find the threshold that maximizes Gini impurity reduction
   - Only consider splits where both sides have ≥ `min_samples_leaf` samples

3. **Gini impurity**: `1 - Σ(p_i²)` where p_i is the proportion of class i.

4. **Threshold**: Midpoint between consecutive distinct values at the best split point.

### Prediction

Traversal from root to leaf in O(tree_depth):

```rust
fn predict_one(row: &[(usize, f32)]) -> usize {
    loop {
        match node {
            Leaf { class } => return class,
            Split { feature, threshold, left, right } => {
                let val = binary_search(row, feature);  // O(log nnz)
                node = if val <= threshold { left } else { right };
            }
        }
    }
}
```

The `sparse_lookup()` function uses binary search on the sorted sparse row, returning 0.0 for absent features.

## Ensemble (Bagging)

### Training

- **100 trees** trained in parallel with rayon
- Each tree gets a **bootstrap sample**: N samples drawn with replacement from N total (each sample appears in ~63% of trees)
- **Deterministic seeds**: A master RNG (seed=42) generates one seed per tree. Each tree's bootstrap is reproducible from its seed.

### Majority Voting

For prediction, all 100 trees vote. The class with the most votes wins:

```
confidence = winner_votes / total_trees
margin = (winner_votes - runner_up_votes) / total_trees
```

Both metrics are reported alongside the predicted class.

## Out-of-Bag (OOB) Accuracy

Each tree's bootstrap sample excludes ~37% of training data. These "out-of-bag" samples provide a free accuracy estimate:

1. For each tree, regenerate its bootstrap using the stored seed
2. Identify OOB samples (not in the bootstrap)
3. Predict OOB samples using only that tree
4. Aggregate: for each training sample, collect votes from all trees where it was OOB
5. Compare majority vote to true label

**Why this matters**: OOB accuracy is nearly unbiased and requires no held-out test set. It's always computed, even when using CV for accuracy estimation.

## Feature Importance

Split-count importance: for each feature (bucket), count how many times it was used as a split across all 100 trees. Features that appear in many splits are more discriminative.

```
importance(feature_j) = Σ_{trees} count_of_splits_on_feature_j
```

The top 500 features are exported with:
- Their k-mer sequences (reverse-mapped from the hashing trick)
- Genomic coordinates (sequence index + position) in the training data

## Hyperparameters

| Parameter | Default | Effect |
|---|---|---|
| `n_trees` | 100 | More trees → smoother predictions, diminishing returns past ~100 |
| `max_depth` | 20 | Limits tree complexity; prevents memorization of training noise |
| `min_samples_leaf` | 5 | Regularization — prevents splits that isolate individual samples |
| `max_features` | √(n_features) | Standard RF: decorrelates trees by forcing different feature subsets |
| `bootstrap_seed` | 42 (master) | Ensures reproducibility; each tree derives its seed from this |

## Complexity

| Operation | Time | Space |
|---|---|---|
| Train one tree | O(n_bootstrap × sqrt(F) × depth) | O(nodes) ≈ O(2^depth) |
| Train ensemble | Above / num_threads (rayon) | O(100 × nodes) |
| Predict one sample | O(100 × depth × log(nnz)) | O(1) |
| OOB accuracy | O(N × 100 × depth) | O(N × n_classes) |
