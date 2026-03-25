# Training Pipeline

## End-to-End Flow

```
Labeled FASTA ──→ Feature Hashing ──→ Accuracy Estimation ──→ Train Final Model ──→ OOB ──→ Save
                  (k-mer → sparse)    (CV or single split)   (all data, 100 trees)         (bincode+zstd)
                                                                      │
                                                                      ↓
                                                              Feature Importance
                                                              (reverse map k-mers + coords)
```

## Step 1: Data Preparation

### Input Format

Labeled multi-FASTA where the **first token** in each header is the class label:

```
>L4 sample_0001
ACTGACTG...
>L2 sample_0002
ACTGACTG...
```

### Vectorization

Each sequence is transformed into a sparse feature vector via [feature hashing](feature-hashing.md):

1. Extract all k-mers (default k=21) using needletail's 2-bit encoding
2. Hash each k-mer into one of 2^20 buckets via bitmask
3. Count occurrences per bucket
4. Output: sorted sparse vector `Vec<(bucket_idx, count)>`

### Label Encoding

Class labels (strings) are mapped to integers via `LabelEncoder`:
- `fit()`: Builds bidirectional mapping (label ↔ integer)
- `transform()`: Converts label strings to integer indices
- Minimum 2 distinct classes required

## Step 2: Accuracy Estimation

Two modes, controlled by `--cv-folds`:

### Single Train/Test Split (default)

- Shuffle indices with seed=42
- Hold out `test_split` fraction (default 20%)
- Train a temporary 100-tree ensemble on the training portion
- Evaluate on the held-out test set
- Report accuracy percentage

### Stratified k-Fold Cross-Validation (`--cv-folds k`)

- Group sample indices by class
- Shuffle within each class (seed=42)
- Assign to folds round-robin → preserves class proportions
- For each fold:
  - Train 100-tree ensemble on k-1 folds
  - Evaluate on the held-out fold
- Report: mean accuracy ± standard deviation (Bessel's correction: divide by n-1)

**Key**: The accuracy estimation step trains **throwaway** ensembles. The final model is always trained on **all** data.

## Step 3: Final Model Training

- Uses **all** samples (no held-out set)
- Trains 100 trees in parallel via rayon (see [random-forest.md](random-forest.md))
- Returns both trees and their bootstrap seeds (for OOB computation)

## Step 4: Out-of-Bag Accuracy

Always computed on the final model (see [random-forest.md](random-forest.md#out-of-bag-oob-accuracy)):

- Regenerates each tree's bootstrap from its seed
- For each sample, collects votes only from trees where it was OOB (~37% of trees)
- Reports majority-vote accuracy across all samples

OOB provides a nearly unbiased accuracy estimate at zero additional cost.

## Step 5: Model Serialization

### Format

```
ModelBundle {
    config: ModelConfig {
        pathotypr_version,   // e.g., "0.2.0"
        kmer_size,           // e.g., 21
        n_trees,             // 100
        format_version,      // 3 (current)
    },
    vectorizer: FeatureHasher { num_buckets },
    label_encoder: LabelEncoder { maps },
    trees: Vec<SparseDecisionTree>,
}
```

### Compression Pipeline

```
ModelBundle ──→ bincode::serialize_into ──→ zstd::Encoder (level 3) ──→ file
```

Streaming serialization: bincode writes directly into the zstd compressor, which writes to a buffered file. Never holds the full uncompressed + compressed model in memory simultaneously.

**Typical sizes**: 5–50 MB compressed for real bacterial genomes; <3 KB for synthetic benchmarks.

## Step 6: Feature Importance Export

### Importance TSV (`model.importance.tsv`)

Top 500 features ranked by split count:

```
rank  bucket  split_count  importance_pct  kmers
1     42      87           2.31            ACTGACTGCTAGCTGATCGATC,GATCGATCGATCGATCGATCG
2     1087    73           1.94            ...
```

`importance_pct = split_count / total_splits_across_all_features × 100`

### Genomic Coordinates TSV (`model.importance.coords.tsv`)

Maps each discriminant k-mer back to its physical location in the training sequences:

```
rank  bucket  split_count  importance_pct  kmer  sequence  lineage  position
1     42      87           2.31            ACTG...  seq_header  L4  1234567
```

This enables researchers to identify which genomic regions drive classification — linking ML features to biology.

## Computational Cost

| Dataset | Training Time | Model Size |
|---|---|---|
| 100 bacterial genomes | ~10 s | ~10 MB |
| 500 bacterial genomes | ~30 s | ~20 MB |
| 4000 synthetic sequences | ~2 s | ~2.5 KB |

Times measured on Apple M4 (4 cores). Training scales approximately linearly with dataset size due to the tree construction dominating.
