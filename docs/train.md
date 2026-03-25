# `pathotypr train`

Train a Random Forest classifier from labeled genome sequences.

## How it works

1. **Read** labeled FASTA (first header token = class label)
2. **Vectorize** sequences via feature hashing (2²⁰ buckets, hashing trick)
3. **Evaluate** accuracy via train/test split or stratified k-fold CV
4. **Train** the final model on **all** data (100 trees, bootstrap sampling)
5. **Compute** out-of-bag (OOB) accuracy from the bootstrap samples
6. **Save** compressed model bundle (bincode + zstd)
7. **Export** feature importance with genomic coordinates

## Usage

```bash
pathotypr train \
  --input labeled_genomes.fasta \
  --output model.pathotypr.zst \
  --kmer-size 21 \
  --test-split 0.2 \
  --threads 8 \
  --excel
```

## Options

| Flag | Default | Description |
|---|---|---|
| `-i, --input` | required | Labeled FASTA file |
| `-o, --output` | required | Output model path (`.pathotypr.zst`) |
| `-k, --kmer-size` | `21` | K-mer size for feature hashing (1–31) |
| `-s, --test-split` | `0.2` | Fraction held out for accuracy estimation |
| `-t, --threads` | all cores | Number of CPU threads |
| `--cv-folds` | — | Stratified k-fold CV instead of single split (e.g., 5 or 10) |
| `--max-depth` | `20` | Maximum tree depth (regularization) |
| `--min-samples-leaf` | `5` | Minimum samples per leaf node (regularization) |

## Output files

| File | Content |
|---|---|
| `model.pathotypr.zst` | Compressed model bundle (vectorizer + encoder + 100 trees) |
| `model.importance.tsv` | Top 500 features: rank, bucket, split count, k-mers |
| `model.importance.coords.tsv` | Genomic coordinates for discriminant k-mers |

## Accuracy estimation

- **Single split** (default): trains on 80%, tests on 20%
- **k-fold CV** (`--cv-folds`): stratified by class, reports mean ± std
- **OOB accuracy**: always computed from bootstrap — free, nearly unbiased

The final model is always trained on 100% of the data regardless of evaluation method.

## Technical details

- **Feature hashing**: 2²⁰ = 1,048,576 buckets via bitmask (no vocabulary)
- **Trees**: 100 decision trees, sqrt(n_features) candidates per split, Gini impurity
- **Bootstrap**: each tree trained on ~63% of samples; ~37% are out-of-bag
- **Serialization**: bincode → zstd streaming compression
- **Memory**: model size depends on tree complexity, typically 5–50 MB compressed

## Algorithm Details

For in-depth documentation of the underlying algorithms:

- [Feature Hashing](algorithms/feature-hashing.md) — The hashing trick, 2-bit encoding, bucket collisions, and reverse mapping
- [Random Forest](algorithms/random-forest.md) — Sparse CART trees, binary search on sparse rows, bootstrap aggregation, and OOB accuracy
- [Training Pipeline](algorithms/training.md) — End-to-end flow: vectorize → evaluate (CV or split) → train final model → serialize with bincode+zstd
