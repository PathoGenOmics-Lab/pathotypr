# `pathotypr train`

**Train a new Random Forest model on a set of labeled genomes and save it as a single compressed bundle.**

`pathotypr train` reads a labeled FASTA, turns each sequence into a fixed-size k-mer feature vector via feature hashing, and fits an ensemble of 100 decision trees. It estimates accuracy (either a single held-out split or stratified cross-validation), then trains the **final** model on **all** of your data and writes it to disk together with feature-importance reports. Use it whenever you have a curated, labeled collection of genomes and want a model that [`pathotypr predict`](predict.md) can later apply to new sequences.

## Synopsis

```bash
pathotypr train -i <labeled.fasta> -o <model.pathotypr.zst> [OPTIONS]
```

The label for each record is the **first whitespace-separated token** of its FASTA header. Everything after the first space is ignored for labeling. See [Input formats](input-formats.md) for details.

```text
>L2.2.1 sample_A collected_2021
ACGTACGT...        →  label = "L2.2.1"
>L4.3.4.2 sample_B
ACGTACGT...        →  label = "L4.3.4.2"
```

## Options

| Option | Default | Required | Description |
|---|---|---|---|
| `-i, --input <FILE>` | — | yes | Labeled training FASTA. The first whitespace-separated token in each header is the class label. |
| `-o, --output <FILE>` | — | yes | Path to write the trained model bundle (bincode + zstd). |
| `-k, --kmer-size <N>` | `21` | no | k-mer size used for feature hashing (1–31). |
| `-s, --test-split <FLOAT>` | `0.2` | no | Fraction of data held out for the single train/test accuracy estimate. |
| `-t, --threads <N>` | all cores | no | Number of CPU threads. |
| `--cv-folds <N>` | off | no | Use stratified k-fold cross-validation for accuracy instead of a single split. The final model is always trained on ALL data. |
| `--max-depth <N>` | `20` | no | Maximum tree depth. Lower reduces overfitting on small datasets. |
| `--min-samples-leaf <N>` | `5` | no | Minimum samples per leaf; higher acts as regularization. |

### Global flags

Available on every subcommand:

| Flag | Description |
|---|---|
| `-v, --verbose` | Increase log verbosity. Repeatable: `-v` = debug, `-vv` = trace. |
| `-h, --help` | Print help and exit. |

!!! note "Validation rules"
    - `--test-split` must be in the range `[0.0, 1.0)`.
    - `--cv-folds`, when set, must be at least `2` and no greater than the number of samples.
    - `--kmer-size` must be between `1` and `31`.
    - The training data must contain **at least two distinct labels**.

## How it works

1. **Read** the labeled FASTA. The first token of each header becomes the class label; records with empty sequences are skipped.
2. **Vectorize** every sequence with feature hashing into a fixed **2²⁰ = 1,048,576**-bucket sparse vector. There is no vocabulary to fit, so train- and predict-time features are identical.
3. **Encode** labels to integer classes.
4. **Estimate accuracy** with one of two strategies:
    - *Single split* (default): hold out `--test-split` of the data, train an evaluation ensemble on the rest, and score it once.
    - *Cross-validation* (`--cv-folds N`): run stratified k-fold CV and report mean ± standard deviation.
5. **Train the final model** on **100%** of the data — 100 trees, each fit on a bootstrap sample (~63% in-bag, ~37% out-of-bag), with `ceil(sqrt(n_features))` candidate features per split.
6. **Compute out-of-bag (OOB) accuracy** from the bootstrap samples — a nearly unbiased estimate that costs nothing extra and is always reported.
7. **Serialize** the model bundle (config + feature hasher + label encoder + 100 trees) with bincode and stream it through zstd (level 3) to the `--output` path.
8. **Export** feature-importance reports and the genomic coordinates of the most discriminant k-mers next to the model file.

!!! warning "Small datasets"
    If the input contains fewer than 10 sequences, training still runs but emits a warning — the resulting model may not be reliable. For small or imbalanced sets, prefer `--cv-folds` for a more stable accuracy estimate, and consider a lower `--max-depth` and higher `--min-samples-leaf` to curb overfitting.

## Examples

=== "Basic training"

    Train with all defaults (k = 21, single 80/20 split, all cores):

    ```bash
    pathotypr train \
      -i references.fasta \
      -o model.pathotypr.zst
    ```

=== "Small / imbalanced dataset"

    Use 5-fold stratified cross-validation and tighter regularization to reduce overfitting:

    ```bash
    pathotypr train \
      -i references.fasta \
      -o model.pathotypr.zst \
      --cv-folds 5 \
      --max-depth 12 \
      --min-samples-leaf 8
    ```

=== "Tuned k-mer and threads"

    Increase k-mer size, hold out 30% for the split estimate, and cap CPU usage:

    ```bash
    pathotypr train \
      -i references.fasta \
      -o model.pathotypr.zst \
      --kmer-size 27 \
      --test-split 0.3 \
      --threads 8 \
      -v
    ```

!!! tip "Reproducibility"
    Data shuffling, fold assignment, and per-tree bootstraps use fixed seeds, so repeated runs on the same input and options produce the same model and accuracy figures.

## Output

`train` writes the model plus two tab-separated report files. The report filenames are derived from `--output` by replacing its final extension. For an output of `model.pathotypr.zst`:

| File | Content |
|---|---|
| `model.pathotypr.zst` | The trained model bundle: config (pathotypr version, k-mer size, tree count, format version), feature hasher, label encoder, and all 100 trees — serialized with bincode and zstd-compressed (level 3). This is the file you pass to [`pathotypr predict`](predict.md). |
| `model.pathotypr.importance.tsv` | The top 500 most-used feature-hash buckets across the ensemble, with the k-mers that map to each. |
| `model.pathotypr.importance.coords.tsv` | Per-occurrence genomic coordinates of the discriminant k-mers behind those top buckets. |

!!! note "Accuracy is reported to the terminal, not written to a file"
    The single-split (or cross-validated) accuracy and the out-of-bag accuracy are logged during the run. Per-fold accuracies print at the default verbosity; `-v` adds further debug detail.

### `*.importance.tsv` columns

Ranked by how often each bucket was chosen as a split across all 100 trees (highest first), limited to the top 500 buckets.

| Column | Description |
|---|---|
| `rank` | 1-based rank of the bucket by split count. |
| `bucket` | Feature-hash bucket index (`0 … 1,048,575`). |
| `split_count` | Number of times this bucket was used as a split point across the ensemble. |
| `importance_pct` | `split_count` as a percentage of the total split count over all buckets. |
| `kmers` | Comma-separated k-mers from the training sequences that hash into this bucket. |

### `*.importance.coords.tsv` columns

One row per occurrence of a discriminant k-mer in the training sequences, ordered by bucket rank, then sequence, then position.

| Column | Description |
|---|---|
| `rank` | Importance rank of the bucket this k-mer maps to. |
| `bucket` | Feature-hash bucket index. |
| `split_count` | Split count for the bucket. |
| `importance_pct` | Percentage of total splits attributable to the bucket. |
| `kmer` | The specific k-mer occurrence. |
| `sequence` | Full FASTA header of the training sequence containing the k-mer. |
| `lineage` | Class label of that sequence. |
| `position` | Position of the k-mer within the sequence. |

## See also

- [`pathotypr predict`](predict.md) — Apply a trained model to new genome sequences.
- [Input formats](input-formats.md) — How labeled FASTA headers and labels are parsed.
- [Feature hashing](algorithms/feature-hashing.md) — The hashing trick and why train/predict features match.
- [Random forest](algorithms/random-forest.md) — Sparse decision trees, bootstrap aggregation, and OOB accuracy.
- [Training pipeline](algorithms/training.md) — End-to-end training internals.
