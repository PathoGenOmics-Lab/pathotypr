# Prediction (streaming batch)

**How the `predict` command turns a FASTA of unknown genomes into lineage calls: a streaming, constant-memory pipeline that vectorizes each batch and aggregates votes across the trained forest.**

The `predict` command classifies new genome sequences using a trained model. It uses **streaming batch processing** to handle arbitrarily large inputs with constant peak memory.

## Pipeline

```text
Input FASTA ──→ Stream in batches of 512 ──→ Vectorize (parallel) ──→ Predict (parallel) ──→ Write TSV
     │                                                                                         │
     └── needletail parser                                                               flush & free
         (handles .gz transparently)                                                     batch memory
```

## Step 1: Model loading

```text
.pathotypr.zst file ──→ zstd::decode_all ──→ bincode::deserialize ──→ ModelBundle
```

**Validation checks** on load:

- Trees array is non-empty
- Label encoder has at least one class
- K-mer size is within 1–31
- Format version matches current (v3)

!!! warning "Version compatibility"
    An incompatible **format version** is a hard error — the model must be retrained with this version of pathotypr. A model trained by a **different pathotypr version** is not rejected: it only produces a warning (not an error), since results may vary.

## Step 2: Streaming batches

Instead of loading all sequences into memory:

1. Open FASTA with `needletail::parse_fastx_file()` (handles plain and gzipped)
2. Read up to 512 records into a batch
3. Process the batch (vectorize + predict)
4. Write results to TSV
5. **Drop the batch** — memory is freed before the next batch
6. Repeat until EOF

!!! note "Why 512?"
    Small enough for constant memory (~50 MB per batch for bacterial genomes), large enough for efficient rayon parallelism.

## Step 3: Vectorization

Each batch of sequences is vectorized in parallel using the same `FeatureHasher` from training:

```rust
let x_sparse = bundle.vectorizer.transform_sparse(&sequences, kmer_size);
```

This produces one sparse vector per sequence, identical to the training representation.

## Step 4: Ensemble voting

For each sparse vector, every tree in the ensemble (100 by default) votes independently:

```text
Tree 1: predict_one(sparse_row) → class 2
Tree 2: predict_one(sparse_row) → class 2
Tree 3: predict_one(sparse_row) → class 0
...
Tree 100: predict_one(sparse_row) → class 2
```

Each `predict_one()` traverses root-to-leaf in O(depth) using binary search on the sparse row.

### Metrics

| Metric | Formula | Meaning |
|---|---|---|
| **Predicted class** | argmax(votes) | Class with most votes |
| **Confidence** | winner_votes / n_trees | Proportion of trees agreeing |
| **Margin** | (winner − runner_up) / n_trees | Separation from second-best class |
| **Other votes** | Top 3 non-winning classes as `label:fraction` (2 decimals) | Alternative classifications |

Here `n_trees` is the number of trees in the loaded model (`bundle.trees.len()`, 100 by default), not a hard-coded constant.

## Output format

### TSV

```tsv
Header                Predicted_Lineage  Confidence  Confidence_Margin  Other_Votes
sample_0001           L4                 0.9800      0.9200             L4.9:0.03,L2:0.01
sample_0002           L2                 1.0000      1.0000
```

### Excel (optional, `--excel`)

Same columns with conditional formatting:

- **Confidence**: green ≥ 0.9, yellow 0.5–0.9, red < 0.5
- **Margin**: green ≥ 0.5, yellow 0.1–0.5, red < 0.1

## Performance

- **Throughput**: ~85,000 genomes/second at 4000 genomes (synthetic)
- **Prediction time per sample**: ~0.03 ms (constant regardless of dataset size)
- **Memory**: Constant peak — batch memory is freed between iterations
- **I/O**: SIMD-accelerated gzip decompression via zlib-ng

## Error handling

- Sequences shorter than the model's k-mer size (including empty ones) are skipped with a warning — no k-mer can be extracted from them
- Non-UTF8 headers/sequences produce parsing errors
- Cancellation is checked between batches (GUI integration)
- On error, partial output files are cleaned up

## See also

- [predict](../predict.md)
- [Random forest](random-forest.md)
