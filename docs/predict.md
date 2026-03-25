# `pathotypr predict`

Predict lineages for new genome sequences using a trained model.

## How it works

1. **Load** model bundle (zstd-compressed)
2. **Validate** model format version and integrity
3. **Stream** input FASTA in batches of 512
4. **Vectorize** each batch via the same feature hasher used during training
5. **Predict** in parallel — each tree votes, majority wins
6. **Write** results immediately (constant memory)

## Usage

```bash
pathotypr predict \
  --input query_genomes.fasta \
  --model model.pathotypr.zst \
  --output predictions.tsv \
  --threads 8 \
  --excel
```

## Options

| Flag | Default | Description |
|---|---|---|
| `-i, --input` | required | FASTA file with query sequences (supports .gz) |
| `-m, --model` | required | Trained model bundle (`.pathotypr.zst`) |
| `-o, --output` | required | Output TSV path |
| `-t, --threads` | all cores | Number of CPU threads |
| `--excel` | off | Also generate .xlsx alongside TSV |

## Output columns

| Column | Description |
|---|---|
| `Header` | Sequence header from input FASTA |
| `Predicted_Lineage` | Most-voted class label |
| `Confidence` | Fraction of trees voting for the winner (0–1) |
| `Confidence_Margin` | Gap between winner and runner-up (0–1) |
| `Other_Votes` | Top 3 alternative predictions with vote fractions |

## Technical details

- **Streaming**: processes 512 sequences at a time — RAM usage is O(batch), not O(N)
- **Parallel**: vectorization and tree voting are parallelized with rayon
- **Gzip support**: handles `.fasta.gz` transparently via needletail
- **Confidence interpretation**: 0.95 = 95% of trees agree; margin 0.40 = winner has 40% more votes than runner-up

## Algorithm Details

For in-depth documentation of the underlying algorithms:

- [Prediction](algorithms/prediction.md) — Streaming batch processing, ensemble voting, confidence and margin metrics
- [Random Forest](algorithms/random-forest.md) — Tree traversal via binary search on sparse rows, majority voting
- [Feature Hashing](algorithms/feature-hashing.md) — How sequences are vectorized identically at train and predict time
