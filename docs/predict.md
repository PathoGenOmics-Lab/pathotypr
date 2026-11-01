# `pathotypr predict`

**Assign a lineage to each genome in a FASTA file using a model produced by `pathotypr train`.**

`pathotypr predict` classifies query sequences against a pre-trained random-forest model bundle. Each record receives a majority-vote lineage call together with a confidence score, a confidence margin, and the runner-up votes. Use it once you have a trained `.pathotypr.zst` model and want to genotype new assemblies or genomes — the input FASTA is streamed in batches, so memory stays constant even for very large collections.

## Synopsis

```bash
pathotypr predict -i <query.fasta> -m <model.pathotypr.zst> -o <predictions.tsv> [OPTIONS]
```

## Inputs

`predict` needs two files.

| What | Flag | Format |
|---|---|---|
| Genomes to classify | `-i, --input` | FASTA, plain or gzipped. Multi-record is fine, and **each record is classified separately** |
| Trained model | `-m, --model` | A `.pathotypr.zst` bundle written by [`train`](train.md) |

### Requirements

- [x] The model must be a **current-format bundle**. An incompatible format version is a hard error and the model has to be retrained. A bundle written by a different pathotypr version only warns.
- [x] Records that yield **no k-mer**, because they are shorter than the model's k or carry no long enough window of A, C, G and T, are reported as `Undetermined` with `Kmers_Used` zero rather than dropped.
- [x] No labels are needed or read. The **full header** is carried through to the output.

!!! warning "The model dictates k, not you"
    `predict` has no `-k`. The k-mer size is whatever `train` used and is stored
    in the bundle, so a model trained at `k=21` always vectorises queries at
    `k=21`. This is also why the minimum usable sequence length depends on the
    model rather than on the query file.

!!! info "Any record with k-mers gets a call"
    A record that yields no k-mer is `Undetermined`. Every other record gets a
    label, because there is no "unknown" class and no confidence threshold: a
    query from a lineage absent from the training set still receives the closest
    label the forest can find, distinguishable only by a **low confidence and
    margin**, and by `Kmers_Used` when it rests on very few k-mers.
    Read those two columns before acting on a call: see
    [Output](#predictionstsv-always-written) for what they mean.

## Options

| Option | Default | Required | Description |
|---|---|---|---|
| `-i, --input <FILE>` | — | yes | Input FASTA of sequences to classify (gzip supported). |
| `-m, --model <FILE>` | — | yes | Model bundle produced by `pathotypr train`. |
| `-o, --output <FILE>` | — | yes | Output predictions TSV. |
| `-t, --threads <N>` | all cores | no | Number of CPU threads. |
| `--excel` | `false` | no | Also write an Excel (`.xlsx`) file alongside the TSV. |

### Global flags

These are available on every `pathotypr` subcommand:

| Flag | Description |
|---|---|
| `-v, --verbose` | Increase log verbosity. Repeatable: `-v` = debug, `-vv` = trace. |
| `-h, --help` | Print help. |

!!! note "The model carries its own settings"
    The k-mer size and feature-hashing parameters are stored inside the model bundle. You do not (and cannot) set them at predict time — the same vectorization used during training is applied automatically, so predictions stay consistent with the model.

## How it works

1. **Load** the model bundle and decompress it (zstd).
2. **Validate** the model bundle: the run stops only if its *format* version differs from the one this build expects. A model trained by a different `pathotypr` release with the same format loads fine, with a warning.
3. **Stream** the input FASTA in batches of 512 records — RAM usage is proportional to the batch, not the number of sequences.
4. **Vectorize** each batch with the exact feature hasher stored in the model.
5. **Predict** in parallel: every tree in the forest votes, and the majority class becomes the call.
6. **Score** each record with a confidence and a confidence margin, then write results immediately and free the batch.

!!! warning "Short sequences are skipped"
    A record that yields no k-mer is reported as `Undetermined` with `Kmers_Used` zero, and the run says at the end how many came out that way. It is not dropped: a sample missing from the output is a sample nobody notices is missing. This covers records shorter than the model's k and, less obviously, N-masked consensus sequences, which are long but hold no window of A, C, G and T: one `N` destroys every k-mer overlapping it, so at k=21 an `N` every twenty bases leaves nothing to read. Those used to come out with an ordinary-looking lineage and confidence, because a tree sends a missing feature down the same branch every time and all the trees then agree. Only a file with no records at all is an error.

## Examples

=== "Basic"

    ```bash
    pathotypr predict \
      -i sample.fasta \
      -m model.pathotypr.zst \
      -o predictions.tsv
    ```

=== "Gzipped input + Excel"

    ```bash
    # .gz is detected and decompressed transparently
    pathotypr predict \
      -i query.fasta.gz \
      -m model.pathotypr.zst \
      -o predictions.tsv \
      --excel
    ```

=== "Limit threads + verbose"

    ```bash
    pathotypr predict \
      -i sample.fasta \
      -m model.pathotypr.zst \
      -o predictions.tsv \
      --threads 8 \
      -v
    ```

!!! tip "Reading confidence"
    `Confidence` is the fraction of trees that voted for the winning lineage (e.g. `0.9500` means 95% of trees agreed). `Confidence_Margin` is how far ahead the winner is over the runner-up (e.g. `0.4000` means the winner led the second-best class by 40% of the trees). A high confidence with a low margin points to two closely competing lineages — inspect `Other_Votes`.

## Output

### `predictions.tsv` (always written)

A tab-separated file with a header row and one row per classified record. Columns are exactly:

| Column | Description |
|---|---|
| `Header` | The input record's full FASTA header line — everything after `>`, description included. |
| `Predicted_Lineage` | Majority-vote lineage label, or `Undetermined` when the record yielded no k-mer and the model therefore said nothing about it. A tie between two lineages is settled by the label, so the answer does not depend on the order the training FASTA happened to be in. |
| `Confidence` | Fraction of trees voting for the winning lineage, in `0.0000`–`1.0000` (4 decimals). |
| `Confidence_Margin` | Winner-minus-runner-up vote fraction, in `0.0000`–`1.0000` (4 decimals). |
| `Other_Votes` | Up to the three next-best lineages as `label:fraction` (2 decimals), comma-separated. Empty when there are no alternatives. |
| `Kmers_Used` | How many k-mers the record yielded. Zero is what `Undetermined` rests on, and a handful tells you a call rests on almost nothing, which no confidence figure shows: every tree walks the same branch when a feature is missing, so they agree, and the agreement reads as confidence. |

Example:

```text
Header	Predicted_Lineage	Confidence	Confidence_Margin	Other_Votes	Kmers_Used
sample_001	L4	0.9600	0.9300	L2:0.03,L1:0.01	4411512
sample_002	L2	0.5400	0.0800	L4:0.46	4398877
masked_003	Undetermined				0
```

### `predictions.xlsx` (only with `--excel`)

When `--excel` is set, an Excel workbook with the **same five columns** is written next to the TSV (the output path with an `.xlsx` extension). Conditional formatting is applied to the `Confidence` and `Confidence_Margin` columns to highlight lower-scoring calls at a glance.

!!! note "Excel is best-effort"
    The TSV is the authoritative output. If the Excel writer cannot be initialized or a row fails to append, `predict` logs a warning and continues writing the TSV normally.

## See also

- [Train](train.md) — build the `.pathotypr.zst` model bundle used here.
- [Input formats](input-formats.md) — accepted FASTA/gzip inputs.
- [Classify](classify.md) — marker-based assembly classification.
- [Prediction algorithm](algorithms/prediction.md) — streaming batches, ensemble voting, confidence and margin metrics.
- [Random forest](algorithms/random-forest.md) — how each tree votes.
- [Feature hashing](algorithms/feature-hashing.md) — how sequences are vectorized identically at train and predict time.
