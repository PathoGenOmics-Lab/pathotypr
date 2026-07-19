# `pathotypr match`

**Find the best-matching reference genome for a set of FASTQ reads using k-mer containment.**

`pathotypr match` counts the k-mers in your sequencing reads and scores them against a database of candidate reference genomes packed into a single multi-FASTA. It reports the single reference that shares the largest fraction of your reads' k-mers, together with a weighted containment fraction. Use it for quick reference selection or species/strain sanity checks before running a full classification.

!!! note "One query set, one best match"
    `match` pools **all** the FASTQ files you provide (whether via `--input` or every row of an `--input-list`) into a single query k-mer set and reports **one** best-matching reference. It does not produce per-sample rows. To keep samples separate, run `match` once per sample.

## Synopsis

```bash
pathotypr match -r <references.fasta> (-i <reads.fq> ... | -l <list.tsv>) [-o <report.tsv>] [OPTIONS]
```

You must supply exactly one input source (`-i/--input` or `-l/--input-list`) and a reference database (`-r/--references`).

## How it works

1. **Count** k-mers across every input FASTQ file in parallel, tallying occurrence counts.
2. **Filter** likely sequencing noise by dropping low-count k-mers (see `--min-kmer-count`).
3. **Stream** the reference multi-FASTA in batches, extracting the unique k-mers of each reference.
4. **Score** each reference by the weighted fraction of query k-mer occurrences it contains.
5. **Report** the highest-scoring reference and its `Shared_Kmer_Fraction`.

!!! note "Bounded memory"
    References are streamed and scored one batch at a time — only a single batch is resident at once. Peak memory is governed by the query k-mer set plus one batch of references, not by the total number of reference genomes in the database. Batch size adapts to the number of CPU threads.

## Options

| Option | Default | Required | Description |
|---|---|---|---|
| `-i, --input <FILE>...` | — | this or `--input-list` | One or more query FASTQ files (gzip supported). |
| `-l, --input-list <FILE>` | — | this or `--input` | TSV listing FASTQ files: `sample_name\treads1.fastq[\treads2.fastq]`. |
| `--paired` | `false` | no | Force paired-end mode. |
| `--no-auto-paired` | `false` | no | Disable auto paired-end detection. |
| `-r, --references <FILE>` | — | **yes** | Single multi-FASTA containing all reference genomes. |
| `-o, --output <FILE>` | stdout | no | Output TSV report. Prints to the console if omitted. |
| `-k, --kmer-size <N>` | `31` | no | k-mer size for comparison (1–31). |
| `-t, --threads <N>` | all cores | no | Number of CPU threads. |
| `--early-stop-confidence <FLOAT>` | `0.0` | no | Reserved; accepted but not read in the current release (no effect). |
| `--early-stop-min-kmers <N>` | `1000000` | no | Reserved; accepted but not read in the current release (no effect). |
| `--strict-percentages` | `true` | no | Legacy flag retained for backward compatibility (default on); it does not change current results. |
| `--min-kmer-count <N>` | `2` | no | Minimum k-mer count to keep (singleton filtering as noise). Set `1` to keep all. Only applied when there are more than 100k unique k-mers and this value is greater than `1`. |
| `--excel` | `false` | no | Also write an Excel (`.xlsx`) file. |

### Global flags

Available on every subcommand:

| Flag | Description |
|---|---|
| `-v, --verbose` | Increase log verbosity. Repeatable: `-v` = debug, `-vv` = trace. |
| `-h, --help` | Print help. |

!!! tip "Paired-end handling"
    When you pass more than one file to `--input`, `match` auto-detects paired-end files from common naming conventions (`_1`/`_2`, `_R1`/`_R2`). Because `match` pools all reads into one query set, this grouping does not change the reported result. Use `--paired` to force pairing (it requires an even number of files) or `--no-auto-paired` to switch detection off.

## Examples

=== "Paired reads"

    ```bash
    pathotypr match \
      -r references.fasta \
      -i reads_R1.fastq.gz reads_R2.fastq.gz \
      -o report.tsv
    ```

    Counts k-mers from both mates, scores them against every genome in `references.fasta`, and writes the best match to `report.tsv`.

=== "Sample list"

    ```bash
    pathotypr match \
      -r references.fasta \
      -l samples.tsv \
      -o report.tsv
    ```

    Reads every FASTQ path listed in `samples.tsv` (`sample_name\treads1.fastq[\treads2.fastq]`) and pools them into a single query. See [Input formats](input-formats.md#sample-list-tsv-split-fastq-and-match-input-list) for the TSV layout.

=== "Excel + tuning"

    ```bash
    pathotypr match \
      -r references.fasta \
      -i reads_R1.fastq.gz reads_R2.fastq.gz \
      -o report.tsv \
      -k 27 -t 8 --excel
    ```

    Uses a k-mer size of 27 on 8 threads and additionally writes `report.xlsx` beside the TSV.

!!! warning "`--excel` needs `-o`"
    The Excel file is written only when an output path is given with `-o/--output`; its name is derived from the TSV path (`report.tsv` → `report.xlsx`). If output goes to the console (no `-o`), `--excel` has no effect.

## Output

### TSV report

The report has a header row followed by exactly one data row. Columns are tab-separated:

```text
Query_Files	Best_Match_Reference	Shared_Kmer_Fraction
```

| Column | Description |
|---|---|
| `Query_Files` | Comma-separated list of the input FASTQ file paths that were pooled into the query. |
| `Best_Match_Reference` | Identifier (first whitespace-delimited token of the FASTA header) of the highest-scoring reference. |
| `Shared_Kmer_Fraction` | Weighted k-mer containment score between `0` and `1`, printed to four decimal places. Higher means more of the query's k-mer content is present in the matched reference. |

Example:

```text
Query_Files	Best_Match_Reference	Shared_Kmer_Fraction
reads_R1.fastq.gz,reads_R2.fastq.gz	NC_000962.3	0.9873
```

If no matching reference is found, the header is still written and a message is logged; no data row is emitted.

### Excel report (optional)

With `--excel` **and** `-o`, `match` also writes an `.xlsx` workbook whose name is the TSV path with `.tsv` replaced by `.xlsx` (otherwise `.xlsx` is appended). It contains the same three columns as the TSV, with a frozen, styled header row and an autofilter.

## See also

- [Input formats](input-formats.md) — FASTQ inputs, the sample-list TSV, and the reference multi-FASTA.
- [Reference matching algorithm](algorithms/reference-matching.md) — k-mer containment scoring, noise filtering, and streaming batch processing.
- [`pathotypr classify`](classify.md) — marker-based lineage classification once a reference is chosen.
- [`pathotypr split-fastq`](split-fastq.md) — marker-based classification straight from FASTQ reads.
- [Benchmarks](benchmarks.md) — runtime and memory characteristics.
