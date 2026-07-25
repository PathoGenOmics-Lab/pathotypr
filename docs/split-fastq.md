# `pathotypr split-fastq`

**Perform split-k-mer typing directly from FASTQ files — alignment-free genotyping straight from raw reads.**

`split-fastq` classifies samples into lineages without assembling or aligning anything. It builds a set of diagnostic reference/alternate k-mer pairs from your marker table, then streams the FASTQ reads past that index in constant memory, counting how often each allele is seen. Reach for it when you have raw sequencing reads (single- or paired-end, gzip is fine) and want a fast, low-footprint lineage call — for example screening many samples in a batch, or genotyping without waiting on an assembly step.

!!! note "SNP-based typing"
    FASTQ mode calls variants from short-read k-mer evidence and does **not** report indels — short reads give unreliable k-mer matches across insertions/deletions. For indel-aware, assembly-based typing use [`classify`](classify.md).

## Synopsis

```bash
pathotypr split-fastq -m <markers.tsv> -r <reference.fasta> (-i <reads.fq> ... | -l <list.tsv>) -o <prefix> [OPTIONS]
```

You must supply markers (`-m`), a reference (`-r`), and exactly one input source: one or more files with `-i/--input`, **or** a sample list with `-l/--input-list`.

## How it works

The marker database is built once and cached in memory, then reused for every sample in the run.

1. **Build** diagnostic k-mer pairs (REF and ALT) for each marker position, using the reference and marker table.
2. **Index** every marker k-mer (forward and reverse complement) behind a Bloom filter that rejects non-matching k-mers before the hash lookup.
3. **Scan** each sample's FASTQ reads in a streaming, constant-memory pass, counting REF vs. ALT k-mer hits per marker.
4. **Call** a variant at a marker when total depth ≥ `--min-depth` and the ALT fraction ≥ `--min-alt-percent`.
5. **Classify** the called markers into a lineage, optionally validating the hierarchical path when `--nested-classification` is set.
6. **Write** a per-sample mutations TSV, then a single run-wide summary TSV (plus `.xlsx` if `--excel`).

!!! tip "Paired-end and multi-chunk samples"
    Paired files are auto-detected from their names (`_1/_2`, `_R1/_R2`). Use `--paired` to force sequential pairing, or `--no-auto-paired` to treat every file as single-end. When a sample lists several R1/R2 files (e.g. multiple lanes), **all** of them are scanned together into one result.

## Options

### Inputs

| Option | Default | Required | Description |
|---|---|---|---|
| `-i, --input <FILE>...` | — | this or `--input-list` | One or more FASTQ files (gzip supported). Paired-end files are auto-detected by name (`_1/_2`, `_R1/_R2`). |
| `-l, --input-list <FILE>` | — | this or `--input` | TSV of samples: `sample_name\tR1.fq[\tR2.fq]`. |
| `--paired` | `false` | no | Force paired-end mode: input files are grouped sequentially in pairs. |
| `--no-auto-paired` | `false` | no | Disable auto paired-end detection; treat every file as single-end. |
| `-r, --reference <FILE>` | — | yes | Reference FASTA the markers refer to. |
| `-m, --markers <FILE>` | — | yes | Marker TSV: `position, REF, ALT, level1, level2...` |

### Calling and classification

| Option | Default | Required | Description |
|---|---|---|---|
| `--min-depth <N>` | `10` | no | Minimum read depth to call a variant at a marker position. |
| `--min-alt-percent <N>` | `95` | no | Minimum alternate-allele frequency (percent) to call a variant. |
| `--nested-classification` | `false` | no | Enable hierarchical (nested) lineage calling. |
| `-k, --kmer-size <N>` | `31` | no | Diagnostic k-mer size (1–31; odd values recommended so the variant base sits at the centre of the k-mer). |

### Output and runtime

| Option | Default | Required | Description |
|---|---|---|---|
| `-o, --output-prefix <PREFIX>` | `split` | no | Prefix for output files. |
| `--excel` | `false` | no | Also write Excel (`.xlsx`) files. |
| `-t, --threads <N>` | all cores | no | Number of CPU threads. |

### Global flags

Available on every subcommand:

| Flag | Description |
|---|---|
| `-v, --verbose` | Increase log verbosity. Repeatable: `-v` = debug, `-vv` = trace. |
| `-h, --help` | Print help. |

## Examples

=== "Paired-end (auto-detected)"

    Two files whose names differ only by `_R1`/`_R2` are recognised as one sample automatically.

    ```bash
    pathotypr split-fastq \
      -m markers.tsv \
      -r reference.fasta \
      -i reads_R1.fastq.gz -i reads_R2.fastq.gz \
      -o genotype
    ```

=== "Batch from a sample list"

    Process many samples in one run, giving each an explicit name. Each row is `sample_name`, then its R1 (and optional R2) path, tab-separated.

    ```bash
    pathotypr split-fastq \
      -m markers.tsv \
      -r reference.fasta \
      -l samples.tsv \
      -o batch_run \
      --nested-classification \
      --excel
    ```

    `samples.tsv`:

    ```tsv
    sampleA	sampleA_R1.fastq.gz	sampleA_R2.fastq.gz
    sampleB	sampleB_R1.fastq.gz	sampleB_R2.fastq.gz
    sampleC	sampleC.fastq.gz
    ```

=== "Force pairing + stricter calls"

    When filenames don't follow a known convention, `--paired` groups the inputs sequentially in pairs. Here we also raise the depth floor and use a smaller k-mer.

    ```bash
    pathotypr split-fastq \
      -m markers.tsv \
      -r reference.fasta \
      -i lane1_1.fq.gz -i lane1_2.fq.gz \
      --paired \
      --min-depth 20 \
      -k 21 \
      -o strict_run
    ```

!!! warning "`--paired` needs an even file count"
    In forced paired mode the files are consumed two at a time, so `-i` must list an even number of files. Mixed single/paired batches are better expressed with `--input-list`.

## Output

All files are prefixed with `--output-prefix` (default `split`). Sample names are sanitised for filesystem safety, and collisions are de-duplicated with a numeric suffix.

| File | When | Content |
|---|---|---|
| `<prefix>_<sample>_mutations.tsv` | one per sample | Every marker that passed the depth and ALT-fraction thresholds. |
| `<prefix>_summary.tsv` | once per run | One row per sample with lineage tallies and the major-lineage call. |
| `<prefix>_<sample>_mutations.xlsx` | with `--excel` | Excel copy of each per-sample mutations table. |
| `<prefix>_summary.xlsx` | with `--excel` | Excel copy of the summary table. |

### Per-sample mutations TSV

`<prefix>_<sample>_mutations.tsv` — one row per called marker.

| Column | Description |
|---|---|
| `pos` | 1-based marker position on the reference. |
| `ref_allele` | Reference allele. |
| `alt_allele` | Alternate allele. |
| `ref_count` | k-mer hits supporting the reference allele. |
| `alt_count` | k-mer hits supporting the alternate allele. |
| `alt_fraction` | ALT percentage, `alt_count / (ref_count + alt_count) × 100`, to 2 decimals. |
| `lineage_path` | Full lineage path for the marker (levels joined by `;`). |
| `extra_annotations…` | Optional trailing column(s): any extra marker annotations, only present when the marker defines them. |

!!! note "Excel conditional formatting"
    In the `.xlsx` copy the `alt_fraction` column is colour-coded: red below 50%, yellow 50–95%, green at 95% and above.

### Run summary TSV

`<prefix>_summary.tsv` — header and one row per sample.

```tsv
genome	lineage:count	major_lineage
```

| Column | Description |
|---|---|
| `genome` | Sample name. |
| `lineage:count` | Space-separated `lineage:count` pairs (each level and how many markers supported it), ordered by count then name. |
| `major_lineage` | The final major-lineage call for the sample. |

## See also

- [Input formats](input-formats.md) — accepted FASTQ inputs and sample-list layout.
- [Marker format](marker_format.md) — how to build the marker TSV (`position, REF, ALT, level1, level2…`).
- [`classify`](classify.md) — assembly-based, indel-aware typing.
- [Marker genotyping algorithm](algorithms/marker-genotyping.md) — diagnostic k-mer construction, the Bloom-filter pre-check, and the streaming scan.
