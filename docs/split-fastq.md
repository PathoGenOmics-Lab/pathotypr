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

## Inputs

| What | Flag | Format |
|---|---|---|
| Marker panel | `-m, --markers` | [Marker TSV](marker_format.md), tab-separated |
| Reference genome | `-r, --reference` | **Single-record** FASTA the marker positions are numbered against |
| Reads | `-i, --input` or `-l, --input-list` | FASTQ, plain or gzipped |

`-i` takes one or more files that together form **one** sample. `-l` takes a
[sample list](input-formats.md#sample-list-tsv-split-fastq-and-match-input-list)
where each row is an **independent** sample with its own report: `sample_name`
followed by one or more FASTQ paths.

### Requirements

The run **fails** unless all of these hold:

- [x] The reference is a **single record**, same rule as [`classify`](classify.md).
- [x] Marker positions are **1-based on that exact reference**.
- [x] `-k` between **1 and 31**.
- [x] Sample names in `-l` are **unique**; a duplicate aborts the run.
- [x] Every FASTQ path exists, checked before scanning starts.

!!! warning "Indel markers are skipped here"
    Any marker whose REF and ALT differ in length is dropped, and the count is
    logged as `indels skipped`. Short reads across repetitive regions (PE/PPE
    families, IS elements) give unreliable k-mer matches for insertions and
    deletions. If your panel leans on indels, type assemblies with
    [`classify`](classify.md) instead, which has full-length context.

!!! warning "Use plain ACGT alleles"
    Marker k-mers are 2-bit encoded. A marker containing an ambiguity code or
    any other non-ACGT character is **never indexed**, so it silently becomes
    unmatchable rather than raising an error.

### Depth and allele fraction

Two thresholds decide whether a marker is called, and they are the ones worth thinking about before a run:

| Option | Default | What it means |
|---|---|---|
| `--min-depth` | `10` | REF + ALT k-mer hits needed at a position before it is called at all |
| `--min-alt-percent` | `95` | Share of those hits that must carry the ALT allele |

!!! info "The default is deliberately strict about mixtures"
    At `95`, a position sitting at 60% ALT is **not called**, so intermediate
    frequencies never reach the output. That is what you want for a clean
    consensus call, and it is why mixed infections surface as **fixed markers on
    incompatible lineage branches** rather than as intermediate fractions. Lower
    it only if you intend to inspect within-sample frequencies, and expect
    sequencing noise to come with them.

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
| `--no-anchor-rescue` | `false` | no | Skip the second pass that reads a marker's locus from its anchors when its k-mers could not. Faster, and those markers are then reported as not determined. |
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

## Reading a marker its k-mers cannot

Both k-mers of a marker span the variant and its flanks, so **one other change
anywhere in the same 31 bases stops them both matching**. The marker is then
unreadable rather than absent, which with a resistance panel is the difference
between "no marker found" and "this sample is susceptible".

The `k-1` bases on each side of the variant are therefore indexed on their own
as well. A marker whose k-mers matched nothing is read again by counting bases
from whichever side still matches:

| What the anchors read | The marker is |
|---|---|
| the alternate allele, at `--min-depth` or more | present, with `evidence` naming the side |
| the reference, or any other base | absent |
| too few reads, or none | not determined, `low_depth` or `not_found` |
| two sides that disagree | not determined, `anchors_disagree` |

An anchor is one base shorter than a k-mer and says nothing about *where* it
matched, so only anchors that occur **exactly once** in the reference are used,
counting both strands. A marker left with no unique anchor is reported as
`no_unique_anchor` rather than guessed at.

For the same reason a marker whose ALT k-mer the reference already carries, or
whose window it carries twice, is never settled on read counts at all: reads
from the other copy count towards it just as well. Five and sixty-seven markers
of the WHO catalogue are in that position, and only their anchors can speak for
them.

!!! example "rpoB S450L behind a neighbour"
    A sample carrying S450L and a synonymous change seven bases away matches
    neither k-mer of the marker: 124 reads cover the position, every one of them
    showing the variant, and the k-mers see none of it. The left anchor sits
    clear of the neighbour, the 112 reads that carry it whole are read, and the
    marker is called with `evidence` = `left_anchor`.

This costs a second pass over the reads, which `--no-anchor-rescue` turns off.
On a 6x whole-genome read set made from H37Rv itself, the run took about 1.5
times as long with the pass, called nothing in either case, and settled all 72
markers whose words the reference duplicates, which no read count can settle:
262 markers were left unread without the pass and 190 with it, all of those at
the edge of the reads' coverage.

## Output

All files are prefixed with `--output-prefix` (default `split`). Sample names are sanitised for filesystem safety, and collisions are de-duplicated with a numeric suffix.

| File | When | Content |
|---|---|---|
| `<prefix>_<sample>_mutations.tsv` | one per sample | Every marker that passed the depth and ALT-fraction thresholds. |
| `<prefix>_summary.tsv` | once per run | One row per sample with lineage tallies and the major-lineage call. |
| `<prefix>_<sample>_undetermined.tsv` | one per sample | Every marker the reads could answer neither way. Only the header when they all were. |
| `<prefix>_<sample>_mutations.xlsx` | with `--excel` | Excel copy of each per-sample mutations table. |
| `<prefix>_<sample>_undetermined.xlsx` | with `--excel` | Excel copy of each per-sample undetermined table. |
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
| `evidence` | What the call rests on: `kmer` when the marker's own k-mer matched, or `left_anchor`, `right_anchor` or `both_anchors` when its locus was read by counting from one side. |
| `extra_annotations…` | Optional trailing column(s): any extra marker annotations, only present when the marker defines them. |

!!! note "Excel conditional formatting"
    In the `.xlsx` copy the `alt_fraction` column is colour-coded: red below 50%, yellow 50–95%, green at 95% and above.

### Per-sample undetermined TSV

`<prefix>_<sample>_undetermined.tsv` holds one row per marker the reads could
not answer, with the same columns as the mutations table plus a `reason`.

| `reason` | What happened |
|---|---|
| `not_found` | Nothing read the locus: neither k-mer matched and no anchor was read either. |
| `low_depth` | The locus was read, by fewer than `--min-depth` reads. |
| `no_unique_anchor` | The k-mers matched nothing and every anchor of the marker occurs more than once in the reference, so none of them can say which copy a read came from. |
| `anchors_disagree` | The two anchors were both read deeply enough and they do not agree, so neither is believed. |
| `alt_kmer_in_reference` | The reference itself carries this marker's ALT k-mer, so reads matching it say nothing about this locus, and its anchors could not settle it either. |
| `ref_kmer_repeated` | The reference carries this marker's window more than once, so reads matching it may come from the other copy. |

The counts in the row are the ones the decision rests on: the k-mer counts when
the anchors were not consulted, and the counts from the anchor that answered
when they were.

A marker with no reads on either allele is not a marker that is absent, and
this is where the difference is written down. Either way the call is "not read",
and the row says so instead of letting the marker fall silently into the
reference pile.

Markers the reads *did* read are not listed, even when nothing was called for
them. Both k-mers come from the same window, so a read matching either one
reads every base of that window: a position covered by a neighbour's reads is
answered for every marker that lies inside it. At rpoB 450 the WHO catalogue
carries ten markers; one of them being called reads the window for all ten, and
the other nine are absent, not unread.

!!! warning "Read it before calling a sample susceptible"
    A drug with no marker called but markers undetermined has not been shown
    susceptible. The summary carries the count so it is visible without opening
    the table.

### Run summary TSV

`<prefix>_summary.tsv` — header and one row per sample.

```tsv
genome	lineage:count	major_lineage	undetermined	undetermined:count
```

| Column | Description |
|---|---|
| `genome` | Sample name. |
| `lineage:count` | Space-separated `lineage:count` pairs (each level and how many markers supported it), ordered by count then name. |
| `major_lineage` | The final major-lineage call for the sample. With a WHO drug-resistance panel, the drugs with a grade 1 or 2 marker instead: the catalogue's columns are a drug and a grade, not a hierarchy, so tallying them as lineage levels answers a question nobody asked. |
| `undetermined` | How many markers the reads could answer neither way. |
| `undetermined:count` | Those markers as `label:count` pairs under the first level of their path (the drug for a resistance panel, the main lineage for a lineage panel). Empty when every marker was read. |

## See also

- [Input formats](input-formats.md) — accepted FASTQ inputs and sample-list layout.
- [Marker format](marker_format.md) — how to build the marker TSV (`position, REF, ALT, level1, level2…`).
- [`classify`](classify.md) — assembly-based, indel-aware typing.
- [Marker genotyping algorithm](algorithms/marker-genotyping.md) — diagnostic k-mer construction, the Bloom-filter pre-check, and the streaming scan.
