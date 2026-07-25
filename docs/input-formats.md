# Input Format Reference

pathotypr reads a small set of plain-text formats. Every file is **UTF-8**; every tabular file is **tab-separated** (`\t`) with **no spaces around the delimiters**. This page is the authoritative reference for each one — which command consumes it, the exact column semantics, and the parsing rules the tools apply.

!!! tip "Looking for one command's inputs?"
    This page is organised by **format**. Each command page also carries an
    **Inputs** section listing exactly what that command needs and what it
    rejects: [train](train.md#inputs) · [predict](predict.md#inputs) ·
    [classify](classify.md#inputs) · [split-fastq](split-fastq.md#inputs) ·
    [match](match.md#inputs).

!!! abstract "Formats at a glance"

    | Format | Consumed by | Flag |
    |---|---|---|
    | [Training FASTA](#training-fasta-train) | [`train`](train.md) | `-i, --input` |
    | [Query FASTA](#training-fasta-train) | [`classify`](classify.md) | `-i, --input` / `--input-files` / `-l, --input-list` |
    | [Query FASTA](#training-fasta-train) | [`predict`](predict.md) | `-i, --input` |
    | [Marker TSV](#marker-tsv-classify-and-split-fastq) | [`classify`](classify.md), [`split-fastq`](split-fastq.md) | `-m, --markers` |
    | [Reference FASTA (single-record)](#reference-fasta) | [`classify`](classify.md), [`split-fastq`](split-fastq.md) | `-r, --reference` |
    | [Reference FASTA (multi-record)](#reference-fasta) | [`match`](match.md) | `-r, --references` |
    | [Sample list — FASTA](#sample-list-tsv-classify-input-list) | [`classify`](classify.md) | `-l, --input-list` |
    | [Sample list — FASTQ](#sample-list-tsv-split-fastq-and-match-input-list) | [`split-fastq`](split-fastq.md), [`match`](match.md) | `-l, --input-list` |
    | [GFF3 annotation](#gff3-annotation-classify) | [`classify`](classify.md) | `--gff` / `--gff-files` / list column 3 |

---

## Training FASTA (`train`)

Standard FASTA. The **first whitespace-delimited token** of each header is used as the class label; everything after the first space is ignored for labeling (but preserved in the record header).

```text
>L4 sample_0001 collected_2021
ACTGATCGATCG...
>L2 sample_0002
ACTGATCGATCG...
```

Labels parsed above: `L4`, `L2`.

- The label token is taken verbatim — underscores, dots, and other characters are preserved (`L4.3.4.2`, `ST258`, `clade_A`, `L4_1_2` are all valid labels).
- Sequences are normalized by the FASTA reader on load.
- Records with an empty sequence are skipped with a warning.

!!! note "Requirements"

    - At least **2 distinct classes**.
    - At least a few sequences per class — **10+ per class recommended** for reliable training.

!!! tip "Query FASTA for `predict` and `classify`"

    [`predict`](predict.md) and [`classify`](classify.md) read the **same FASTA format** but do **not** parse a label from the header — the full header is carried through to the output. Multi-record FASTA is accepted. With `-i, --input` and `--input-files`, `classify` treats **each record as its own genome** (one output row per contig); with `-l, --input-list` every contig of a sample is aggregated under its sample name.

---

## Marker TSV (`classify` and `split-fastq`)

Tab-separated. Blank lines and lines beginning with `#` are ignored. A header row is **optional** and tolerated by both commands (see the admonition below).

### Columns

| Column | Content | Required |
|---|---|---|
| 1 | Genomic position (**1-based**) | ✅ |
| 2 | REF allele | ✅ |
| 3 | ALT allele | ✅ |
| 4 … first empty cell | Lineage hierarchy — one level per column | ✅ (at least one) |
| after the first empty cell | Annotation columns: **gene**, then **mutation** | Optional |

### Examples

=== "Lineage only"

    ```text
    #pos	ref	alt	level1	level2
    761155	C	T	L4	L4.9
    2155168	G	A	L2	L2.2
    ```

=== "With gene/mutation annotation"

    ```text
    #pos	ref	alt	level1	level2		gene	mutation
    761155	C	T	L4	L4.9		gyrA	Ser95Thr
    2155168	G	A	L2	L2.2		katG	Ser315Thr
    ```

    The blank column (two consecutive tabs) after `L4.9` / `L2.2` separates the lineage levels from the annotation columns.

=== "Small indel (classify only)"

    ```text
    #pos	ref	alt	level1	level2
    4247431	CC	C	L1	L1.1
    ```

    Indel markers are used by [`classify`](classify.md) (FASTA) and **silently skipped** by [`split-fastq`](split-fastq.md) (FASTQ).

!!! warning "Annotations require an empty separator cell"

    Lineage levels are read **left-to-right until the first empty cell**. Everything after that empty cell is treated as annotation. If you place `gene`/`mutation` **directly after** the last lineage level with no empty column in between, they are parsed as **extra lineage levels**, not annotations. To attach annotations, leave exactly one empty column (i.e. two consecutive tabs) after the deepest lineage level.

### Parsing rules

- **Positions are 1-based** relative to the [reference genome](#reference-fasta).
- Lineage levels are combined into a hierarchical, semicolon-joined path internally: `L4` + `L4.9` → `L4;L4.9`.
- Annotation columns: `classify` uses the **first** annotation cell as the gene and the **second** as the mutation; `split-fastq` retains all trailing annotation cells. When a GFF is supplied to `classify`, gene/amino-acid values derived from the GFF take precedence, falling back to these columns.
- A row is skipped if it has fewer than 4 columns, has a non-numeric position, or has an empty REF, ALT, or lineage path.
- **Allele-length limit:** an allele must be shorter than the k-mer size. For `classify`, the allele must also leave room for the flanks, so effectively `max(len(REF), len(ALT)) ≤ kmer_size − 2 × min_flank_bases` (default `min_flank_bases` = 10); markers that exceed this are skipped.

!!! info "Supported variant types differ by command"

    | Command | Input | SNPs | MNVs | Small indels |
    |---|---|:--:|:--:|:--:|
    | [`classify`](classify.md) | FASTA assemblies | ✅ | ✅ | ✅ |
    | [`split-fastq`](split-fastq.md) | FASTQ reads | ✅ | ✅ | ❌ |

    Indels are intentionally skipped in the FASTQ workflow: short reads across repetitive regions (e.g. PE/PPE genes in *M. tuberculosis*) produce unreliable k-mer matches. Indels are reliably resolved only through the full-length FASTA `classify` workflow.

!!! note "Header row handling"

    A header is not required, but both commands tolerate one:

    - `classify` skips any first-column value that is not a number (a `pos`/`position` header fails the numeric parse and is dropped).
    - `split-fastq` detects a header when the first data line (lower-cased) contains `pos`, `ref`, or `lineage`, and skips it.

    A header commented with a leading `#` works for both. See the [Marker format](marker_format.md) page for curation guidance and quality criteria.

---

## Sample List TSV (`classify --input-list`)

Tab-separated. One sample per row: a name, a FASTA path, and an optional GFF path.

```text
sample_A	/data/genomes/sample_A.fasta	/data/gff/sample_A.gff3
sample_B	/data/genomes/sample_B.fasta
sample_C	genomes/sample_C.fasta
```

| Column | Content | Required |
|---|---|---|
| 1 | Sample name (**must be unique**) | ✅ |
| 2 | Path to the query FASTA file | ✅ |
| 3 | Path to a per-sample [GFF3 annotation](#gff3-annotation-classify) | Optional |

- Sample names must be unique — a duplicate name aborts the run.
- Every FASTA path (and every GFF path, when given) must exist, or the run aborts with a "file does not exist" error.
- Blank first cells and rows beginning with `#` are skipped; a header row whose first cell is `sample`/`sample_name`/`genome` **and** whose second cell contains `fasta` or `path` is auto-detected and skipped.

!!! tip "Alternatives to `--input-list`"

    Instead of a list, `classify` also accepts a single genome with `-i, --input` (plus an optional `--gff`), or several files with `--input-files` (with `--gff-files` matched to each FASTA by filename stem). The third list column is simply the per-sample equivalent of `--gff`.

---

## Sample List TSV (`split-fastq` and `match --input-list`)

Tab-separated. One sample per row: a name followed by one or more FASTQ paths (plain or gzipped).

```text
sample_A	/data/reads/sample_A_R1.fastq.gz	/data/reads/sample_A_R2.fastq.gz
sample_B	/data/reads/sample_B.fastq.gz
```

| Column | Content | Required |
|---|---|---|
| 1 | Sample name | ✅ |
| 2+ | FASTQ path(s): **one** for single-end, **two** for paired-end (R1 then R2) | ✅ (at least one) |

- Paths are validated up front — a missing FASTQ aborts the run.
- Blank lines and lines beginning with `#` are skipped; rows with no FASTQ path are skipped with a warning.

!!! warning "`split-fastq` vs `match` treat the list differently"

    - **`split-fastq`** processes each row as an **independent sample** with its own genotype report (plus a combined summary). Sample names must be unique — a duplicate aborts the run.
    - **`match`** **pools every FASTQ file from all rows into a single combined query** and reports the single best-matching reference. The sample-name column is read but not used to separate results — list only the reads for the one query you want to match.

---

## Reference FASTA

The reference defines the coordinate system for marker positions (`classify`, `split-fastq`) or the pool of candidate genomes to match against (`match`).

=== "classify / split-fastq — single-record"

    `-r, --reference` — a **single-record** FASTA containing the genome the marker positions are numbered against. The sequence is normalized on load. A FASTA with **more than one record is rejected** with:

    ```text
    Reference FASTA '<path>' contains multiple records; provide a single-record FASTA.
    ```

=== "match — multi-record"

    `-r, --references` — a **multi-record** (multi-FASTA) file where **each record is one candidate reference genome**. `match` scores the query reads against every record and reports the best containment match.

---

## GFF3 Annotation (`classify`)

Standard 9-column, tab-delimited GFF3. Used only by [`classify`](classify.md) to translate variants into amino-acid changes. `split-fastq` and `match` do not consume GFF.

```text
##gff-version 3
NC_000962.3	RefSeq	CDS	7302	9818	.	+	0	ID=cds-GyrA;gene=gyrA;locus_tag=Rv0006
NC_000962.3	RefSeq	CDS	2153889	2156111	.	-	0	ID=cds-KatG;gene=katG;locus_tag=Rv1908c
```

- **Only `CDS` features are used**; all other feature types are ignored.
- **Gene name** is taken from the attributes column: `gene=` always wins when present; otherwise the **first** of `locus_tag=`, `Name=` or `ID=` **in the order they appear on the line** is used:

    1. `gene=…`
    2. `locus_tag=…`
    3. `Name=…`
    4. `ID=…`

    A CDS with none of these is labeled `Unknown`.

- **Strand** (`+` / `-`) sets the translation direction; any other value is treated as unknown and not translated.
- **Coordinates are 1-based inclusive** (GFF3 standard) and converted internally to 0-based for interval-tree lookups.
- SNPs and MNVs overlapping a CDS are translated to 3-letter amino-acid changes; indel-length differences are reported as a frameshift at the affected codon.
- Comment lines (`#…`) and blank lines are ignored.

!!! tip "Three ways to supply a GFF to `classify`"

    | Input mode | GFF flag | Matching |
    |---|---|---|
    | `-l, --input-list` | list column 3 | per sample (row) |
    | `-i, --input` | `--gff <file>` | applied to the single genome |
    | `--input-files` | `--gff-files <files…>` | matched to each FASTA by filename stem; a **single** GFF is applied to every FASTA |

---

## See also

- [Marker format](marker_format.md) — curating marker sets, hierarchical lineages, and quality criteria.
- [`pathotypr train`](train.md) — building a model from a labeled training FASTA.
- [`pathotypr predict`](predict.md) — applying a model to query FASTA files.
- [`pathotypr classify`](classify.md) — marker genotyping on assemblies (marker TSV + single-record reference + optional GFF).
- [`pathotypr split-fastq`](split-fastq.md) — marker genotyping directly from FASTQ reads.
- [`pathotypr match`](match.md) — best-reference search against a multi-record reference FASTA.
- [Installation](installation.md) — getting the `pathotypr` binary.
