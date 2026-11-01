# Input Format Reference

pathotypr reads a small set of plain-text formats. Every file is **UTF-8**; every tabular file is **tab-separated** (`\t`) with **no spaces around the delimiters**. No output may be named after an input: a command whose output would land on one of its input files stops before writing anything, and names both. This page is the authoritative reference for each one: which command consumes it, the exact column semantics, and the parsing rules the tools apply.

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

## Example dataset

The examples on this page are invented, and they make up one small dataset.
[Output files](output-files.md) runs every command on it, so each output shown
there can be traced back to the inputs shown here.

- `reference.fasta`: `toy_ref`, a random 4,000 bp sequence holding two invented genes, `geneR` (3001 to 3600, `+` strand) and `geneI` (3701 to 3898, `-` strand).
- `lineage_markers.tsv`: seven lineage markers, for L1, L1.1, L2, L2.2, L4, and two for L4.9, one of them a deletion.
- `resistance_markers.tsv`: four markers in the WHO shape, three in `geneR` for RIF and one in `geneI` for INH.
- `annotation.gff`: the two genes.
- `samples.tsv` and `genomes/`: three assemblies and the list that names them.
- `reads.tsv` and `reads/`: error-free 100 bp reads of those assemblies at 30x, from both strands, and their list.
- `training.fasta`: ten genomes of each of L1, L2 and L4.
- `references.fasta`: one genome per lineage, the candidates for `match`.
- `queries.fasta`: the three assemblies and `masked_1`, a record of 300 `N`, for `predict`.

Every genome carries 15 SNPs its whole lineage shares, and each training genome
three of its own. The three samples differ in this:

| Sample | Lineage | Resistance changes | Also |
|---|---|---|---|
| `sample_A` | L2.2 | `geneR` S50L | |
| `sample_B` | L4.9 | `geneI` S15T, `geneR` A70V | a second change at 2111, six bases from its L4.9 SNP at 2105 |
| `sample_C` | L4 | `geneR` c.-15C>T | it lost the 201 bases from 2050 to 2250, where both L4.9 markers sit |

---

## Training FASTA (`train`)

Standard FASTA. The **first whitespace-delimited token** of each header is used as the class label; everything after the first space is ignored for labeling (but preserved in the record header).

```text
>L1 toy_L1_01
AATAAACGGCCTTACAGATTATACAGCCGGACCTGGCTTAATTCGTCTTATCTTAACCTC
TTCCACGTTGAACGGGAGGGTTCATGCTCTAACATCAGGGATCCATGCATTGTTTTCCAA
...
>L2 toy_L2_01
AATAAACGGCCTTACAGATTATACAGCCGGACCTGGCTTAATTCGTCTTATCTTAACCTC
TTCCACGTTGACCGGGAGGGTTCATGCCCTAACATCAGGGATCCATGCATTGTTTTCCAA
...
```

The start of `training.fasta` in the [example dataset](#example-dataset). The
labels parsed are `L1` and `L2`; `toy_L1_01` stays in the header and plays no
part in the label.

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

Tab-separated, with LF or CRLF line endings (the deposited panels use CRLF; both are read the same). Blank lines and lines beginning with `#` are ignored. A header row is **optional** and tolerated by both commands (see the admonition below).

The same file format carries two kinds of panel, and the tools treat them differently:

- a **lineage panel**, where the columns after the alleles are lineage levels (`L4`, `L4.9`);
- a **resistance panel** in the shape of the WHO catalogue, where they are a drug, a phenotype, a marker name, a WHO grade, a gene and a mutation. See [Resistance panels](#resistance-panels-who-catalogue) below.

### Columns

| Column | Content | Required |
|---|---|---|
| 1 | Genomic position (**1-based**) | ✅ |
| 2 | REF allele | ✅ |
| 3 | ALT allele | ✅ |
| 4 … first empty cell | Lineage hierarchy — one level per column | ✅ (at least one) |
| after the first empty cell | Annotation columns: **gene**, then **mutation** | Optional |

### Examples

=== "Lineage panel"

    ```text
    position	ref_base	alt_base	lineage1	lineage2
    301	T	G	L1
    702	T	G	L1	L1.1
    1103	G	T	L2
    1504	C	T	L2	L2.2
    2105	G	C	L4	L4.9
    2207	CA	C	L4	L4.9
    2606	A	G	L4
    ```

    `lineage_markers.tsv` of the [example dataset](#example-dataset), with the
    header the deposited lineage panel uses. A sublineage marker lists its
    parent first, so `1504` is an L2.2 marker and counts for L2 as well.

=== "With gene/mutation annotation"

    ```text
    #pos	ref	alt	level1	level2		gene	mutation
    3149	C	T	RIF	RIF-R		geneR	S50L
    3855	C	G	INH	INH-R		geneI	S15T
    ```

    The blank column (two consecutive tabs) after `RIF-R` / `INH-R` separates the levels from the annotation columns: the path of the first marker is `RIF;RIF-R`, its gene `geneR` and its mutation `S50L`. A panel written this way is read as a **lineage** panel, with `RIF` and `INH` as its top levels: `classify` calls `sample_A` `RIF`, as it would call a lineage. To be answered as a resistance panel it has to use the WHO shape below.

=== "Small indel (classify only)"

    ```text
    2207	CA	C	L4	L4.9
    ```

    The deletion of the example panel: `CA` at 2207 and 2208 becomes `C`, so
    the `A` at 2208 is gone. Indel markers are used by [`classify`](classify.md)
    (FASTA) and skipped by [`split-fastq`](split-fastq.md) (FASTQ), which says
    how many it skipped in its log.

!!! warning "Annotations require an empty separator cell"

    Lineage levels are read **left-to-right until the first empty cell**. Everything after that empty cell is treated as annotation. If you place `gene`/`mutation` **directly after** the last lineage level with no empty column in between, they are parsed as **extra lineage levels**, not annotations. To attach annotations, leave exactly one empty column (i.e. two consecutive tabs) after the deepest lineage level.

### Parsing rules

- **Positions are 1-based** relative to the [reference genome](#reference-fasta).
- Lineage levels are combined into a hierarchical, semicolon-joined path internally: `L4` + `L4.9` → `L4;L4.9`.
- Annotation columns: `classify` uses the **first** annotation cell as the gene and the **second** as the mutation; `split-fastq` retains all trailing annotation cells. When a GFF is supplied to `classify`, gene/amino-acid values derived from the GFF take precedence, falling back to these columns.
- A row is skipped, with a warning that names it, if it has fewer than 4 columns, has a non-numeric position, or has an empty REF, ALT, or lineage path. The rest of the file is still read.
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

### Resistance panels (WHO catalogue)

A resistance panel has nine columns and **no empty separator cell**, so all six
columns after the alleles form the marker's path. This is
`resistance_markers.tsv` of the [example dataset](#example-dataset):

```text
#pos	ref	alt	drug	resistance	marker_name	grade	gene	mutation
2986	C	T	RIF	RIF-R	geneR_c.-15C>T	5) Not assoc w R	geneR	c.-15C>T
3149	C	T	RIF	RIF-R	geneR_S50L	1) Assoc w R	geneR	S50L
3209	C	T	RIF	RIF-R	geneR_A70V	3) Uncertain significance	geneR	A70V
3855	C	G	INH	INH-R	geneI_S15T	1) Assoc w R	geneI	S15T
```

The path of `geneR_S50L` is `RIF;RIF-R;geneR_S50L;1) Assoc w R;geneR;S50L`,
and it is read by position: the **first** level is the drug and the **fourth**
is the WHO grade. The fifth and sixth are the gene and mutation the app shows on
its resistance card; in the TSV outputs they stay inside the path, and the
`Gene` column of `classify` is filled from the GFF when one is given.

The deposited panel has the same shape. These three rows are taken from it:

```text
#pos	ref	alt	drug	resistance	marker_name	grade	gene	mutation
761155	C	T	RIF	RIF-R	rpoB_S450L	1) Assoc w R	rpoB	S450L
2155168	C	G	INH	INH-R	katG_S315T	1) Assoc w R	katG	S315T
759545	TA	CG	RIF	RIF-R	rpoB_c.-261A>G	5) Not assoc w R	rpoB	c.-261A>G
```

A panel is recognised as a resistance panel when at least half of its markers
have a path of four levels or more whose fourth level reads as a WHO grade, so a
lineage panel with one odd row is never mistaken for one. The grade decides what
a marker found in a genome means:

| Grade cell | Read as |
|---|---|
| `1) Assoc w R`, `2) Assoc w R - Interim` | associated with resistance: the drug is called resistant |
| `3) Uncertain significance` | uncertain: reported, never called resistant |
| `4) Not assoc w R - Interim`, `5) Not assoc w R` | not associated |
| anything else, such as `TBProfiler_literature` | left out of the resistance call; the marker is still found and listed |

The number before the parenthesis is what is read; without one, a cell
containing "not assoc", "uncertain" or "assoc" is read by those words. The
deposited panel has 402 rows graded `TBProfiler_literature`, which the call
leaves out. With a
resistance panel the summary's call is the list of drugs with a grade 1 or 2
marker, for example `Resistant: INH,RIF`, instead of a lineage; see
[classify](classify.md#resistance-panels).

!!! note "Several rows per catalogue entry"
    The deposited panel spells one catalogue mutation several ways: the single
    base change, and longer alleles covering the same codon with its
    neighbours. Each row is a marker of its own, and a sample can match more
    than one of them for a single change.

---

## Sample List TSV (`classify --input-list`)

Tab-separated. One sample per row: a name, a FASTA path, and an optional GFF path.

```text
sample_A	genomes/sample_A.fasta	annotation.gff
sample_B	genomes/sample_B.fasta	annotation.gff
sample_C	genomes/sample_C.fasta	annotation.gff
```

`samples.tsv` of the [example dataset](#example-dataset). A sample without an
annotation leaves the third column out, and any path may be absolute.

| Column | Content | Required |
|---|---|---|
| 1 | Sample name (**must be unique**) | ✅ |
| 2 | Path to the query FASTA file | ✅ |
| 3 | Path to a per-sample [GFF3 annotation](#gff3-annotation-classify) | Optional |

- Sample names must be unique — a duplicate name aborts the run.
- Every FASTA path (and every GFF path, when given) must exist, or the run aborts with a "file does not exist" error.
- A relative path is read from the directory pathotypr runs in, not from the list's own folder: the example list works from the folder that holds `genomes/`.
- Blank first cells and rows beginning with `#` are skipped; a header row whose first cell is `sample`/`sample_name`/`genome` **and** whose second cell contains `fasta` or `path` is auto-detected and skipped.

!!! tip "Alternatives to `--input-list`"

    Instead of a list, `classify` also accepts a single genome with `-i, --input` (plus an optional `--gff`), or several files with `--input-files` (with `--gff-files` matched to each FASTA by filename stem). The third list column is simply the per-sample equivalent of `--gff`.

---

## Sample List TSV (`split-fastq` and `match --input-list`)

Tab-separated. One sample per row: a name followed by one or more FASTQ paths (plain or gzipped).

```text
sample_A	reads/sample_A.fastq
sample_B	reads/sample_B.fastq
sample_C	reads/sample_C.fastq
```

`reads.tsv` of the [example dataset](#example-dataset): single-end reads, one
file per sample. A paired-end sample lists R1 and then R2 on its row, plain or
gzipped:

```text
sample_D	/data/reads/sample_D_R1.fastq.gz	/data/reads/sample_D_R2.fastq.gz
```

| Column | Content | Required |
|---|---|---|
| 1 | Sample name | ✅ |
| 2+ | FASTQ path(s): **one** for single-end, **two** for paired-end (R1 then R2) | ✅ (at least one) |

- Paths are validated up front — a missing FASTQ aborts the run.
- A relative path is read from the directory pathotypr runs in, not from the list's own folder.
- Blank lines and lines beginning with `#` are skipped; rows with no FASTQ path are skipped with a warning.

!!! warning "`split-fastq` vs `match` treat the list differently"

    - **`split-fastq`** processes each row as an **independent sample** with its own genotype report (plus a combined summary). Sample names must be unique — a duplicate aborts the run.
    - **`match`** **pools every FASTQ file from all rows into a single combined query** and reports the single best-matching reference. The sample-name column is read but not used to separate results — list only the reads for the one query you want to match.

---

## Reference FASTA

The reference defines the coordinate system for marker positions (`classify`, `split-fastq`) or the pool of candidate genomes to match against (`match`).

=== "classify / split-fastq — single-record"

    `-r, --reference`: a **single-record** FASTA containing the genome the marker positions are numbered against. The sequence is normalized on load. The reference of the [example dataset](#example-dataset):

    ```text
    >toy_ref invented example reference
    AATAAACGGCCTTACAGATTATACAGCCGGACCTGGCTTAATTCGTCTTATCTTAACCTC
    TTCCACGTTGAACGGGAGGGTTCATGCTCTAACATCAGGGATCCATGCATTGTTTTCCAA
    ...
    ```

    A FASTA with **more than one record is rejected** with:

    ```text
    Reference FASTA '<path>' contains multiple records; provide a single-record FASTA.
    ```

=== "match — multi-record"

    `-r, --references`: a **multi-record** (multi-FASTA) file where **each record is one candidate reference genome**. `match` scores the query reads against every record and reports the best containment match. The candidates of the [example dataset](#example-dataset):

    ```text
    >ref_L1 invented L1 representative
    AATAAACGGCCTTACAGATTATACAGCCGGACCTGGCTTAATTCGTCTTATCTTAACCTC
    ...
    >ref_L2 invented L2 representative
    AATAAACGGCCTTACAGATTATACAGCCGGACCTGGCTTAATTCGTCTTATCTTAACCTC
    ...
    >ref_L4 invented L4 representative
    AATAAACGGCCTTACAGATTATACAGCCGGACCTGGCTTAATTCGTCTTATCTTAACCTC
    ...
    ```

---

## GFF3 Annotation (`classify`)

Standard 9-column, tab-delimited GFF3. Used only by [`classify`](classify.md) to translate variants into amino-acid changes. `split-fastq` and `match` do not consume GFF.

```text
##gff-version 3
toy_ref	invented	CDS	3001	3600	.	+	0	ID=cds-geneR;gene=geneR;locus_tag=TOY_0001
toy_ref	invented	CDS	3701	3898	.	-	0	ID=cds-geneI;gene=geneI;locus_tag=TOY_0002
```

`annotation.gff` of the [example dataset](#example-dataset). The RefSeq
annotation of H37Rv has the same shape:

```text
##gff-version 3
NC_000962.3	RefSeq	CDS	7302	9818	.	+	0	ID=cds-GyrA;gene=gyrA;locus_tag=Rv0006
NC_000962.3	RefSeq	CDS	2153889	2156111	.	-	0	ID=cds-KatG;gene=katG;locus_tag=Rv1908c
```

- **Only `CDS` features are used**; all other feature types are ignored.
- **Gene name** is taken from the attributes column. `gene=` wins wherever it appears; without it, whichever of `locus_tag=`, `Name=` or `ID=` comes **first on the line** is used, so `ID=cds-geneR;Name=GeneR;locus_tag=TOY_0001` is labeled `cds-geneR` and `locus_tag=TOY_0002;ID=cds-geneI` is labeled `TOY_0002`. A CDS with none of these is labeled `Unknown`.

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
