# Output files

Every `pathotypr` subcommand writes tab-separated (TSV) result files, and most can additionally emit an Excel (`.xlsx`) copy. This page is the authoritative reference for **what each command writes, where, and what every column means**. For flags, examples, and algorithm background, follow the per-command links.

!!! abstract "Conventions used throughout"
    - **TSV**: UTF-8, tab-delimited, exactly one header row, one record per line. Cells with no value are written **empty** (nothing between the tabs).
    - **`--excel`**: when supported and enabled, the same table is written to a `.xlsx` alongside the TSV. The Excel name is derived from the TSV path by replacing a trailing `.tsv` with `.xlsx`; if the name does not end in `.tsv`, `.xlsx` is appended. Excel sheets add an autofilter, auto-sized columns, and (on some columns) conditional-formatting colours.
    - **Coordinates** are **1-based** unless stated otherwise.
    - Numeric precision is fixed per column (e.g. `0.9800`) — see each table.

!!! note "Global flags do not change file contents"
    `-v, --verbose` (repeatable: `-v` = debug, `-vv` = trace), `-h, --help`, and `-V, --version` are available on every subcommand. Verbosity only affects **terminal logging** (progress bars, resolved output paths, accuracy figures) — it never alters the files described below.

## At a glance

| Command | Always writes | Also writes |
|---|---|---|
| [`train`](train.md) | model bundle, `.importance.tsv`, `.importance.coords.tsv` | — (no `--excel`; accuracy is logged, not filed) |
| [`predict`](predict.md) | predictions TSV | `.xlsx` with `--excel` |
| [`classify`](classify.md) | detailed TSV, `_summary.tsv` | `.xlsx` with `--excel`; `_masked.fasta` with `--output-masked-fasta` |
| [`split-fastq`](split-fastq.md) | per-sample `_mutations.tsv`, `_summary.tsv` | `.xlsx` with `--excel` |
| [`match`](match.md) | best-match TSV (or stdout) | `.xlsx` with `--excel` **and** `-o` |

---

## `train`

Writes the trained model plus two feature-importance reports. Filenames are derived from `--output` by replacing its **final** extension. For `-o model.pathotypr.zst`:

| File | When | Content |
|---|---|---|
| `model.pathotypr.zst` | always | Model bundle (config + feature hasher + label encoder + 100 trees), bincode-serialized and zstd-compressed. This is the file [`predict`](predict.md) consumes. |
| `model.pathotypr.importance.tsv` | always | Top 500 feature-hash buckets ranked by ensemble split usage. |
| `model.pathotypr.importance.coords.tsv` | always | Per-occurrence genomic coordinates of the discriminant k-mers behind those buckets. |

!!! note "Accuracy is not written to a file"
    Single-split (or cross-validated) accuracy and out-of-bag accuracy are **logged to the terminal** during the run, not saved. Use `-v` for per-fold detail. `train` has no `--excel` option.

### `*.importance.tsv` columns

Ranked by how often each bucket was chosen as a split across all 100 trees (highest first), limited to the top 500.

| Column | Description |
|---|---|
| `rank` | 1-based rank of the bucket by split count. |
| `bucket` | Feature-hash bucket index (`0 … 1,048,575`). |
| `split_count` | Number of times this bucket was used as a split across the ensemble. |
| `importance_pct` | `split_count` as a percentage of the total split count over all buckets (4 decimals). |
| `kmers` | Comma-separated training k-mers that hash into this bucket. |

### `*.importance.coords.tsv` columns

One row per occurrence of a discriminant k-mer, ordered by bucket rank, then sequence, then position.

| Column | Description |
|---|---|
| `rank` | Importance rank of the bucket this k-mer maps to. |
| `bucket` | Feature-hash bucket index. |
| `split_count` | Split count for the bucket. |
| `importance_pct` | Percentage of total splits attributable to the bucket (4 decimals). |
| `kmer` | The specific k-mer occurrence. |
| `sequence` | Full FASTA header of the training sequence containing the k-mer. |
| `lineage` | Class label of that sequence. |
| `position` | Position of the k-mer within the sequence. |

See also: [Training pipeline](algorithms/training.md), [Feature hashing](algorithms/feature-hashing.md), [Input formats](input-formats.md).

---

## `predict`

Writes one predictions table to `--output`, streamed as sequences are classified.

| File | When | Content |
|---|---|---|
| `<output>` | always | One row per input sequence. |
| `<output → .xlsx>` | with `--excel` | Same columns, with conditional-formatting colours (see below). |

**Header:** `Header  Predicted_Lineage  Confidence  Confidence_Margin  Other_Votes`

| Column | Description |
|---|---|
| `Header` | Sequence header from the input FASTA. |
| `Predicted_Lineage` | Most-voted class label across the ensemble (`Unknown` if no tree produced a valid vote). |
| `Confidence` | Fraction of trees voting for the winner, `0`–`1` (4 decimals). |
| `Confidence_Margin` | Winner votes minus runner-up votes, over total trees, `0`–`1` (4 decimals). |
| `Other_Votes` | Up to the **top 3** runner-up labels as `label:fraction` (fraction to 2 decimals), comma-separated; empty when there are no other votes. |

```text
Header	Predicted_Lineage	Confidence	Confidence_Margin	Other_Votes
sample_A	L2.2.1	0.9800	0.9200	L2.2:0.04,L4:0.02
sample_B	L4.3.4.2	0.7500	0.3000	L4.3:0.20,L4:0.05
```

!!! tip "Excel colour coding"
    In the `.xlsx`, **Confidence** is shaded red `< 0.5`, yellow `0.5–0.9`, green `≥ 0.9`; **Confidence_Margin** is shaded red `< 0.1`, yellow `0.1–0.5`, green `≥ 0.5` — so low-confidence calls stand out at a glance.

!!! note "Sequences shorter than *k* are skipped"
    Input records shorter than the model's k-mer size produce no row (a warning is logged). If **no** records qualify, the command errors instead of writing an empty file.

See also: [Prediction](algorithms/prediction.md), [`train`](train.md).

---

## `classify`

Marker-based variant calling on assemblies. Writes a detailed per-marker table and a per-genome summary, plus optional Excel and masked-FASTA outputs. The `--output-prefix` (`-o`) drives all names.

| File | When | Content |
|---|---|---|
| `<prefix>.tsv` | always | One row per marker match per genome (the **detailed** table). |
| `<prefix>_summary.tsv` | always | One row per genome with its major-lineage call. |
| `<prefix>.xlsx`, `<prefix>_summary.xlsx` | with `--excel` | Excel copies of the two tables. |
| `<input-stem>_masked.fasta` | with `--output-masked-fasta` | One file **per input assembly**, marker sites replaced by `N`, written to the prefix's parent directory. Only produced for reference-coordinate sequences (records whose length matches the reference); other inputs are skipped with the reason logged. See [classify](classify.md#output). |

!!! note "Prefix handling"
    If `--output-prefix` already ends in `.tsv` (e.g. `-o run.tsv`), the detailed file is used verbatim (`run.tsv`) and the summary drops the trailing `.tsv` before appending (`run_summary.tsv`). Otherwise `.tsv` / `_summary.tsv` are appended to the prefix.

### Detailed table — `<prefix>.tsv`

**Header:** `genome  k-mer  k-merPOS  SNPgenome  SNPreference  REF  ALT  lineage  Gene  Gene_Start  Gene_End  AA_Pos  AA_Change`

| Column | Description |
|---|---|
| `genome` | Genome/sample name (from the input list or FASTA record ID). In `--input-files` mode it is prefixed with `[filename]`. |
| `k-mer` | The diagnostic marker k-mer matched in the query. |
| `k-merPOS` | 1-based start position of that k-mer in the query contig. |
| `SNPgenome` | 1-based position of the variant allele within the query contig. |
| `SNPreference` | 1-based marker position in the reference genome. |
| `REF` | Reference allele. |
| `ALT` | Alternate (marker) allele observed. |
| `lineage` | Lineage path assigned to the marker. |
| `Gene` | Gene ID — from the GFF if one was provided, otherwise the marker's own gene field; empty if neither. |
| `Gene_Start` | 1-based gene start from the GFF; empty if no annotation. |
| `Gene_End` | Gene end from the GFF; empty if no annotation. |
| `AA_Pos` | Amino-acid position of the change; empty if not resolvable. |
| `AA_Change` | Amino-acid change — from GFF translation, otherwise the marker's mutation field; empty if neither. |

!!! note "Genomes with no marker hits still appear"
    A genome that matches no markers gets a single row containing only its name, with all remaining columns empty — so every input is accounted for.

### Summary table — `<prefix>_summary.tsv`

**Header:** `genome  lineage:count  major_lineage`

| Column | Description |
|---|---|
| `genome` | Genome/sample name. |
| `lineage:count` | Space-separated `lineage:count` pairs at **every hierarchical level** (each ancestor path prefix is counted), sorted by count descending, then by name. |
| `major_lineage` | Final lineage call. With `--nested-classification`, the deepest fully-supported path wins; otherwise the most abundant lineage wins. |

=== "Detailed"

    ```text
    genome	k-mer	k-merPOS	SNPgenome	SNPreference	REF	ALT	lineage	Gene	Gene_Start	Gene_End	AA_Pos	AA_Change
    sampleA	ACGT…GTCA	1523	1533	745	C	T	L2;L2.2	rv0667	759807	763325	450	S450L
    sampleA	TTGA…CCAG	9021	9031	7362	G	A	L2	katG	2153889	2156111	315	
    sampleB												
    ```

=== "Summary"

    ```text
    genome	lineage:count	major_lineage
    sampleA	L2:2 L2.2:1	L2.2
    sampleB		
    ```

See also: [Assembly classification](algorithms/assembly-classification.md), [Marker format](marker_format.md).

---

## `split-fastq`

Alignment-free genotyping from FASTQ reads. Writes one detailed table **per sample** plus a combined summary. The `--output-prefix` (`-o`, default `split`) drives all names.

| File | When | Content |
|---|---|---|
| `<prefix>_<sample>_mutations.tsv` | always | Per-marker calls for one sample (one file per sample). |
| `<prefix>_summary.tsv` | always | One row per sample with its major-lineage call. |
| `<prefix>_<sample>_mutations.xlsx`, `<prefix>_summary.xlsx` | with `--excel` | Excel copies. |

!!! note "Sample names in filenames are sanitised"
    Characters outside `A–Z a–z 0–9 _ - .` become `_`, and leading dots / surrounding underscores are trimmed. If two samples sanitise to the same name, the second gets a `_2` suffix, the third `_3`, and so on. The `genome` column in the summary keeps the **original** sample name.

### Detailed table — `<prefix>_<sample>_mutations.tsv`

**Header:** `pos  ref_allele  alt_allele  ref_count  alt_count  alt_fraction  lineage_path  extra_annotations...`

Only markers meeting `--min-depth` coverage **and** `--min-alt-percent` are written.

| Column | Description |
|---|---|
| `pos` | 1-based marker position. |
| `ref_allele` | Reference allele. |
| `alt_allele` | Alternate allele. |
| `ref_count` | Reads supporting the REF k-mer. |
| `alt_count` | Reads supporting the ALT k-mer. |
| `alt_fraction` | ALT as a percentage of `ref_count + alt_count` (2 decimals, `0`–`100`). |
| `lineage_path` | `;`-joined lineage path for the marker. |
| `extra_annotations…` | Zero or more **additional tab-separated columns** copied verbatim from the marker's annotation fields — present only when the marker defines them. |

!!! tip "Excel colour coding"
    In the detailed `.xlsx`, the annotation columns collapse to a single `annotations` header, and **alt_fraction** is shaded red `< 50`, yellow `50–95`, green `≥ 95`.

### Summary table — `<prefix>_summary.tsv`

**Header:** `genome  lineage:count  major_lineage`

| Column | Description |
|---|---|
| `genome` | Original sample name. |
| `lineage:count` | Space-separated `lineage:count` pairs across all hierarchical levels, sorted by count descending, then by name. |
| `major_lineage` | Final lineage call (`--nested-classification` enables hierarchical path validation; otherwise most-abundant wins). |

!!! warning "No indels in FASTQ mode"
    `split-fastq` calls SNPs/MNVs only. Insertions and deletions are detected reliably via the assembly-based [`classify`](classify.md) workflow.

See also: [Marker genotyping](algorithms/marker-genotyping.md), [Marker format](marker_format.md).

---

## `match`

Finds the single best-matching reference for a set of FASTQ reads. All input reads (from `--input` or every file in `--input-list`) are pooled into **one** query, so the report has exactly one data row.

| File | When | Content |
|---|---|---|
| `<output>` | when `-o` is given (otherwise printed to **stdout**) | Best-match report. |
| `<output → .xlsx>` | with `--excel` **and** `-o` | Excel copy. |

**Header:** `Query_Files  Best_Match_Reference  Shared_Kmer_Fraction`

| Column | Description |
|---|---|
| `Query_Files` | Comma-separated list of all input FASTQ paths that were pooled into the query. |
| `Best_Match_Reference` | Header of the top-scoring reference genome. |
| `Shared_Kmer_Fraction` | Weighted k-mer containment score, `0`–`1` (4 decimals). |

```text
Query_Files	Best_Match_Reference	Shared_Kmer_Fraction
reads_R1.fastq.gz,reads_R2.fastq.gz	NC_000962.3 H37Rv	0.9873
```

!!! note "No match, no data row"
    If no reference shares any k-mers with the query, only the header row is written and a message is logged. Because the Excel writer needs the resolved output path, `--excel` produces a file **only when `-o` is also set** — with stdout output, no `.xlsx` is created.

See also: [Reference matching](algorithms/reference-matching.md).

---

## See also

- [Input formats](input-formats.md) — FASTA/FASTQ inputs and how labels are parsed.
- [Marker format](marker_format.md) — Marker TSV columns and lineage-path definitions used by `classify` and `split-fastq`.
- [Desktop GUI](gui.md) — The same outputs, produced from the graphical interface.
