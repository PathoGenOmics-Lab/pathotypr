# Output files

Every `pathotypr` subcommand writes tab-separated (TSV) result files, and most can additionally emit an Excel (`.xlsx`) copy. This page is the authoritative reference for **what each command writes, where, and what every column means**. For flags, examples, and algorithm background, follow the per-command links.

!!! abstract "Conventions used throughout"
    - **TSV**: UTF-8, tab-delimited, exactly one header row, one record per line. Cells with no value are written **empty** (nothing between the tabs).
    - **`--excel`**: when supported and enabled, the same table is written to a `.xlsx` alongside the TSV. The Excel name is derived from the TSV path by replacing a trailing `.tsv` with `.xlsx`; if the name does not end in `.tsv`, `.xlsx` is appended. Excel sheets add an autofilter, auto-sized columns, and (on some columns) conditional-formatting colours.
    - **Coordinates** are **1-based** unless stated otherwise.
    - **No output over an input**: a command whose output would be written over one of its own input files (for example `classify -l run.tsv -o run`, whose detail table is `run.tsv`) stops before writing anything and says which two files collide.
    - Numeric precision is fixed per column (e.g. `0.9800`) — see each table.

!!! note "Global flags do not change file contents"
    `-v, --verbose` (repeatable: `-v` = debug, `-vv` = trace), `-h, --help` are available on every subcommand (`-V, --version` only on the top-level `pathotypr` command). Verbosity only affects **terminal logging** (progress bars, resolved output paths, accuracy figures) — it never alters the files described below.

## At a glance

| Command | Always writes | Also writes |
|---|---|---|
| [`train`](train.md) | model bundle, `.importance.tsv`, `.importance.coords.tsv` | — (no `--excel`; accuracy is logged, not filed) |
| [`predict`](predict.md) | predictions TSV | `.xlsx` with `--excel` |
| [`classify`](classify.md) | detailed TSV, `_summary.tsv` | `.xlsx` with `--excel`; `_masked.fasta` with `--output-masked-fasta` |
| [`split-fastq`](split-fastq.md) | per-sample `_mutations.tsv` and `_undetermined.tsv`, `_summary.tsv` | `.xlsx` with `--excel` |
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
    Single-split (or cross-validated) accuracy and out-of-bag accuracy are **logged to the terminal** during the run, not saved. Per-fold accuracies are printed at the default verbosity when `--cv-folds` is used; `-v` only adds debug-level detail. `train` has no `--excel` option.

### `*.importance.tsv` columns

Ranked by how often each bucket was chosen as a split across all 100 trees (highest first), limited to the top 500.

| Column | Description |
|---|---|
| `rank` | 1-based rank of the bucket by split count. |
| `bucket` | Feature-hash bucket index (`0 … 1,048,575`). |
| `split_count` | Number of times this bucket was used as a split across the ensemble. |
| `importance_pct` | `split_count` as a percentage of the total split count over all buckets (4 decimals). |
| `kmers` | Comma-separated training k-mers that fall into this bucket. A bucket is the last ten bases of the canonical k-mer, so these are every k-mer sharing them, and the one the splits used is among them without being singled out: read them as candidates. |

### `*.importance.coords.tsv` columns

One row per occurrence of a discriminant k-mer, ordered by bucket rank, then sequence, then position.

| Column | Description |
|---|---|
| `rank` | Importance rank of the bucket this k-mer maps to. |
| `bucket` | Feature-hash bucket index. |
| `split_count` | Split count for the bucket. |
| `importance_pct` | Percentage of total splits attributable to the bucket (4 decimals). |
| `kmer` | The bases found at `position`, read on the sequence's own strand. That can be the reverse complement of the k-mer listed for the bucket in `*.importance.tsv`, which is the canonical one. |
| `sequence` | Full FASTA header of the training sequence containing the k-mer. |
| `lineage` | Class label of that sequence. |
| `position` | 1-based start of those bases in the sequence. |

Excerpts of a real run on four genomes of two lineages:

```text
rank	bucket	split_count	importance_pct	kmers
1	603290	2	2.1739	CAGACCACCGGGCATCAGCGG,AACCGGATTGAGCATCAGCGG,AGTCTGCCACCGCATCAGCGG,…
```

```text
rank	bucket	split_count	importance_pct	kmer	sequence	lineage	position
1	603290	2	2.1739	CCGCTGATGCCTACAAAACTG	L1 train_L1_0	L1	8020
1	603290	2	2.1739	CAGACCACCGGGCATCAGCGG	L1 train_L1_0	L1	242297
```

Every k-mer of bucket 603290 ends in `GCATCAGCGG`: those are the ten bases the
bucket is made of. The first coordinate row reads `CCGCTGATGCCTACAAAACTG` at
position 8020, the reverse complement of `CAGTTTTGTAGGCATCAGCGG`, which ends in
the same ten bases.

See also: [Training pipeline](algorithms/training.md), [Feature hashing](algorithms/feature-hashing.md), [Input formats](input-formats.md).

---

## `predict`

Writes one predictions table to `--output`, streamed as sequences are classified.

| File | When | Content |
|---|---|---|
| `<output>` | always | One row per input sequence. |
| `<output → .xlsx>` | with `--excel` | Same columns, with conditional-formatting colours (see below). |

**Header:** `Header  Predicted_Lineage  Confidence  Confidence_Margin  Other_Votes  Kmers_Used`

| Column | Description |
|---|---|
| `Header` | Sequence header from the input FASTA. |
| `Predicted_Lineage` | Most-voted class label across the ensemble, or `Undetermined` when the record yielded no k-mer; `Confidence` and `Confidence_Margin` are then empty and `Kmers_Used` is zero. |
| `Confidence` | Fraction of trees voting for the winner, `0`–`1` (4 decimals). |
| `Confidence_Margin` | Winner votes minus runner-up votes, over total trees, `0`–`1` (4 decimals). |
| `Other_Votes` | Up to the **top 3** runner-up labels as `label:fraction` (fraction to 2 decimals), comma-separated; empty when there are no other votes. |
| `Kmers_Used` | k-mers the record yielded. Zero means the model said nothing about it, and `Predicted_Lineage` is then `Undetermined` with the confidence columns empty. |

```text
Header	Predicted_Lineage	Confidence	Confidence_Margin	Other_Votes	Kmers_Used
sample_A	L2.2.1	0.9400	0.9000	L2.2:0.04,L4:0.02	4398877
sample_B	L4.3.4.2	0.7500	0.5500	L4.3:0.20,L4:0.05	4411512
masked_consensus_02	Undetermined				0
```

`Confidence` and the fractions in `Other_Votes` are shares of the same 100
trees, so they add up to at most 1, and `Confidence_Margin` is `Confidence`
minus the first fraction in `Other_Votes`. The last row is from a real run on an
N-masked sequence: no window of A, C, G and T long enough for the model's k, so
no k-mer, no vote, and the row says so instead of carrying a lineage.

!!! tip "Excel colour coding"
    In the `.xlsx`, **Confidence** is shaded red `< 0.5`, yellow `0.5–0.9`, green `≥ 0.9`; **Confidence_Margin** is shaded red `< 0.1`, yellow `0.1–0.5`, green `≥ 0.5` — so low-confidence calls stand out at a glance.

!!! note "Every input record gets a row"
    A record that yields no k-mer, because it is shorter than the model's k or N-masked, is reported as `Undetermined` with `Kmers_Used` zero, and the run ends by saying how many came out that way. Only a file with no records at all is an error.

See also: [Prediction](algorithms/prediction.md), [`train`](train.md).

---

## `classify`

Marker-based variant calling on assemblies. Writes a detailed per-marker table, a per-genome summary and a table of the markers that could not be read, plus optional Excel and masked-FASTA outputs. The `--output-prefix` (`-o`) drives all names.

| File | When | Content |
|---|---|---|
| `<prefix>.tsv` | always | One row per marker match per genome (the **detailed** table). |
| `<prefix>_summary.tsv` | always | One row per genome with its major-lineage call and its count of undetermined markers. |
| `<prefix>_undetermined.tsv` | always | One row per marker per genome that could be read neither as present nor as absent. Only the header when every marker was read. |
| `<prefix>.xlsx`, `<prefix>_summary.xlsx`, `<prefix>_undetermined.xlsx` | with `--excel` | Excel copies of the three tables. |
| `<input-stem>_masked.fasta` | with `--output-masked-fasta` | One file **per input assembly**, marker sites replaced by `N`, written to the prefix's parent directory. Only produced for reference-coordinate sequences (records whose length matches the reference); other inputs are skipped with the reason logged. Two inputs sharing a file stem do not overwrite each other — the second becomes `<stem>_2_masked.fasta`. See [classify](classify.md#output). |

!!! note "Prefix handling"
    If `--output-prefix` already ends in `.tsv` (e.g. `-o run.tsv`), the detailed file is used verbatim (`run.tsv`) and the summary and undetermined tables drop the trailing `.tsv` before appending (`run_summary.tsv`, `run_undetermined.tsv`). Otherwise `.tsv`, `_summary.tsv` and `_undetermined.tsv` are appended to the prefix.

### Detailed table — `<prefix>.tsv`

**Header:** `genome  k-mer  k-merPOS  SNPgenome  SNPreference  REF  ALT  lineage  Gene  Gene_Start  Gene_End  AA_Pos  AA_Change  alleles_found  evidence`

| Column | Description |
|---|---|
| `genome` | With `-i`, the record's whole FASTA header line; with `--input-files`, that line prefixed with `[filename]`; with `-l`, the sample name from the list. |
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
| `AA_Change` | Amino-acid change, translated from the sample's own codon when a GFF is given, otherwise the marker's mutation field; empty if neither. The sample's codon is used rather than the reference with this marker's allele dropped into it, so a second change in the same codon does not produce an amino acid nobody carries. |
| `alleles_found` | `ALT` when only the marker's alternative k-mer was found; `ALT,REF` when both turned up, meaning two copies of the region disagree or the alternative k-mer also occurs elsewhere, and an anchor then read this locus as carrying the variant. Both k-mers turning up is never on its own reported as the marker being present. |
| `evidence` | `kmer` when the whole k-mer matched, or `left_anchor` / `right_anchor` / `both_anchors` when the allele was read from the flanks because a nearby change hid the k-mer. An anchor read establishes no genome coordinate, so `k-mer`, `k-merPOS` and `SNPgenome` are empty on those rows. |

!!! note "Genomes with no marker hits still appear"
    A genome that matches no markers gets a single row containing only its name, with all remaining columns empty — so every input is accounted for.

### Summary table — `<prefix>_summary.tsv`

**Header:** `genome  lineage:count  major_lineage  undetermined  undetermined:count`

| Column | Description |
|---|---|
| `genome` | Genome/sample name. |
| `lineage:count` | Space-separated `lineage:count` pairs at **every hierarchical level** (each ancestor path prefix is counted), sorted by count descending, then by name. With a resistance panel, markers found per drug instead: a catalogue path is not a hierarchy, and its grade has spaces in it. |
| `major_lineage` | Final call. For a lineage panel: with `--nested-classification` the deepest fully-supported path wins, otherwise the most abundant lineage. For a WHO resistance catalogue it is the resistance verdict instead, for example `Resistant: INH,RIF`. See [classify](classify.md#resistance-panels). |
| `undetermined` | Number of markers that could be read neither way in this genome. |
| `undetermined:count` | Those markers as `label:count` pairs under the **first level** of their path only (the drug for a resistance panel, the main lineage for a lineage panel), sorted like `lineage:count`. Empty when every marker was read. |

### Undetermined table (`<prefix>_undetermined.tsv`)

Every marker that neither the k-mers nor the anchors could read: a change on
both sides of the variant, a region absent from the assembly, a contig break or
a repeated anchor hid the position, so the marker is neither present nor ruled
out. Markers sharing a position and reference allele are read together, so
carrying one of them does not leave the others undetermined. See
[classify](classify.md#present-absent-and-undetermined) for how to read these
rows in a resistance panel.

**Header:** `genome  SNPreference  REF  ALT  lineage  Gene  mutation  reason`

| Column | Description |
|---|---|
| `genome` | Genome/sample name. |
| `SNPreference` | 1-based marker position in the reference genome. |
| `REF`, `ALT` | Alleles of the marker. |
| `lineage` | Lineage (or drug) path of the marker. |
| `Gene`, `mutation` | The marker's own gene and mutation annotation columns; empty when the marker TSV has none, which includes WHO-shaped panels, where the gene and mutation are part of the path. |
| `reason` | `not_found`: nothing of the marker occurs and no anchor could read it. `anchors_disagree`: the two anchors read different alleles. `ambiguous_both_alleles`: both k-mers turned up and no anchor could say which allele the locus carries. `alt_kmer_in_reference`: the reference itself carries this marker's ALT k-mer, so finding it means nothing, and no anchor could tell the real locus from the copy. `ref_kmer_repeated`: the reference carries the window more than once and the genome did not show every copy. `no_reference_kmer`: the marker has no usable reference k-mer (a reference allele outside A, C, G and T) and the anchors could not read it either. |

Every example below is an excerpt of a real run.

=== "Detailed, lineage panel"

    H37Rv against the lineage panel, with `-i`: the `genome` column is the
    record's whole header line. No GFF was given, so the gene and amino-acid
    columns are empty.

    ```text
    genome	k-mer	k-merPOS	SNPgenome	SNPreference	REF	ALT	lineage	Gene	Gene_Start	Gene_End	AA_Pos	AA_Change	alleles_found	evidence
    NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome	CGCGTCGGTGAGTATCGTGCACGCGTGTGCG	15102	15117	15117	G	C	L4						ALT	kmer
    NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome	CGTTCGTAGCCGCAGGACGGCGGTGTCCACG	311598	311613	311613	T	G	L4;L4.9						ALT	kmer
    ```

=== "Detailed, resistance panel with a GFF"

    Two samples from an `--input-list`, against the WHO-shaped resistance panel
    and the H37Rv GFF. `double_435` carries both base changes of one codon, so
    three catalogue rows fire: the two-base one by its k-mer, and the two
    single-base ones by an anchor each, which is why their `k-mer`, `k-merPOS`
    and `SNPgenome` cells are empty. All three report the amino acid the
    sample's own codon reads. `h37rv` matches nothing, and gets the row with its
    name only.

    ```text
    genome	k-mer	k-merPOS	SNPgenome	SNPreference	REF	ALT	lineage	Gene	Gene_Start	Gene_End	AA_Pos	AA_Change	alleles_found	evidence
    double_435				761109	G	T	RIF;RIF-R;rpoB_D435Y;1) Assoc w R;rpoB;D435Y	rpoB	759807	763325	435	Asp435Phe(GAC>TTC)	ALT	left_anchor
    double_435	TGAGCCAATTCATGTTCCAGAACAACCCGCT	761095	761109	761109	GA	TT	RIF;RIF-R;rpoB_D435F;1) Assoc w R;rpoB;D435F	rpoB	759807	763325	435	Asp435Phe(GAC>TTC)	ALT	kmer
    double_435				761110	A	T	RIF;RIF-R;rpoB_D435V;1) Assoc w R;rpoB;D435V	rpoB	759807	763325	435	Asp435Phe(GAC>TTC)	ALT	right_anchor
    h37rv														
    ```

=== "Summary"

    ```text
    genome	lineage:count	major_lineage	undetermined	undetermined:count
    NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome	L4:62 L4;L4.9:6	L4	0	
    ```

    By default the most abundant level wins, and every marker also counts
    towards its ancestors, so H37Rv is called `L4` (62 markers) rather than
    `L4.9` (6). With `--nested-classification` the same genome is called
    `L4;L4.9`, the deepest level with support.

=== "Summary, resistance panel"

    ```text
    genome	lineage:count	major_lineage	undetermined	undetermined:count
    double_435	RIF:3	Resistant: RIF	23	RIF:23
    h37rv		No resistance markers found	0	
    ```

    With a resistance panel `lineage:count` counts the markers found per drug.
    The 23 undetermined markers of `double_435` are indels next to the changed
    codon: indels have no anchors, and the change breaks both of their k-mers.

=== "Undetermined"

    ```text
    genome	SNPreference	REF	ALT	lineage	Gene	mutation	reason
    double_435	761095	T	TGAG	RIF;RIF-R;rpoB_Leu430_Ser431insArg;2) Assoc w R - Interim;rpoB;Leu430_Ser431insArg			not_found
    double_435	761107	TGGACCA	T	RIF;RIF-R;rpoB_Asp435_Gln436del;2) Assoc w R - Interim;rpoB;Asp435_Gln436del			ref_kmer_repeated
    ```

    With a WHO-shaped panel the `Gene` and `mutation` cells are empty: the gene
    and mutation are levels five and six of the path, not annotation columns.

See also: [Assembly classification](algorithms/assembly-classification.md), [Marker format](marker_format.md).

---

## `split-fastq`

Alignment-free genotyping from FASTQ reads. Writes one detailed table **per sample** plus a combined summary. The `--output-prefix` (`-o`, default `split`) drives all names.

| File | When | Content |
|---|---|---|
| `<prefix>_<sample>_mutations.tsv` | always | Per-marker calls for one sample (one file per sample). |
| `<prefix>_<sample>_undetermined.tsv` | always | Markers the reads could answer neither way, one row each, with the `reason`. Only the header when every marker was read. See [split-fastq](split-fastq.md#per-sample-undetermined-tsv). |
| `<prefix>_summary.tsv` | always | One row per sample with its major-lineage call and its count of undetermined markers. |
| `<prefix>_<sample>_mutations.xlsx`, `<prefix>_<sample>_undetermined.xlsx`, `<prefix>_summary.xlsx` | with `--excel` | Excel copies of the three tables. |

!!! note "Sample names in filenames are sanitised"
    Characters outside `A–Z a–z 0–9 _ - .` become `_`, and leading dots / surrounding underscores are trimmed. If two samples sanitise to the same name, the second gets a `_2` suffix, the third `_3`, and so on. The `genome` column in the summary keeps the **original** sample name.

### Detailed table — `<prefix>_<sample>_mutations.tsv`

**Header:** `pos  ref_allele  alt_allele  ref_count  alt_count  alt_fraction  lineage_path  evidence  extra_annotations...`

Only markers meeting `--min-depth` coverage **and** `--min-alt-percent` are written.

| Column | Description |
|---|---|
| `pos` | 1-based marker position. |
| `ref_allele` | Reference allele. |
| `alt_allele` | Alternate allele. |
| `ref_count` | Reads reading the reference allele. |
| `alt_count` | Reads reading the alternate allele. |
| `alt_fraction` | ALT as a percentage of `ref_count + alt_count` (2 decimals, `0`–`100`). |
| `lineage_path` | `;`-joined lineage path for the marker. |
| `evidence` | `kmer` when the marker's own k-mer matched, or `left_anchor` / `right_anchor` / `both_anchors` when its locus was read by counting from one side because a nearby change hid both k-mers. The counts on that row are then the anchor's. |
| `extra_annotations…` | Zero or more **additional tab-separated columns** copied verbatim from the marker's annotation fields — present only when the marker defines them. |

!!! tip "Excel colour coding"
    In the detailed `.xlsx`, the annotation columns collapse to a single `annotations` header, and **alt_fraction** is shaded red `< 50`, yellow `50–95`, green `≥ 95`.

### Undetermined table (`<prefix>_<sample>_undetermined.tsv`)

**Header:** `pos  ref_allele  alt_allele  ref_count  alt_count  lineage_path  reason`

One row per marker the reads could answer neither way, with `reason` one of
`not_found`, `low_depth`, `no_unique_anchor`, `anchors_disagree`,
`alt_kmer_in_reference` or `ref_kmer_repeated`. A marker
read through a neighbour's window, or read from one of its own anchors, is not
listed: it was answered. See
[split-fastq](split-fastq.md#per-sample-undetermined-tsv).

### Summary table — `<prefix>_summary.tsv`

**Header:** `genome  lineage:count  major_lineage  undetermined  undetermined:count`

| Column | Description |
|---|---|
| `genome` | Original sample name. |
| `lineage:count` | Space-separated `lineage:count` pairs across all hierarchical levels, sorted by count descending, then by name. With a resistance panel, markers found per drug instead. |
| `major_lineage` | Final lineage call (`--nested-classification` enables hierarchical path validation; otherwise most-abundant wins). With a WHO drug-resistance panel, the drugs with a grade 1 or 2 marker instead. |
| `undetermined` | How many markers the reads could answer neither way. |
| `undetermined:count` | Those markers as `label:count` pairs under the first level of their path. Empty when every marker was read. |

Excerpts of a real run of two samples against the resistance panel. Both carry
the same rifampicin variant; in the second, another change seven bases away
hides both of its k-mers, and the marker is read from its left anchor instead.

=== "Mutations"

    ```text
    pos	ref_allele	alt_allele	ref_count	alt_count	alt_fraction	lineage_path	evidence	extra_annotations...
    761155	C	T	0	119	100.00	RIF;RIF-R;rpoB_S450L;1) Assoc w R;rpoB;S450L	kmer
    ```

    ```text
    pos	ref_allele	alt_allele	ref_count	alt_count	alt_fraction	lineage_path	evidence	extra_annotations...
    761155	C	T	0	112	100.00	RIF;RIF-R;rpoB_S450L;1) Assoc w R;rpoB;S450L	left_anchor
    ```

    The header always ends in `extra_annotations...`; the cells only appear
    when the panel has annotation columns after an empty separator, which the
    WHO-shaped panel does not.

=== "Undetermined"

    ```text
    pos	ref_allele	alt_allele	ref_count	alt_count	lineage_path	reason
    759545	TA	CG	0	0	RIF;RIF-R;rpoB_c.-261A>G;5) Not assoc w R;rpoB;c.-261A>G	not_found
    759546	A	G	0	0	RIF;RIF-R;rpoB_c.-261A>G;5) Not assoc w R;rpoB;c.-261A>G	not_found
    ```

    These reads cover a single region, so most of the panel is listed here as
    `not_found`: nothing read it. A whole-genome run lists only the markers
    something actually hid.

=== "Summary"

    ```text
    genome	lineage:count	major_lineage	undetermined	undetermined:count
    variant_and_neighbour	RIF:1	Resistant: RIF	87049	EMB:13401 OTHER:11171 RIF:9525 INH:9014 …
    variant_only	RIF:1	Resistant: RIF	87059	EMB:13401 OTHER:11171 RIF:9535 INH:9014 …
    ```

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

**Header:** `Query_Files  Best_Match_Reference  Shared_Kmer_Fraction  Query_Kmers`

| Column | Description |
|---|---|
| `Query_Files` | Comma-separated list of all input FASTQ paths that were pooled into the query. |
| `Best_Match_Reference` | Header of the top-scoring reference genome. |
| `Shared_Kmer_Fraction` | Weighted k-mer containment score, `0`–`1` (4 decimals). |
| `Query_Kmers` | The denominator: k-mer occurrences kept from the query after `--min-kmer-count`. |

```text
Query_Files	Best_Match_Reference	Shared_Kmer_Fraction	Query_Kmers
sample_reads.fastq	reference_A	0.8886	101500
```

A real run: 1,450 reads at 0.15x against two references 250 SNPs apart. The
score is the share of the query's 101,500 k-mer occurrences that `reference_A`
contains; `reference_B` scored 0.8812 and is not reported, since only the best
match is. A whole-genome run gives far more: about 35 million k-mer occurrences at 10x with 150 bp reads, and 106 million at 30x.

!!! note "No match, no data row"
    If no reference shares any k-mers with the query, the first reference is still reported with a `Shared_Kmer_Fraction` of `0.0000` — read a near-zero score as "no usable match". Only an empty reference FASTA yields a header-only report. Because the Excel writer needs the resolved output path, `--excel` produces a file **only when `-o` is also set** — with stdout output, no `.xlsx` is created.

See also: [Reference matching](algorithms/reference-matching.md).

---

## See also

- [Input formats](input-formats.md) — FASTA/FASTQ inputs and how labels are parsed.
- [Marker format](marker_format.md) — Marker TSV columns and lineage-path definitions used by `classify` and `split-fastq`.
- [Desktop GUI](gui.md) — The same outputs, produced from the graphical interface.
