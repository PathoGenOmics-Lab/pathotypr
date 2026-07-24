# Marker File Format

`pathotypr` classifies genomes and reads against a set of **lineage-diagnostic variants** supplied as a tab-separated (TSV) marker file. The same file format is consumed by two subcommands:

- [`pathotypr classify`](classify.md) — assembly / FASTA workflow (`get_positions` parser)
- [`pathotypr split-fastq`](split-fastq.md) — alignment-free FASTQ workflow (`build_markers` parser)

Both read the file column-for-column in the same order; the few behavioural differences (header detection, indel support) are called out explicitly below.

!!! note "One file, two engines"
    A single marker TSV works for both subcommands. Write it once in the portable form described here and it will classify assemblies **and** raw reads identically.

!!! warning "Panels for other organisms are untested"
    The format carries no assumption about species — the panel you write defines
    the organism. But pathotypr has only been validated on the *M. tuberculosis*
    complex, so a panel for another organism should be checked against a truth
    set you trust before you rely on its calls.

## Column layout

Fields are separated by a single **tab** (`\t`). Positions are **1-based** relative to the reference genome. Rows are read left-to-right:

| Column | Field | Required | Notes |
|---|---|---|---|
| 1 | **position** | ✅ | 1-based genomic coordinate; must parse as an integer |
| 2 | **REF** | ✅ | Reference allele — `A`, `C`, `G`, `T` (may be multi-base) |
| 3 | **ALT** | ✅ | Alternate allele — `A`, `C`, `G`, `T` (may be multi-base) |
| 4 … *m* | **lineage hierarchy** | ✅ (≥ 1) | One level per column, root → leaf (e.g. `L2` `L2.2` `L2.2.1`) |
| *m* + 1 | **(empty)** | conditional | Blank separator column — **required** if annotation columns follow |
| *m* + 2 | **gene** | optional | Gene name for reporting/annotation |
| *m* + 3 | **mutation** | optional | Amino-acid or nucleotide change for reporting |

A row must have **at least 4 columns** (position + REF + ALT + one lineage level) or it is skipped with a warning.

### How lineage and annotation columns are split

This is the part most people get wrong, so it is worth stating exactly.

Starting at **column 4**, `pathotypr` reads cells as lineage levels **until it hits the first empty cell**. That empty cell terminates the hierarchy. The non-empty cells *after* it are annotations — the first is treated as **gene**, the second as **mutation**.

```text
position  REF  ALT   L2   L2.2   L2.2.1  ⟨empty⟩  katG   Ser315Thr
   │       │    │    └──── lineage ────┘    │      └──── annotation ───┘
   1       2    3         4 … m             m+1          m+2, m+3
```

- In `classify`, the lineage levels are joined internally with `;` → `L2;L2.2;L2.2.1`.
- In `split-fastq`, the lineage levels are kept as an ordered list → `["L2", "L2.2", "L2.2.1"]`.

Both represent the **same** hierarchy, which is why one-level-per-column is the portable form.

!!! danger "The empty separator column is mandatory for annotations"
    Gene and mutation are only recognised when an **empty column** sits between the last lineage level and the first annotation. Without it, every trailing cell is swallowed into the lineage hierarchy and the gene/mutation stay unset. See the invalid example below.

!!! tip "Use one column per lineage level, not semicolons"
    Write the hierarchy across tab-separated columns (`L2` `⇥` `L2.2` `⇥` `L2.2.1`). `classify` also tolerates a single cell that already contains `L2;L2.2;L2.2.1`, but `split-fastq` treats that cell as **one** literal level — so the two workflows would disagree. Column-per-level is unambiguous everywhere.

## Header detection

A header row is **optional**. Both parsers tolerate one, but they detect it differently:

=== "`classify` (FASTA)"

    No dedicated header logic. Any row whose **first field is not an integer** is silently skipped, which naturally absorbs a `position ref alt …` header. Blank lines and lines beginning with `#` are also skipped.

=== "`split-fastq` (FASTQ)"

    Blank lines and `#` lines are skipped first. Then the **first remaining line** is treated as a header and skipped if — case-insensitively — it contains any of `pos`, `lineage`, or `ref`. Only that first line is checked.

!!! tip "Prefix the header with `#` for maximum portability"
    Starting the header line with `#` makes **both** parsers skip it unconditionally, regardless of the wording you choose:

    ```text
    #position	ref	alt	level1	level2	level3
    ```

## Variant types

| Variant | REF / ALT | `classify` (FASTA) | `split-fastq` (FASTQ) |
|---|---|---|---|
| SNP | 1 bp / 1 bp | ✅ | ✅ |
| MNV | *n* bp / *n* bp (equal length) | ✅ | ✅ |
| Indel | REF length ≠ ALT length | ✅ | ❌ skipped |

!!! warning "Indels are FASTA-only"
    `split-fastq` skips any marker whose REF and ALT differ in length (logged as `indels skipped`). Short reads across repetitive regions (PE/PPE families, IS elements in *M. tuberculosis*) produce unreliable k-mer matches for insertions/deletions. Detect indels with the assembly-based [`classify`](classify.md) workflow, which uses full-length context.

### Allele-length limits

A marker is only usable if its longest allele fits inside the k-mer window (`--kmer-size`, default `31`):

- **Both workflows:** `max(len(REF), len(ALT))` must be **< k** (< 31 by default).
- **`classify` only:** `--min-flank-bases` (default `10`) further requires
  `max(len(REF), len(ALT)) ≤ k − 2 × min_flank_bases` — i.e. **≤ 11 bp** at defaults. Longer alleles are skipped and reported in the run log.

!!! note "Use plain ACGT alleles"
    `pathotypr` 2-bit-encodes k-mers for speed. Ambiguity codes and other non-`ACGT` characters cannot be encoded: `classify` falls back to a slower string index, and in `split-fastq` the affected k-mers are simply never indexed (the marker becomes unmatchable). Keep alleles to `A`, `C`, `G`, `T`.

## Valid examples

=== "Lineage markers"

    One lineage level per column. Rows may carry different hierarchy depths.

    ```text
    #position	ref	alt	level1	level2	level3
    615938	A	G	L1
    1799921	C	G	L2
    801959	A	G	L2	L2.2
    2831482	C	T	L2	L2.2	L2.2.1
    ```

=== "With gene / mutation annotations"

    Note the **empty column** between the lineage levels and `gene` (shown here as the double gap — it is a real empty tab-separated cell).

    ```text
    #position	ref	alt	lineage		gene	mutation
    761155	C	T	RIF		rpoB	S450L
    2155168	G	A	INH		katG	Ser315Thr
    ```

    Parsed as: lineage `RIF`, gene `rpoB`, mutation `S450L` (and lineage `INH`, gene `katG`, mutation `Ser315Thr`).

=== "MNV"

    Equal-length REF/ALT — works in both workflows.

    ```text
    #position	ref	alt	lineage
    3820	AC	GT	L4
    ```

## Invalid / skipped rows

!!! failure "Annotations without the empty separator column"
    ```text
    #position	ref	alt	lineage	gene	mutation
    761155	C	T	RIF	rpoB	S450L
    ```
    With no empty cell after `RIF`, all three trailing cells become lineage levels → lineage `RIF;rpoB;S450L`, and **gene/mutation are never set**. Insert an empty column before `gene`.

Other rows that are dropped (with a log message):

| Row | Why it is skipped |
|---|---|
| `615938	A	G` | Fewer than 4 columns (no lineage level) |
| `foo	A	G	L1` | Position does not parse as an integer (treated as header/comment) |
| `615938		G	L1` | Empty REF allele |
| `615938	A		L1` | Empty ALT allele |
| `615938	A	G	` | No non-empty lineage level |
| `4247431	CC	C	L1` | Indel — skipped by `split-fastq` (accepted by `classify`) |

!!! info "Duplicate positions and shared markers"
    Multiple rows may share a position (e.g. different ALTs, or a marker diagnostic for more than one lineage). Distinct REF/ALT pairs generate distinct k-mers and coexist. If two rows generate the **identical** k-mer, `classify` keeps the last one and emits a duplicate-marker warning, while `split-fastq` records both hits.

## Reference genome

Markers are interpreted against a **single-record** FASTA reference supplied with `-r, --reference`. A multi-record FASTA is rejected. The marker positions must be 1-based coordinates on that exact sequence.

```bash
pathotypr classify \
  -m markers.tsv \
  -r reference.fasta \
  -i sample.fasta \
  -o classify_run
```

!!! tip "Debug skipped markers with verbose logging"
    Every subcommand accepts the global flags `-v, --verbose` (repeatable: `-v` = debug, `-vv` = trace), `-h, --help`, and `-V, --version`. Run with `-v` (or `-vv`) to see per-row messages explaining exactly which markers were skipped and why — invalid position, empty allele, allele too long for the k-mer window, or indel skipped in FASTQ mode:

    ```bash
    pathotypr split-fastq -m markers.tsv -r reference.fasta \
      -i reads_R1.fq.gz -i reads_R2.fq.gz -o run -vv
    ```

## Validation checklist

Before publishing a marker set, sanity-check it:

- [ ] Tab-separated, UTF-8, Unix line endings.
- [ ] Column order is `position` → `REF` → `ALT` → lineage levels → (empty) → gene → mutation.
- [ ] Positions are 1-based and land on the intended base of your single-record reference.
- [ ] Every row has at least one lineage level (≥ 4 columns).
- [ ] Annotation rows include the **empty separator column** before `gene`.
- [ ] Allele lengths respect the k-mer / flank limits (≤ 11 bp at default `classify` settings).
- [ ] If you rely on `split-fastq`, all markers are SNPs or equal-length MNVs (no indels).
- [ ] A dry run with `-vv` reports no unexpected skips.

## See also

- [`pathotypr classify`](classify.md) — assembly/FASTA classification and its options (`--kmer-size`, `--min-flank-bases`, `--gff`).
- [`pathotypr split-fastq`](split-fastq.md) — alignment-free FASTQ genotyping (`--min-depth`, `--min-alt-percent`, `--kmer-size`).
- [Input Format Reference](input-formats.md) — all input file formats (training FASTA, marker TSV, sample lists) in one place.
- [Marker Genotyping algorithm](algorithms/marker-genotyping.md) — diagnostic k-mer construction, Bloom-filter pre-check, and parallel scanning.
- [Assembly Classification algorithm](algorithms/assembly-classification.md) — how FASTA markers are matched and lineage paths validated.
