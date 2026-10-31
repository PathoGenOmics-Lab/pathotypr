# `pathotypr classify`

**Assigns each input genome to a lineage by scanning it for a curated set of diagnostic marker k-mers.**

`pathotypr classify` performs marker-based genotyping on assembled genomes. It builds diagnostic k-mers from your marker definitions and a reference sequence, scans every query assembly for exact matches, optionally annotates each hit with gene and amino-acid information from a GFF, and reports both a detailed per-marker table and a per-genome lineage summary. Use it when you have one or many draft or complete assemblies (draft, multi-contig assemblies included) and a marker panel that discriminates the lineages you care about.

!!! note "Assemblies in, lineages out"
    `classify` works on nucleotide FASTA assemblies, not raw reads. If you are starting from sequencing reads, assemble them first (or use the read-based workflow) and pass the resulting FASTA here.

## Synopsis

```bash
pathotypr classify \
  -m <markers.tsv> \
  -r <reference.fasta> \
  (-i <sample.fasta> | -l <list.tsv> | --input-files <a.fasta> ...) \
  -o <prefix> \
  [OPTIONS]
```

Exactly one input source is required: a single FASTA (`-i`), a TSV list of samples (`-l`), or several FASTA files (`--input-files`). With `-i` and `--input-files` each FASTA **record** is treated as its own genome, so a multi-contig assembly yields one row per contig; `-l` aggregates every contig of a sample under its sample name.

## Inputs

`classify` needs three things, plus an optional annotation.

| What | Flag | Format |
|---|---|---|
| Marker panel | `-m, --markers` | [Marker TSV](marker_format.md), tab-separated |
| Reference genome | `-r, --reference` | **Single-record** FASTA that the marker positions are numbered against |
| Genomes to type | one of `-i` / `-l` / `--input-files` | FASTA, plain or gzipped |
| Annotation *(optional)* | `--gff` / `--gff-files` / list column 3 | [GFF3](input-formats.md#gff3-annotation-classify); adds gene and amino-acid columns |

### Choosing the input source

The three sources are **not** interchangeable, and the difference decides how your results are grouped:

| Source | One row per | Use when |
|---|---|---|
| `-i, --input` | FASTA **record** | A single genome, or you want per-contig results |
| `--input-files` | FASTA **record**, named `[file] record` | Several assemblies, per-contig results |
| `-l, --input-list` | **sample** | Draft assemblies: every contig of a sample is aggregated under one name |

!!! warning "Draft assemblies belong in `--input-list`"
    With `-i` or `--input-files`, a 200-contig draft produces **200 rows**, each
    typed on the handful of markers that happen to fall on that contig, and none
    of them representing the isolate. Use `-l` and give the sample a name, so its
    contigs are counted together.

### Requirements

The run **fails** unless all of these hold:

- [x] The reference is a **single record**. A multi-FASTA is rejected with `contains multiple records; provide a single-record FASTA`.
- [x] Marker positions are **1-based coordinates on that exact reference**. Nothing cross-checks this, so a panel built against a different assembly produces silently wrong calls.
- [x] `--kmer-size` between **1 and 31**.
- [x] `--min-flank-bases` below **half** of `--kmer-size`, or no marker k-mer can be built.
- [x] Sample names in `-l` are **unique**; a duplicate aborts the run.
- [x] Every path in the list exists, checked before any work starts.

Rows in the marker TSV are skipped, with a log message, when they have fewer than 4 columns, a non-numeric position, an empty REF/ALT/lineage, or an allele too long for the k-mer window (`max(len(REF), len(ALT)) ≤ kmer_size − 2 × min_flank_bases`, so **11 bp** at defaults).

!!! tip "Both strands are scanned"
    Contig orientation in a draft assembly is arbitrary, so each contig is
    scanned forward and reverse-complemented. A marker found on several contigs,
    or on both strands, still counts **once**.

## Options

| Option | Default | Required | Description |
|---|---|---|---|
| `-m, --markers <FILE>` | — | yes | Marker definition TSV (position, REF, ALT, lineage levels...). |
| `-r, --reference <FILE>` | — | yes | Single-record reference FASTA the marker positions refer to. |
| `-i, --input <FILE>` | — | one of `-i`/`-l`/`--input-files` | A single sample FASTA. Each FASTA record is treated as a separate genome, so a multi-contig assembly yields one summary row per contig. |
| `-l, --input-list <FILE>` | — | one of `-i`/`-l`/`--input-files` | TSV of samples: `sample_name\tfasta_path[\tgff_path]`. All contigs of a sample's FASTA are aggregated under its sample name. |
| `--input-files <FILE>...` | — | one of `-i`/`-l`/`--input-files` | Multiple sample FASTA files. Each record stays a separate genome, named `[filename] record_id`. |
| `--gff <FILE>` | — | no | GFF annotation for `--input` (adds gene + amino-acid change columns). Requires `--input`. |
| `--gff-files <FILE>...` | — | no | Multiple GFFs matched by filename, for `--input-files`. |
| `-o, --output-prefix <PREFIX>` | — | yes | Prefix for all output files. |
| `--kmer-size <N>` | `31` | no | Marker k-mer size (1–31; 31 is both the default and the maximum). |
| `-t, --threads <N>` | all cores | no | Number of CPU threads. |
| `--nested-classification` | `false` | no | Enable hierarchical (nested) lineage calling using multi-level marker columns. |
| `--min-flank-bases <N>` | `10` | no | Minimum flanking bases on each side of the allele in a marker k-mer. Must be below half of `--kmer-size`, or no k-mer can be built and the run errors. |
| `--output-masked-fasta` | `false` | no | Also write masked FASTA files with marker positions replaced by `N`. Only applies to reference-coordinate sequences (see below). |
| `--excel` | `false` | no | Also write Excel (`.xlsx`) files alongside the TSVs. |

### Global flags

Available on every subcommand:

| Flag | Description |
|---|---|
| `-v, --verbose` | Increase log verbosity. Repeatable: `-v` = debug, `-vv` = trace. |
| `-h, --help` | Print help. |

!!! warning "`--gff` is single-sample only"
    `--gff` applies to `--input` mode only. For `--input-files` batches use `--gff-files`, where each GFF is paired to a FASTA by filename. For `--input-list` mode, add the GFF path as an optional third column in the list TSV.

## How it works

1. The reference FASTA and marker TSV are loaded, and two k-mers are built per marker from the same window: the ALT allele flanked on both sides by reference context (at least `--min-flank-bases` bases per side, total length `--kmer-size`), and the REF allele in the same flanks. The first says the variant is there; the second says it is not. Each substitution marker also gets two **anchors**, the `k-1` bases on either side of the variant, used when neither k-mer matches.
2. Every contig of every query genome is scanned on **both strands** — the contig and its reverse complement — because contig orientation in a draft assembly is arbitrary. A marker found more than once (on several contigs, or on both strands) counts as one observation.
3. If a GFF is supplied for a sample, each match is annotated with the overlapping gene and the resulting amino-acid change.
4. Each marker ends up in one of four states (see [Present, absent and undetermined](#present-absent-and-undetermined)). Only markers whose ALT k-mer was found count towards the lineage call.
5. Per genome, the observed markers are tallied per lineage and a `major_lineage` is called. With `--nested-classification`, multi-level marker columns are used to validate the hierarchical lineage path; without it, the most abundant lineage is reported directly.
6. Results are streamed to a detailed per-marker TSV, a per-genome summary TSV and a TSV of the markers that could not be read (plus optional `.xlsx` and masked FASTA outputs).

## Examples

=== "Single genome + GFF"

    Classify one assembly and annotate hits with gene and amino-acid changes, also emitting Excel:

    ```bash
    pathotypr classify \
      -m markers.tsv \
      -r reference.fasta \
      -i sample.fasta \
      --gff sample.gff3 \
      -o classify_run \
      --excel
    ```

=== "Batch from a list"

    Classify many samples described in a TSV (`sample_name`, `fasta_path`, optional `gff_path`), using nested lineage calling:

    ```bash
    pathotypr classify \
      -m markers.tsv \
      -r reference.fasta \
      -l samples.tsv \
      -o classify_batch \
      --nested-classification
    ```

    Example `samples.tsv`:

    ```tsv
    isolate_001	assemblies/isolate_001.fasta	annotations/isolate_001.gff3
    isolate_002	assemblies/isolate_002.fasta
    isolate_003	assemblies/isolate_003.fasta	annotations/isolate_003.gff3
    ```

=== "Multiple FASTA files"

    Classify several assemblies at once and also write masked FASTA files, with paired GFFs matched by filename:

    ```bash
    pathotypr classify \
      -m markers.tsv \
      -r reference.fasta \
      --input-files sample1.fasta sample2.fasta sample3.fasta \
      --gff-files sample1.gff3 sample2.gff3 sample3.gff3 \
      -o batch_run \
      --output-masked-fasta \
      -t 8
    ```

!!! tip "Tune the k-mer size to your panel"
    `--kmer-size` (default `31`) and `--min-flank-bases` (default `10`) together control how much reference context surrounds each allele. Larger k-mers are more specific but less tolerant of nearby variation in the query; keep enough flanking bases on each side so markers remain uniquely anchored.

## Output

All files share the `-o/--output-prefix` value. Given `-o classify_run`:

| File | When | Content |
|---|---|---|
| `classify_run.tsv` | always | Detailed table: one row per marker match per genome. |
| `classify_run_summary.tsv` | always | Summary table: one row per genome with lineage counts, the major lineage and how many markers could not be read. |
| `classify_run_undetermined.tsv` | always | One row per marker per genome that could be read neither as present nor as absent. Header only when every marker was read. |
| `classify_run.xlsx`, `classify_run_summary.xlsx`, `classify_run_undetermined.xlsx` | `--excel` | Excel (`.xlsx`) versions of the three tables, written alongside the TSVs. |
| `<input>_masked.fasta` | `--output-masked-fasta` | One masked FASTA per input FASTA file (named after the input file stem), with marker positions replaced by `N`. Colliding stems are disambiguated rather than overwritten. |

!!! warning "Masking needs reference coordinates"
    Mask positions come from the marker TSV, which is defined against the
    reference. They can therefore only be applied to a sequence that lives in
    the **reference coordinate system** — a reference-based consensus genome, or
    an alignment whose records all have the reference length.

    The contigs of a de-novo assembly are unrelated to those coordinates, so
    masking one would replace arbitrary positions with `N`. When a record's
    length does not match the reference, `classify` logs the reason and writes
    no masked file for that input rather than emitting a misleading one. The
    detailed and summary tables are unaffected.

!!! note
    If the output prefix already ends in `.tsv`, the detailed file uses it as-is and the summary is derived from the trimmed base name. Genomes with no marker hits still appear in the detailed TSV as a single row with empty match columns, and in the summary with an empty `lineage:count` and `Unclassified` as the major lineage — so every input is accounted for.

### Present, absent and undetermined

A k-mer is an exact match, so any other change inside a marker's window,
however harmless (a synonymous substitution, a lineage polymorphism), stops its
ALT k-mer from matching. Looking for the REF k-mer as well is what tells that
case apart from a genome that simply lacks the variant:

| ALT k-mer | REF k-mer | State | Where it is reported |
|---|---|---|---|
| found | not found | **present** | detailed TSV, `alleles_found` = `ALT` |
| found | found | **both alleles** | detailed TSV, `alleles_found` = `ALT,REF` |
| not found | found | **absent**: the genome carries the reference allele | nowhere, as before |
| not found | not found | the **anchors** decide (below) | |

Markers at the same position with the same reference allele share their REF
k-mer, so reading any one allele there reads the position for all of them: a
genome carrying one of them is not reported undetermined for the others.

#### The anchors

When neither k-mer matches, `classify` falls back to the two flanks on their
own. Each substitution marker carries a **left anchor**, the `k-1` bases ending
where the variant starts, and a **right anchor**, the `k-1` bases beginning
where it ends. Neither spans the variant, so one change on one side leaves the
other anchor whole; the allele is then read off the genome by counting bases
from the anchor instead of being recognised as a whole k-mer.

That reads alleles no catalogue lists, which an index of expected k-mers
cannot do:

| What the anchors read | State | Where it is reported |
|---|---|---|
| the ALT allele | **present** | detailed TSV, `evidence` names the anchor used |
| the REF allele | **absent** | nowhere |
| anything else | **undetermined**, `reason` = `other_allele` | undetermined TSV, with the bases in `observed_allele` |
| nothing usable | **undetermined**, `reason` = `not_found` | undetermined TSV |
| the two anchors disagree | **undetermined**, `reason` = `anchors_disagree` | undetermined TSV |

An allele is read only from an anchor that occurs **exactly once** in the
genome; a repeated anchor says nothing about which copy was read. Only
substitutions get anchors: with an indel the allele's own length is what is in
question, so bases counted from one anchor cannot settle it, and those markers
keep the k-mer answer.

The pass runs only when the k-mers left something unread, so a genome the panel
reads whole costs nothing extra.

**Both alleles** means two copies of the region disagree, or the ALT k-mer also
occurs somewhere else in the genome. Treat such a call with care: a marker
whose ALT k-mer matches a paralogous region is found in every genome.

**Undetermined** means the position could not be read even from the anchors: a
change on both sides of the variant, a region missing from the assembly, a
contig break, or a repeated anchor. It does not mean the variant is absent, and
it does not mean it is present.

!!! warning "Reading undetermined markers in a resistance panel"
    A drug with at least one marker **present** is called the same way whatever
    else is undetermined: a found variant hides its neighbours inside the same
    window, so a resistant genome usually lists a few undetermined markers next
    to the one it carries.

    A drug with **no marker present but some undetermined** is the case to look
    at. Absence of those markers was not established, so the genome cannot be
    reported susceptible to that drug on the strength of this run alone. The
    `undetermined:count` column of the summary gives the count per drug, and
    the undetermined TSV lists every marker involved.

### Detailed TSV columns

Header: `genome	k-mer	k-merPOS	SNPgenome	SNPreference	REF	ALT	lineage	Gene	Gene_Start	Gene_End	AA_Pos	AA_Change	alleles_found`

| Column | Description |
|---|---|
| `genome` | Sample/genome name. |
| `k-mer` | The matched marker k-mer sequence. |
| `k-merPOS` | 1-based start position of the k-mer on the strand that produced the match (reverse-strand hits are numbered on the reverse complement). |
| `SNPgenome` | 1-based position of the allele within the query genome. |
| `SNPreference` | Position the marker refers to in the reference. |
| `REF` | Reference allele of the marker. |
| `ALT` | Alternate (marker) allele. |
| `lineage` | Lineage assigned to this marker. |
| `Gene` | Overlapping gene (from the GFF, otherwise from the marker definition if present). |
| `Gene_Start` | 1-based start coordinate of that gene. |
| `Gene_End` | End coordinate of that gene. |
| `AA_Pos` | Amino-acid position of the change within the gene. |
| `AA_Change` | Amino-acid change (from the GFF annotation, otherwise from the marker definition if present). |
| `alleles_found` | `ALT` when only the variant was found; `ALT,REF` when the reference allele was found as well. |
| `evidence` | `kmer` when the whole k-mer matched, or `left_anchor`, `right_anchor` or `both_anchors` when the allele was read from the flanks. The `k-mer`, `k-merPOS` and `SNPgenome` columns are empty for an anchor read, which establishes no genome coordinate. |

`Gene_Start`, `Gene_End` and `AA_Pos` are only populated when a GFF annotation is provided and the marker falls inside an annotated feature. `Gene` and `AA_Change` fall back to the gene and mutation values from the marker TSV when present, so they can be filled even without a GFF; columns with no available value are left empty.

### Summary TSV columns

Header: `genome	lineage:count	major_lineage	undetermined	undetermined:count`

| Column | Description |
|---|---|
| `genome` | Sample/genome name. |
| `lineage:count` | Observed lineages with the number of supporting markers for each. |
| `major_lineage` | The single lineage called for the genome (respecting `--nested-classification` when enabled). |
| `undetermined` | Number of markers that could be read neither way in this genome. |
| `undetermined:count` | The same markers counted under the first level of their path: the drug for a resistance panel, the main lineage for a lineage panel. Empty when every marker was read. |

### Undetermined TSV columns

Header: `genome	SNPreference	REF	ALT	lineage	Gene	mutation	reason`

| Column | Description |
|---|---|
| `genome` | Sample/genome name. |
| `SNPreference` | Position the marker refers to in the reference. |
| `REF`, `ALT` | Alleles of the marker. |
| `lineage` | Lineage (or drug) path of the marker. |
| `Gene`, `mutation` | Gene and mutation from the marker TSV, empty when it has none. |
| `reason` | `not_found`: nothing of this marker occurs in the genome and no anchor could read it. `other_allele`: the anchors read a base other than either allele of the marker. `anchors_disagree`: the two anchors read different alleles. `no_reference_kmer`: the marker has no usable REF k-mer (its reference allele holds a base other than A, C, G or T) and the anchors could not read it either. |
| `observed_allele` | With `other_allele`, the bases actually found where the allele should be. Empty otherwise. |

## See also

- [Input formats](input-formats.md) — marker TSV, reference FASTA, sample list, and GFF layouts.
- [Assembly classification](algorithms/assembly-classification.md) — marker k-mer generation, encoded vs. text index, GFF annotation, indel support, and lineage-calling modes in depth.
- [Command reference](index.md) — all `pathotypr` subcommands.
