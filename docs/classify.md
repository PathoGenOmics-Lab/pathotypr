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

Exactly one input source is required: a single FASTA (`-i`), a TSV list of samples (`-l`), or several FASTA files (`--input-files`).

## Options

| Option | Default | Required | Description |
|---|---|---|---|
| `-m, --markers <FILE>` | — | yes | Marker definition TSV (position, REF, ALT, lineage levels...). |
| `-r, --reference <FILE>` | — | yes | Single-record reference FASTA the marker positions refer to. |
| `-i, --input <FILE>` | — | one of `-i`/`-l`/`--input-files` | A single sample FASTA. Each FASTA record is treated as a separate genome (one output row per contig). |
| `-l, --input-list <FILE>` | — | one of `-i`/`-l`/`--input-files` | TSV of samples: `sample_name\tfasta_path[\tgff_path]`. All contigs of a sample's FASTA are aggregated under its sample name. |
| `--input-files <FILE>...` | — | one of `-i`/`-l`/`--input-files` | Multiple sample FASTA files (GUI batch mode). |
| `--gff <FILE>` | — | no | GFF annotation for `--input` (adds gene + amino-acid change columns). Requires `--input`. |
| `--gff-files <FILE>...` | — | no | Multiple GFFs matched by filename, for `--input-files`. |
| `-o, --output-prefix <PREFIX>` | — | yes | Prefix for all output files. |
| `--kmer-size <N>` | `31` | no | Marker k-mer size. |
| `-t, --threads <N>` | all cores | no | Number of CPU threads. |
| `--nested-classification` | `false` | no | Enable hierarchical (nested) lineage calling using multi-level marker columns. |
| `--min-flank-bases <N>` | `10` | no | Minimum flanking bases on each side of the allele in a marker k-mer. |
| `--output-masked-fasta` | `false` | no | Also write masked FASTA files with marker positions replaced by `N`. |
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

1. The reference FASTA and marker TSV are loaded, and one diagnostic k-mer is built per marker: the ALT allele flanked on both sides by reference context (at least `--min-flank-bases` bases per side, total length `--kmer-size`).
2. Every contig of every query genome is scanned for exact matches to those marker k-mers.
3. If a GFF is supplied for a sample, each match is annotated with the overlapping gene and the resulting amino-acid change.
4. Per genome, the observed markers are tallied per lineage and a `major_lineage` is called. With `--nested-classification`, multi-level marker columns are used to validate the hierarchical lineage path; without it, the most abundant lineage is reported directly.
5. Results are streamed to a detailed per-marker TSV and a per-genome summary TSV (plus optional `.xlsx` and masked FASTA outputs).

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
| `classify_run_summary.tsv` | always | Summary table: one row per genome with lineage counts and the major lineage. |
| `classify_run.xlsx`, `classify_run_summary.xlsx` | `--excel` | Excel (`.xlsx`) versions of the detailed and summary tables, written alongside the TSVs. |
| `<input>_masked.fasta` | `--output-masked-fasta` | One masked FASTA per input FASTA file (named after the input file stem), with marker positions replaced by `N`. |

!!! note
    If the output prefix already ends in `.tsv`, the detailed file uses it as-is and the summary is derived from the trimmed base name. Genomes with no marker hits still appear in the detailed TSV as a single row with empty match columns, so every input is accounted for.

### Detailed TSV columns

Header: `genome	k-mer	k-merPOS	SNPgenome	SNPreference	REF	ALT	lineage	Gene	Gene_Start	Gene_End	AA_Pos	AA_Change`

| Column | Description |
|---|---|
| `genome` | Sample/genome name. |
| `k-mer` | The matched marker k-mer sequence. |
| `k-merPOS` | 1-based start position of the k-mer within the query genome. |
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

`Gene_Start`, `Gene_End` and `AA_Pos` are only populated when a GFF annotation is provided and the marker falls inside an annotated feature. `Gene` and `AA_Change` fall back to the gene and mutation values from the marker TSV when present, so they can be filled even without a GFF; columns with no available value are left empty.

### Summary TSV columns

Header: `genome	lineage:count	major_lineage`

| Column | Description |
|---|---|
| `genome` | Sample/genome name. |
| `lineage:count` | Observed lineages with the number of supporting markers for each. |
| `major_lineage` | The single lineage called for the genome (respecting `--nested-classification` when enabled). |

## See also

- [Input formats](input-formats.md) — marker TSV, reference FASTA, sample list, and GFF layouts.
- [Assembly classification](algorithms/assembly-classification.md) — marker k-mer generation, encoded vs. text index, GFF annotation, indel support, and lineage-calling modes in depth.
- [Command reference](index.md) — all `pathotypr` subcommands.
