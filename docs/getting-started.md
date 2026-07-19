# Getting started with MTBC

This end-to-end tutorial takes you from a fresh install to a lineage call for the *Mycobacterium tuberculosis* complex (MTBC), using the published pathotypr marker panel and pre-trained model. You will classify an assembled genome, genotype raw reads, run the machine-learning predictor, and read every output file.

!!! abstract "What you'll do"
    1. [Install pathotypr](#1-install-pathotypr)
    2. [Get the MTBC panel, model, and reference](#2-get-the-mtbc-panel-model-and-reference)
    3. [Classify an assembly with `pathotypr classify`](#3-classify-an-assembly)
    4. [Genotype raw reads with `pathotypr split-fastq`](#4-genotype-raw-reads)
    5. [Predict a lineage with `pathotypr predict`](#5-predict-a-lineage-with-the-model)
    6. [Read the output](#6-read-the-output)

The MTBC panel resolves lineages **L1–L10** and animal-adapted **A1–A4**. Two complementary routes are shown: **marker-based** calling (`classify` / `split-fastq`, exact known SNPs) and **model-based** prediction (`predict`, a Random Forest over k-mer features). Run either or both.

| Your data | Command | Route |
|---|---|---|
| Assembled genome (FASTA) | [`classify`](classify.md) | Marker-based |
| Raw reads (FASTQ) | [`split-fastq`](split-fastq.md) | Marker-based |
| Assembled genome (FASTA) | [`predict`](predict.md) | Model-based |

---

## 1. Install pathotypr

Pick any install route from the [**Installation**](installation.md) guide — the desktop app, Bioconda, or a build from source. The fastest for the CLI is Bioconda:

```bash
conda create -n pathotypr -c bioconda pathotypr
conda activate pathotypr
```

Confirm the binary is on your `PATH`:

```bash
pathotypr --version
pathotypr --help
```

!!! tip "Verbose logging"
    Every subcommand accepts `-v` (debug) or `-vv` (trace). Add it to any command below if you need to see what pathotypr is doing under the hood.

---

## 2. Get the MTBC panel, model, and reference

The marker panel and pre-trained model live on Zenodo ([record 19210044](https://zenodo.org/records/19210044), DOI [10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)). The marker positions are defined against the **H37Rv** reference genome (`GCF_000195955.2`), which you fetch from NCBI.

=== "curl"

    ```bash
    # Lineage SNP panel (3,707 markers, L1–L10 + A1–A4)
    curl -L -o pathotypr_lineage_markers_v1.0.0.tsv \
      "https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1"

    # Drug-resistance panel (WHO catalogue 2021)
    curl -L -o pathotypr_dr_markers_v1.0.0.tsv \
      "https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1"

    # Pre-trained Random Forest model (k=31, 100 trees)
    curl -L -o pathotypr_rf_model_v1.0.0.pathotypr \
      "https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1"

    # H37Rv reference genome (single-record)
    curl -L -o h37rv.fasta.gz \
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
    gunzip h37rv.fasta.gz
    ```

=== "wget"

    ```bash
    # Lineage SNP panel (3,707 markers, L1–L10 + A1–A4)
    wget -O pathotypr_lineage_markers_v1.0.0.tsv \
      "https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1"

    # Drug-resistance panel (WHO catalogue 2021)
    wget -O pathotypr_dr_markers_v1.0.0.tsv \
      "https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1"

    # Pre-trained Random Forest model (k=31, 100 trees)
    wget -O pathotypr_rf_model_v1.0.0.pathotypr \
      "https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1"

    # H37Rv reference genome (single-record)
    wget -O h37rv.fasta.gz \
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
    gunzip h37rv.fasta.gz
    ```

!!! info "Why H37Rv?"
    `classify` and `split-fastq` need a **single-record** reference FASTA whose coordinates match the marker positions. H37Rv (`GCF_000195955.2`, chromosome `NC_000962.3`) is a single contig with no plasmids, so it satisfies this directly. The model in `pathotypr_rf_model_v1.0.0.pathotypr` is self-contained and does **not** need the reference.

After this step your working directory holds the three Zenodo assets plus `h37rv.fasta`. You'll also need your own sample(s): an **assembled genome** (`sample.fasta`) and, for the reads route, a **FASTQ pair** (`sample_R1.fastq.gz`, `sample_R2.fastq.gz`). See [Input formats](input-formats.md) for the exact specs.

---

## 3. Classify an assembly

Call the known lineage SNPs in an assembled genome. `--nested-classification` turns on hierarchical calling so the multi-level marker columns resolve sub-lineages, not just the top level.

=== "Lineage panel"

    ```bash
    pathotypr classify \
      -m pathotypr_lineage_markers_v1.0.0.tsv \
      -r h37rv.fasta \
      -i sample.fasta \
      -o sample_lineage \
      --nested-classification \
      --excel
    ```

=== "Drug-resistance panel"

    ```bash
    # Same command, swap in the DR panel to scan for resistance mutations
    pathotypr classify \
      -m pathotypr_dr_markers_v1.0.0.tsv \
      -r h37rv.fasta \
      -i sample.fasta \
      -o sample_dr \
      --excel
    ```

This writes `sample_lineage.tsv` (one row per marker hit), `sample_lineage_summary.tsv` (one row for the genome with its major lineage), and — because of `--excel` — the matching `.xlsx` files.

!!! tip "Gene and amino-acid change columns"
    Pass a GFF3 annotation with `--gff` (for `-i`) to populate the `Gene`, `Gene_Start`, `Gene_End`, `AA_Pos`, and `AA_Change` columns — especially useful with the DR panel to name the affected gene and mutation. Without a GFF those columns stay empty. See [`classify`](classify.md) for GFF handling.

!!! note "Many genomes at once"
    Swap `-i sample.fasta` for `-l genomes.tsv` (a `sample_name<TAB>fasta_path[<TAB>gff_path]` list) or `--input-files a.fasta b.fasta …` to run a whole batch under one output prefix.

---

## 4. Genotype raw reads

No assembly? Genotype straight from FASTQ. `split-fastq` counts REF vs ALT diagnostic k-mers at each marker — no alignment step. Paired files named `_R1/_R2` (or `_1/_2`) are auto-detected; `--paired` makes the pairing explicit.

```bash
pathotypr split-fastq \
  -m pathotypr_lineage_markers_v1.0.0.tsv \
  -r h37rv.fasta \
  -i sample_R1.fastq.gz -i sample_R2.fastq.gz \
  --paired \
  -o sample_genotype \
  --nested-classification \
  --excel
```

Outputs are `sample_genotype_<sample>_mutations.tsv` (per-marker detail: position, alleles, counts, ALT fraction, lineage) and `sample_genotype_summary.tsv` (lineage counts plus the major call).

!!! tip "Tuning calls for coverage"
    A variant is called when depth ≥ `--min-depth` (default `10`) **and** the ALT fraction ≥ `--min-alt-percent` (default `95`). Lower `--min-depth` for shallow sequencing, or lower `--min-alt-percent` if you expect mixed / heteroresistant samples.

!!! warning "Indels are skipped from reads"
    Short reads give unreliable k-mer matches across MTBC repeats (PE/PPE, IS elements), so `split-fastq` calls SNPs and MNVs only. For indel-bearing markers, use the assembly route in [Step 3](#3-classify-an-assembly).

---

## 5. Predict a lineage with the model

`predict` is the model-based alternative: it vectorizes the assembly with the same feature hasher used at training time and lets the Random Forest vote. It needs only your FASTA and the downloaded model — no markers, no reference.

```bash
pathotypr predict \
  -i sample.fasta \
  -m pathotypr_rf_model_v1.0.0.pathotypr \
  -o sample_prediction.tsv \
  --excel
```

This writes `sample_prediction.tsv` (plus `sample_prediction.xlsx`) with one row per input sequence.

!!! question "classify or predict?"
    `classify` reports exactly **which** known markers were hit — auditable, and it also flags resistance mutations with the DR panel. `predict` needs no reference or markers and returns a lineage with a calibrated confidence. Running both is a good cross-check.

---

## 6. Read the output

### `classify` — detailed TSV

`sample_lineage.tsv` has one row per marker match, with these columns:

| Column | Meaning |
|---|---|
| `genome` | Sample / contig the hit came from |
| `k-mer` | Diagnostic k-mer that matched |
| `k-merPOS` | Position of the match within the genome |
| `SNPgenome` | Allele observed in the genome |
| `SNPreference` | Allele in the reference |
| `REF` / `ALT` | Marker reference and alternate alleles |
| `lineage` | Lineage label carried by the marker |
| `Gene`, `Gene_Start`, `Gene_End` | Gene context (populated only with `--gff`) |
| `AA_Pos`, `AA_Change` | Amino-acid position and change (only with `--gff`) |

The companion `sample_lineage_summary.tsv` collapses this to a single row: the genome and its **major lineage** call (with sub-lineage when `--nested-classification` is on).

### `split-fastq` — mutations + summary

`sample_genotype_<sample>_mutations.tsv` lists each called marker with its position, REF/ALT alleles, REF and ALT read counts, ALT fraction, and lineage. `sample_genotype_summary.tsv` gives one row per sample: the per-lineage marker counts and the final major-lineage call.

### `predict` — predictions TSV

`sample_prediction.tsv` has one row per input sequence:

| Column | Meaning |
|---|---|
| `Header` | Sequence header from the input FASTA |
| `Predicted_Lineage` | Most-voted lineage label |
| `Confidence` | Fraction of trees voting for the winner (0–1) |
| `Confidence_Margin` | Gap between winner and runner-up (0–1) |
| `Other_Votes` | Top 3 alternatives with their vote fractions |

!!! example "Reading a confidence score"
    A `Confidence` of `0.95` means 95% of the trees agreed on the call; a `Confidence_Margin` of `0.40` means the winner drew 40% more votes than the runner-up. High values on both signal a clean, unambiguous prediction.

!!! tip "Excel everywhere"
    Every command above used `--excel`, so each `.tsv` has a side-by-side `.xlsx`. Drop the flag for TSV-only output.

---

## Next steps

<div class="grid cards" markdown>

-   :material-console:{ .lg .middle } &nbsp; [**`classify`**](classify.md)

    ---

    Full options, GFF annotation, and the lineage-calling logic.

-   :material-dna:{ .lg .middle } &nbsp; [**`split-fastq`**](split-fastq.md)

    ---

    Read-based genotyping, depth/frequency thresholds, and pairing modes.

-   :material-brain:{ .lg .middle } &nbsp; [**`predict`**](predict.md)

    ---

    Model-based prediction, confidence metrics, and output details.

-   :material-school:{ .lg .middle } &nbsp; [**`train`**](train.md)

    ---

    Build your own Random Forest model from labeled genomes.

</div>

- Bringing your own organism or panel? See the [marker format](marker_format.md) and [input formats](input-formats.md) references.
- Prefer clicking to typing? The same workflow runs in the [Desktop GUI](gui.md).
