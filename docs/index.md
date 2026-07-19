# pathotypr

**Lineage classification and marker-driven genotyping — from assemblies or raw reads.**

pathotypr is a Rust toolkit that classifies microbial genomes into lineages and
genotypes them against user-defined marker panels. It works with both assembled
genomes (FASTA) and raw sequencing reads (FASTQ), runs on a single laptop, and
ships with a native desktop GUI.

![pathotypr workflow](assets/scheme.svg#only-light){ width="900" }
![pathotypr workflow](assets/scheme-dark.png#only-dark){ width="900" }

## Five commands, one binary

| Command | What it does | Input |
|---|---|---|
| [**`train`**](train.md) | Build a Random Forest classifier from labeled genomes | FASTA |
| [**`predict`**](predict.md) | Assign lineages using a trained model | FASTA + model |
| [**`classify`**](classify.md) | Call known SNP markers in assemblies | FASTA + markers |
| [**`split-fastq`**](split-fastq.md) | Alignment-free genotyping from reads | FASTQ + markers |
| [**`match`**](match.md) | Find the closest reference genome | FASTQ + references |

## Highlights

<div class="grid cards" markdown>

-   :material-bacteria-outline:{ .lg .middle } &nbsp; **Organism-agnostic**

    ---

    Bring your own markers for any pathogen — nothing is hard-coded to a
    single organism.

-   :material-lightning-bolt-outline:{ .lg .middle } &nbsp; **Fast**

    ---

    Rust + SIMD gzip + parallel k-mers process a sample in roughly
    **1–2 seconds**.

-   :material-monitor-dashboard:{ .lg .middle } &nbsp; **Desktop GUI**

    ---

    A native app built with Tauri — drag-and-drop files, no server required.

-   :material-table-large:{ .lg .middle } &nbsp; **Excel + TSV output**

    ---

    Tabular results everywhere, with interactive visualizations in the GUI.

</div>

## Quick start

=== "Assemblies (FASTA)"

    ```bash
    # Call known markers in an assembled genome
    pathotypr classify \
      -m markers.tsv -r reference.fasta \
      -i sample.fasta -o results --excel

    # Or predict a lineage with a trained model
    pathotypr predict -i query.fasta -m model.pathotypr.zst -o predictions.tsv
    ```

=== "Raw reads (FASTQ)"

    ```bash
    # Alignment-free genotyping straight from reads
    pathotypr split-fastq \
      -m markers.tsv -r reference.fasta \
      -i reads_R1.fastq.gz -i reads_R2.fastq.gz --paired \
      -o genotype

    # Find the closest reference genome
    pathotypr match -i reads_R1.fastq.gz reads_R2.fastq.gz \
      -r references.fasta -o match.tsv
    ```

=== "Train a model"

    ```bash
    # Train a Random Forest lineage classifier from labeled genomes
    pathotypr train -i labeled_genomes.fasta -o model.pathotypr.zst
    ```

!!! tip "Excel output"
    Add `--excel` to any command to also generate `.xlsx` files alongside the TSVs.

## Where to next

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } &nbsp; [**Installation**](installation.md)

    ---

    Pre-built desktop apps, Bioconda, or build the CLI/GUI from source.

-   :material-console:{ .lg .middle } &nbsp; [**Commands**](train.md)

    ---

    Detailed usage, options, and output formats for each subcommand.

-   :material-file-code-outline:{ .lg .middle } &nbsp; [**Input formats**](input-formats.md)

    ---

    Specifications for training FASTA, marker TSVs, and input lists.

-   :material-cog-outline:{ .lg .middle } &nbsp; [**Algorithms**](algorithms/README.md)

    ---

    The design and data structures behind each module.

</div>

## Citation

If you use pathotypr, please cite:

> Ruiz-Rodriguez P, Coscollá M. **Pathotypr: harmonised MTBC lineage assignment
> and resistance-associated variant detection for genomic surveillance.**
> *bioRxiv* (2026). doi: [10.64898/2026.03.24.714002](https://doi.org/10.64898/2026.03.24.714002)

Software & markers DOI: [10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)

pathotypr is released under the [GNU Affero General Public License v3.0](https://github.com/PathoGenOmics-Lab/pathotypr/blob/main/LICENSE).
