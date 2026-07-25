<p align="center">
  <img src="logo/pathotypr.svg" alt="pathotypr logo" width="750" />
</p>

<div align="center">

[![Documentation](https://img.shields.io/badge/docs-pathotypr-%23b01000?style=flat-square)](https://pathogenomics-lab.github.io/pathotypr/)
[![License: AGPL v3](https://img.shields.io/badge/license-AGPL%20v3-%23af64d1?style=flat-square)](LICENSE)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pathotypr?style=flat-square&color=%23009E73&label=bioconda)](https://anaconda.org/bioconda/pathotypr)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/pathotypr?style=flat-square&color=%23009E73&label=downloads)](https://anaconda.org/bioconda/pathotypr)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-%23b31b1b?style=flat-square)](https://doi.org/10.64898/2026.03.24.714002)
[![Markers](https://img.shields.io/badge/markers-Zenodo%2019210044-%230072B2?style=flat-square)](https://zenodo.org/records/19210044)

**Lineage classification and marker-driven genotyping — from assemblies or raw reads.**

[Quick Start](#quick-start) · [Commands](#commands) · [GUI](#gui) · [Documentation](https://pathogenomics-lab.github.io/pathotypr/) · [Citation](#citation)

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## What is pathotypr?

pathotypr is a Rust toolkit that classifies microbial genomes into lineages and genotypes them against user-defined marker panels. It works with both assembled genomes (FASTA) and raw sequencing reads (FASTQ), runs on a single laptop, and ships with a native desktop GUI.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/scheme-dark.webp">
    <source media="(prefers-color-scheme: light)" srcset="docs/assets/scheme.svg">
    <img src="docs/assets/scheme.svg" alt="pathotypr workflow schema" width="1200" />
  </picture>
</p>

**Five commands, one binary:**

| Command | What it does | Input |
|---|---|---|
| **`train`** | Build a Random Forest classifier from labeled genomes | FASTA |
| **`predict`** | Assign lineages using a trained model | FASTA + model |
| **`classify`** | Call known SNP markers in assemblies | FASTA + markers |
| **`split-fastq`** | Alignment-free genotyping from reads | FASTQ + markers |
| **`match`** | Find the closest reference genome | FASTQ + references |

**Key features:**
- 🦠 **Organism-agnostic** — bring your own markers for any pathogen
- ⚡ **Fast** — Rust + SIMD gzip + parallel k-mers (~1–2 s per sample)
- 🖥️ **Desktop GUI** — native app via Tauri, no server required
- 📊 **Excel + TSV** output with interactive visualizations in the GUI

## Installation

### Desktop GUI (pre-built)

**[⬇ Download the latest release](https://github.com/PathoGenOmics-Lab/pathotypr/releases/latest)**, then pick the file for your platform:

| Platform | File | Notes |
|---|---|---|
| 🍎 macOS (Apple Silicon) | `Pathotypr_<version>_aarch64.dmg` | M1 / M2 / M3 / M4 Macs |
| 🍎 macOS (Intel) | `Pathotypr_<version>_x64.dmg` | Pre-2020 Macs |
| 🐧 Linux (.deb) | `Pathotypr_<version>_amd64.deb` | Debian / Ubuntu |
| 🐧 Linux (.rpm) | `Pathotypr-<version>-1.x86_64.rpm` | Fedora / RHEL |
| 🐧 Linux (AppImage) | `Pathotypr_<version>_amd64.AppImage` | Any distro, no install needed |
| 🪟 Windows (installer) | `Pathotypr_<version>_x64-setup.exe` | Windows 10+ |
| 🪟 Windows (.msi) | `Pathotypr_<version>_x64_en-US.msi` | Windows 10+ (MSI) |

> [!NOTE]
> **macOS users**: The app is not signed with an Apple Developer certificate. On first launch, right-click the app → **Open** → click **Open** in the dialog. See [Apple support](https://support.apple.com/en-us/HT202491) for details.
>
> **Windows users**: Windows SmartScreen may show a warning for unrecognized apps. Click **More info** → **Run anyway** to proceed.

Older versions are on the [releases page](https://github.com/PathoGenOmics-Lab/pathotypr/releases).

### CLI (Bioconda)

```bash
conda create -n pathotypr -c bioconda pathotypr
conda activate pathotypr
pathotypr --help
```

### CLI (from source)

```bash
git clone https://github.com/PathoGenOmics-Lab/pathotypr.git
cd pathotypr
cargo build --release -p pathotypr-core --bin pathotypr
./target/release/pathotypr --help
```

### GUI (from source)

See [docs/gui.md](https://pathogenomics-lab.github.io/pathotypr/gui/) for building the Tauri desktop app from source.

## MTBC Marker Files & Pre-trained Model

Ready-to-use marker panels and a pre-trained Random Forest model for *Mycobacterium tuberculosis* complex (MTBC) are available on Zenodo:

| File | Description | Download |
|---|---|---|
| `pathotypr_lineage_markers_v1.0.0.tsv` | 3,707 lineage SNPs (L1–L10, A1–A4) | [⬇ Download](https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1) |
| `pathotypr_dr_markers_v1.0.0.tsv` | 102,213 DR mutations (WHO catalogue v2, 2023) | [⬇ Download](https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1) |
| `pathotypr_rf_model_v1.0.0.pathotypr` | Pre-trained RF model (k=31, 100 trees) | [⬇ Download](https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1) |

> **DOI:** [10.5281/zenodo.19210044](https://zenodo.org/records/19210044)

## Quick Start

```bash
# Train a lineage model
pathotypr train -i labeled_genomes.fasta -o model.pathotypr.zst

# Predict lineages
pathotypr predict -i query.fasta -m model.pathotypr.zst -o predictions.tsv

# Classify markers in assemblies
pathotypr classify -m markers.tsv -r reference.fasta -i sample.fasta -o results

# Genotype from FASTQ reads
pathotypr split-fastq -m markers.tsv -r reference.fasta \
  -i reads_R1.fastq.gz -i reads_R2.fastq.gz --paired -o genotype

# Find best reference match
pathotypr match -i reads_R1.fastq.gz reads_R2.fastq.gz \
  -r references.fasta -o match.tsv
```

Add `--excel` to any command to also generate `.xlsx` files.

## Commands

Each command has its own detailed documentation:

| Command | Docs | Summary |
|---|---|---|
| `train` | [docs/train.md](https://pathogenomics-lab.github.io/pathotypr/train/) | Random Forest on k-mer feature-hashed vectors |
| `predict` | [docs/predict.md](https://pathogenomics-lab.github.io/pathotypr/predict/) | Streaming batch prediction with confidence scores |
| `classify` | [docs/classify.md](https://pathogenomics-lab.github.io/pathotypr/classify/) | Marker k-mer matching + GFF annotation + masked FASTA |
| `split-fastq` | [docs/split-fastq.md](https://pathogenomics-lab.github.io/pathotypr/split-fastq/) | Alignment-free genotyping with Bloom filter acceleration |
| `match` | [docs/match.md](https://pathogenomics-lab.github.io/pathotypr/match/) | K-mer containment scoring against reference databases |

Run `pathotypr <command> --help` for all options.

### Algorithm Details

For in-depth descriptions of the algorithms, data structures, and design decisions behind each module, see [docs/algorithms/](https://pathogenomics-lab.github.io/pathotypr/algorithms/):

| Document | Topic |
|---|---|
| [Feature Hashing](https://pathogenomics-lab.github.io/pathotypr/algorithms/feature-hashing/) | The hashing trick: k-mers → fixed-size sparse vectors |
| [Random Forest](https://pathogenomics-lab.github.io/pathotypr/algorithms/random-forest/) | Sparse CART trees with bootstrap aggregation |
| [Training Pipeline](https://pathogenomics-lab.github.io/pathotypr/algorithms/training/) | Vectorize → evaluate → train → OOB → export |
| [Prediction](https://pathogenomics-lab.github.io/pathotypr/algorithms/prediction/) | Streaming batch prediction with majority voting |
| [Marker Genotyping](https://pathogenomics-lab.github.io/pathotypr/algorithms/marker-genotyping/) | Diagnostic k-mers + Bloom filter for FASTQ scanning |
| [Reference Matching](https://pathogenomics-lab.github.io/pathotypr/algorithms/reference-matching/) | K-mer containment scoring with streaming batches |
| [Assembly Classification](https://pathogenomics-lab.github.io/pathotypr/algorithms/assembly-classification/) | Marker calling on FASTA with GFF annotation |

## Input Formats

### Training FASTA

The first token in each header is the class label:

```
>L4 sample_0001
ACTG...
>L2 sample_0002
ACTG...
```

### Marker TSV

Tab-separated: `position  REF  ALT  level1  [level2  ...]`

```
#pos    ref    alt    level1    level2
761155  C      T      L4        L4.9
2155168 G      A      L2        L2.2
```

Lineage columns are read until the first empty cell. Columns after the empty cell are treated as annotations.

> See [docs/input-formats.md](https://pathogenomics-lab.github.io/pathotypr/input-formats/) for full format specifications.

## GUI

The desktop app includes all five workflows with drag-and-drop file selection, interactive result tables, and real-time progress indicators.

```bash
# Development
cargo tauri dev

# Production build
cargo tauri build
```

> See [docs/gui.md](https://pathogenomics-lab.github.io/pathotypr/gui/) for system dependencies and build instructions.

## Performance

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="benchmarks/figures/dashboard-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="benchmarks/figures/dashboard.png">
    <img src="benchmarks/figures/dashboard.png" alt="pathotypr performance benchmarks" width="900" />
  </picture>
</p>

Benchmarked on real *M. tuberculosis* genomes (~4.4 Mb, k=21), Mac mini M4, 4 threads:

| Module | Time | Peak RAM | Key property |
|---|---|---|---|
| **train** (10 genomes) | 0.6 s | 302 MB | Scales with dataset size |
| **train** (50 genomes) | 55 s | 1.4 GB | |
| **predict** (5 genomes) | 0.25 s | 198 MB | ~50 ms/genome, constant |
| **classify** (5 genomes) | 0.10 s | 92 MB | ~20 ms/genome |
| **split-fastq** (65× PE) | 10.5 s | 26 MB | **Constant memory** |
| **match** (20 refs) | 78 s | 4.6 GB | Streaming batches |

- **SIMD-accelerated** gzip decompression (zlib-ng)
- **Streaming** I/O — split-fastq holds 26 MB regardless of input size
- **85,000+ genomes/second** prediction throughput (synthetic benchmarks)

> See [docs/benchmarks.md](https://pathogenomics-lab.github.io/pathotypr/benchmarks/) for detailed charts, scaling plots, and pathotypr vs fastlin comparison.

## Comparison

| | pathotypr | fastlin | TB-Profiler | Mykrobe | SNP-IT | KvarQ |
|---|---|---|---|---|---|---|
| Alignment-free (FASTQ) | ✅ | ✅ | ❌ | ✅ | ❌ | ✅ |
| Assemblies (FASTA) | ✅ | ✅ | ❌ | ❌ | VCF only | ❌ |
| Custom markers | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| ML training | ✅ | ❌ | ❌ | ❌ | ✅ | ❌ |
| DR prediction | ✅ | ❌ | ✅ | ✅ | ❌ | ✅ |
| Desktop GUI | ✅ | ❌ | Web | ✅ | ❌ | ✅ |
| Standalone binary | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| Organism-agnostic | ✅ | TB only | TB only | Limited | TB only | TB only |
| Speed (per sample) | ~1 s | <5 s | 3–10 min | ~3 min | 1–2 min | ~2 min |

## Project Structure

```
pathotypr/
├── pathotypr-core/           # Core library + CLI
│   └── src/
│       ├── main.rs           # CLI entry point
│       ├── train.rs          # Random Forest training + OOB + CV
│       ├── predict.rs        # Streaming batch prediction
│       ├── classify/         # Assembly-based marker classification
│       │   ├── mod.rs        #   Orchestration + genome analysis
│       │   ├── markers.rs    #   Marker parsing + k-mer generation
│       │   ├── annotation.rs #   GFF parsing + AA translation
│       │   └── masking.rs    #   FASTA masking at marker sites
│       ├── classify_split_fastq.rs  # FASTQ genotyping orchestration
│       ├── split_kmer.rs     # Diagnostic k-mer engine + Bloom filter
│       ├── match/            # Reference matching
│       │   ├── mod.rs        #   Scoring + coarse-to-fine matching
│       │   └── index.rs      #   Compact inverted index + cache
│       ├── sparse_tree.rs    # Custom CART on sparse vectors
│       ├── vectorizer.rs     # Feature hashing (hashing trick)
│       ├── model.rs          # Model bundle + label encoder
│       ├── lineage.rs        # Hierarchical lineage classification
│       ├── fasta_io.rs       # FASTA reading (needletail)
│       ├── paired_end.rs     # Paired-end FASTQ detection
│       ├── excel.rs          # Streaming Excel export
│       ├── errors.rs         # Error types + cancellation
│       └── common.rs         # Thread pool + shared utilities
├── src-tauri/                # Desktop app backend (Tauri)
├── frontend/                 # GUI (HTML/CSS/JS)
├── docs/                     # Detailed documentation
└── logo/                     # Branding assets
```

## Citation

If you use pathotypr, please cite:

> Ruiz-Rodriguez P, Coscollá M. **Pathotypr: harmonised MTBC lineage assignment and resistance-associated variant detection for genomic surveillance.** *bioRxiv* (2026). doi: [10.64898/2026.03.24.714002](https://doi.org/10.64898/2026.03.24.714002)

```bibtex
@article{ruiz-rodriguez_pathotypr_2026,
  title     = {Pathotypr: harmonised {MTBC} lineage assignment and resistance-associated variant detection for genomic surveillance},
  author    = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  journal   = {bioRxiv},
  year      = {2026},
  doi       = {10.64898/2026.03.24.714002},
  url       = {https://www.biorxiv.org/content/10.64898/2026.03.24.714002v1}
}
```

> Software & markers DOI: [10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)

## License

[GNU Affero General Public License v3.0](LICENSE)

---

<h2 id="contributors" align="center">

✨ [Contributors](https://github.com/PathoGenOmics-Lab/pathotypr/graphs/contributors)
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
pathotypr is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
---
