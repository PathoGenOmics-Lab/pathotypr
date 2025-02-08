
<p align="center">
  <img src="logo/pathotypr.png" title="pathotypr.png logo" style="width:750px; height: auto;">
</p>

<div align="center">
  

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>  

# pathotypr

**pathotypr** is a powerful command-line tool for genome classification using machine learning and SNP markers. It provides three main functionalities:

- 🎓 **Train**: Build ML models from FASTA sequences
- 🔮 **Predict**: Classify new genomes using trained models
- 🧬 **Classify**: Process genomic markers (SNPs) against a reference for lineage/drug resistance classification

## Installation

```bash
# Install via Conda
conda create -n pathotypr
conda activate pathotypr
conda install -c bioconda pathotypr

# Or using Mamba (faster)
mamba create -n pathotypr
mamba activate pathotypr
mamba install -c bioconda pathotypr

# Or build from source
git clone https://github.com/PathoGenOmics-Lab/pathotypr.git
cd pathotypr
cargo build --release
```

## Quick Start

```bash
# Train a model
pathotypr train --fasta input.fasta --output my_model --kmer_size 6

# Predict classifications
pathotypr predict --fasta input.fasta --model_base my_model --output predictions.txt

# Classify using SNP markers
pathotypr classify --tsv_pos markers.tsv --ref_fasta ref.fasta --fasta_genomes genomes.fasta --output results.txt
```

## Features

- 🚀 Fast parallel processing using Rayon
- 📊 Random Forest classification for ML-based prediction
- 🧪 K-mer based sequence analysis
- 💾 Compressed model storage
- 🔍 Reference-based SNP marker detection
- 🧬 Flexible marker positions for closed/mapped genomes
- 📈 Progress tracking and logging

## Documentation

### Train Mode
Trains a model from FASTA files with headers in `Lineage_sequenceID` format:
```bash
pathotypr train --fasta input.fasta --output my_model --kmer_size 21
```

### Predict Mode
Classifies genomes using a trained model:
```bash
pathotypr predict --fasta input.fasta --model_base my_model --output predictions.txt
```

### Classify Mode
Process genomes using SNP markers against a reference sequence. Supports both lineage-defining SNPs and drug resistance markers. Works with:
- Mapped genomes (SNP positions relative to reference)
- Closed genomes (SNP positions may vary)

```bash
# Using FASTA input
pathotypr classify --tsv_pos markers.tsv --ref_fasta ref.fasta --fasta_genomes genomes.fasta --output results.txt

# Using TSV input
pathotypr classify --tsv_pos markers.tsv --ref_fasta ref.fasta --tsv_genomes genomes.tsv --output results.txt
```

## Project Structure

```
pathotypr/
├── src/
│   ├── main.rs     # CLI handling
│   ├── train.rs    # Model training
│   ├── predict.rs  # Classification
│   └── classify.rs # Marker processing
└── Cargo.toml
```

## Key Dependencies

- 🎯 clap: CLI parsing
- 🤖 smartcore: Machine learning
- ⚡ rayon: Parallel processing
- 🧬 bio: Bioinformatics tools
- 📊 serde: Serialization

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Acknowledgments

- SmartCore team
- Rust Bioinformatics community
- Crate maintainers

---
<h2 id="contributors" align="center">

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
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
