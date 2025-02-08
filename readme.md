# Pathotypr

**Pathotypr** is a powerful command-line tool for genome classification using machine learning and SNP markers. It provides three main functionalities:

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
git clone https://github.com/yourusername/pathotypr.git
cd pathotypr
cargo build --release
```

## Quick Start

```bash
# Train a model
pathotypr train --fasta input.fasta --output my_model --kmer_size 21

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

## Contact & Support

- **Author**: Paula Ruiz Rodriguez
- **Email**: paula.ruiz.rodriguez@csic.es
- **License**: GPL-3
- **Issues**: Submit via GitHub issue tracker

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Acknowledgments

- SmartCore team
- Rust Bioinformatics community
- Crate maintainers