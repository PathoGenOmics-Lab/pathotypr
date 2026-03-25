# Documentation

## Command Reference

Usage guides, CLI options, and output formats for each command.

| Command | Description |
|---|---|
| [train.md](train.md) | Train a Random Forest classifier from labeled genomes |
| [predict.md](predict.md) | Predict lineages using a trained model |
| [classify.md](classify.md) | Call known SNP markers in assembled genomes |
| [split-fastq.md](split-fastq.md) | Alignment-free genotyping from raw FASTQ reads |
| [match.md](match.md) | Find the closest reference genome for a sample |

## Guides

| Document | Description |
|---|---|
| [input-formats.md](input-formats.md) | Format specifications for FASTA, marker TSV, and input lists |
| [gui.md](gui.md) | Building and using the Tauri desktop application |
| [benchmarks.md](benchmarks.md) | Performance benchmarks: speed, throughput, memory, and comparisons |

## Algorithm Details

In-depth descriptions of the algorithms, data structures, and design decisions behind each module.

| Document | Topic |
|---|---|
| [algorithms/feature-hashing.md](algorithms/feature-hashing.md) | The hashing trick: k-mers → fixed-size sparse vectors |
| [algorithms/random-forest.md](algorithms/random-forest.md) | Sparse CART trees with bootstrap aggregation and OOB accuracy |
| [algorithms/training.md](algorithms/training.md) | End-to-end pipeline: vectorize → CV/split → train → serialize |
| [algorithms/prediction.md](algorithms/prediction.md) | Streaming batch prediction with majority voting |
| [algorithms/marker-genotyping.md](algorithms/marker-genotyping.md) | Diagnostic k-mers + Bloom filter for FASTQ scanning |
| [algorithms/reference-matching.md](algorithms/reference-matching.md) | K-mer containment scoring with streaming batches |
| [algorithms/assembly-classification.md](algorithms/assembly-classification.md) | Marker calling on FASTA assemblies with GFF annotation |
