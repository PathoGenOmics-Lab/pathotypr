# Algorithm Documentation

This section describes the algorithms and data structures behind each pathotypr
module. For CLI usage and options, see the command guides:
[train](../train.md), [predict](../predict.md), [classify](../classify.md),
[split-fastq](../split-fastq.md), and [match](../match.md).

| Module | Document | Core Idea |
|---|---|---|
| Feature Hashing | [feature-hashing.md](feature-hashing.md) | The hashing trick: k-mers → fixed-size sparse vectors |
| Random Forest | [random-forest.md](random-forest.md) | Sparse CART trees with bootstrap aggregation |
| Training Pipeline | [training.md](training.md) | End-to-end: vectorize → evaluate → train → OOB → export |
| Prediction | [prediction.md](prediction.md) | Streaming batch prediction with majority voting |
| Marker Genotyping | [marker-genotyping.md](marker-genotyping.md) | Diagnostic k-mers + Bloom filter for FASTQ scanning |
| Reference Matching | [reference-matching.md](reference-matching.md) | K-mer containment scoring with streaming batches |
| Assembly Classification | [assembly-classification.md](assembly-classification.md) | Marker calling on FASTA assemblies with GFF annotation |
