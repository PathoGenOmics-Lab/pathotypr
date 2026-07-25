# Algorithm documentation

**This section explains the algorithms and data structures behind each pathotypr module — the core ideas, not the command-line switches.**

Each page below focuses on a single module: the concept it implements, the data structures it relies on, and how the pieces fit together end to end. For command-line usage and options, see the command guides in [See also](#see-also).

!!! tip
    The table is ordered roughly along the pathotypr pipeline — from turning k-mers into feature vectors, through training and prediction, to marker- and reference-based calling. Reading it top to bottom is a good first pass.

| Module | Document | Core Idea |
|---|---|---|
| Feature Hashing | [feature-hashing.md](feature-hashing.md) | The hashing trick: k-mers → fixed-size sparse vectors |
| Random Forest | [random-forest.md](random-forest.md) | Sparse CART trees with bootstrap aggregation |
| Training Pipeline | [training.md](training.md) | End-to-end: vectorize → evaluate → train → OOB → export |
| Prediction | [prediction.md](prediction.md) | Streaming batch prediction with majority voting |
| Marker Genotyping | [marker-genotyping.md](marker-genotyping.md) | Diagnostic k-mers + Bloom filter for FASTQ scanning |
| Reference Matching | [reference-matching.md](reference-matching.md) | K-mer containment scoring with streaming batches |
| Assembly Classification | [assembly-classification.md](assembly-classification.md) | Marker calling on FASTA assemblies with GFF annotation |

## See also

- [train](../train.md)
- [predict](../predict.md)
- [classify](../classify.md)
- [split-fastq](../split-fastq.md)
- [match](../match.md)
- [Home](../index.md)
