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

## Where the code lives

Each document above describes one part of this tree, so the two read together.

```text
pathotypr-core/src/
├── main.rs                   # CLI entry point
├── lib.rs                    # Library root
├── defaults.rs               # Default resource URLs and filenames
├── train.rs                  # Random Forest training + OOB + CV
├── predict.rs                # Streaming batch prediction
├── classify/                 # Assembly-based marker classification
│   ├── mod.rs                #   Orchestration + genome analysis
│   ├── markers.rs            #   Marker parsing + k-mer generation
│   ├── annotation.rs         #   GFF parsing + AA translation
│   └── masking.rs            #   FASTA masking at marker sites
├── classify_split_fastq.rs   # FASTQ genotyping orchestration
├── split_kmer.rs             # Diagnostic k-mer engine + Bloom filter
├── match/mod.rs              # Reference matching: scoring + coarse-to-fine
├── sparse_tree.rs            # Custom CART on sparse vectors
├── vectorizer.rs             # Feature hashing (hashing trick)
├── model.rs                  # Model bundle + label encoder
├── lineage.rs                # Hierarchical lineage classification
├── fasta_io.rs               # FASTA reading (needletail)
├── paired_end.rs             # Paired-end FASTQ detection
├── excel.rs                  # Streaming Excel export
├── errors.rs                 # Error types + cancellation
└── common.rs                 # Thread pool + shared utilities
```

The desktop app lives alongside it in `src-tauri/` (Rust backend) and
`frontend/` (HTML/CSS/JS); see [Desktop GUI](../gui.md) for that side.

## See also

- [train](../train.md)
- [predict](../predict.md)
- [classify](../classify.md)
- [split-fastq](../split-fastq.md)
- [match](../match.md)
- [Home](../index.md)
