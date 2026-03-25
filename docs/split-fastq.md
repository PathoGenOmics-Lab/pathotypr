# `pathotypr split-fastq`

Alignment-free genotyping directly from FASTQ reads.

## How it works

1. **Build** diagnostic k-mer pairs (REF + ALT) from marker definitions against the reference
2. **Index** all marker k-mers (forward + reverse complement) in a hash map with Bloom filter pre-check
3. **Scan** FASTQ reads in parallel batches of 8,192
4. **Count** REF vs ALT k-mer hits per marker
5. **Call** variants where ALT fraction exceeds threshold at sufficient depth
6. **Classify** into lineage using hierarchical path validation

## Usage

```bash
# Paired-end reads
pathotypr split-fastq \
  -m markers.tsv -r reference.fasta \
  -i reads_R1.fastq.gz -i reads_R2.fastq.gz \
  --paired -o genotype --excel

# Sample list
pathotypr split-fastq \
  -m markers.tsv -r reference.fasta \
  -l samples.tsv \
  -o batch_run --nested-classification

# Auto-detect pairing
pathotypr split-fastq \
  -m markers.tsv -r reference.fasta \
  -i sample_R1.fq.gz -i sample_R2.fq.gz \
  -o auto_run
```

## Options

| Flag | Default | Description |
|---|---|---|
| `-m, --markers` | required | Marker definitions TSV |
| `-r, --reference` | required | Reference genome FASTA (single record) |
| `-i, --input` | — | One or more FASTQ files |
| `-l, --input-list` | — | TSV: sample_name → FASTQ path(s) |
| `--paired` | off | Manual paired-end mode (files in pairs) |
| `--no-auto-paired` | off | Disable automatic paired-end detection |
| `--min-depth` | `10` | Minimum REF+ALT coverage to call a variant |
| `--min-alt-percent` | `95` | Minimum ALT fraction (%) to call a variant |
| `-k, --kmer-size` | `31` | K-mer size for diagnostic k-mers |
| `-o, --output-prefix` | `split` | Prefix for output files |
| `--nested-classification` | off | Use hierarchical lineage path validation |
| `-t, --threads` | all cores | Number of CPU threads |
| `--excel` | off | Also generate .xlsx files |

## Output files

| File | Content |
|---|---|
| `<prefix>_<sample>_mutations.tsv` | Per-marker detail: pos, alleles, counts, fraction, lineage |
| `<prefix>_summary.tsv` | One row per sample: lineage counts + major lineage call |
| `*.xlsx` | Excel versions (if `--excel`) |

## Why no indels?

Short reads produce unreliable k-mer matches for insertions and deletions, especially across repetitive regions (PE/PPE families in *M. tuberculosis*, IS elements). Indels are detected reliably via the assembly-based `classify` workflow.

## Technical details

- **Bloom filter**: 1M-bit filter rejects ~99% of non-matching k-mers before hash table lookup
- **Parallel scanning**: reads are buffered in a contiguous `ReadBatch` (zero per-read allocations) and processed in rayon chunks
- **Marker caching**: marker database is cached in memory across samples (fingerprinted by file path + size + mtime)
- **Paired-end auto-detection**: recognizes `_R1/_R2`, `_1/_2`, `.1/.2` naming patterns
- **Canonical k-mers**: the index stores both forward and RC k-mers explicitly, so scanning uses non-canonical mode for correctness

## Algorithm Details

For in-depth documentation of the underlying algorithms:

- [Marker Genotyping](algorithms/marker-genotyping.md) — Diagnostic k-mer construction, Bloom filter pre-check, ReadBatch buffer, parallel scanning, why indels are skipped in FASTQ mode
