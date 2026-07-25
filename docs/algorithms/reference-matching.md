# Reference matching (k-mer containment)

**How pathotypr picks the closest reference genome for a read set using weighted k-mer containment scoring.**

The `match` command finds the closest reference genome for a set of FASTQ reads. It uses **k-mer containment scoring** — measuring what fraction of the query's k-mers are present in each reference genome.

```text
FASTQ reads ──→ Count query k-mers ──→ Filter noise ──→ Score references (streaming) ──→ Best match
```

## Algorithm

### Step 1: Query k-mer counting

All FASTQ files are processed in parallel (one rayon task per file):

1. Parse FASTQ with needletail (handles .gz transparently)
2. Extract canonical k-mers (k=31 by default) via `bit_kmers(k, canonical=true)`
3. Count occurrences in a thread-local `FxHashMap<u64, u32>`
4. Merge thread-local maps via parallel reduction

### Step 2: Noise filtering

Singleton k-mers (count=1) are likely sequencing errors:

- Only applied when total unique k-mers > 100,000
- Default threshold: `min_kmer_count = 2`
- Typical reduction: ~90-93% of unique k-mers removed
- Effect: faster scoring and fewer false matches

### Step 3: Streaming reference scoring

References are read from a multi-FASTA file in **adaptive batches** (batch_size = num_threads, clamped to 4–64):

```text
Multi-FASTA ──→ Read batch of N refs ──→ Score batch in parallel ──→ Update best
                                              │                         │
                                              └── Free batch memory ────┘
                                                                        │
                ──→ Read next batch ──────────────────────────────────────┘
```

!!! info "Why streaming?"
    With 500 reference genomes of ~4.4 Mb each, loading all references would require ~2.2 GB for sequences alone, plus ~17 GB for k-mer sets. Streaming processes `num_threads` references at a time, keeping memory bounded.

### Scoring a single reference

For each reference in the batch (parallelized with rayon):

1. Extract all canonical k-mers from the reference sequence
2. Sort and deduplicate (only unique k-mers matter)
3. For each unique reference k-mer, look it up in the query k-mer map
4. If found, add the **query count** to `shared_kmer_count`

### Score calculation

```text
Shared_Kmer_Fraction = Σ(query_count for shared k-mers) / Σ(all query k-mer counts)
```

This is a **weighted containment score**: k-mers that appear more frequently in the reads contribute more. A perfect match (all query k-mers present in the reference) gives a score of 1.0.

!!! tip "Weighted, not just presence/absence"
    Using weighted counts (instead of just presence/absence) makes the score robust to coverage variation — highly covered regions contribute proportionally.

## Memory model

| Component | Memory | Notes |
|---|---|---|
| Query k-mers | ~50–200 MB | All retained k-mers with counts |
| Reference batch | batch_size × ~4.4 MB | Raw sequences |
| Reference k-mers | batch_size × ~34 MB | Sorted unique k-mer arrays (temporary) |
| **Total** | ~300–500 MB typical | Independent of total number of references |

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--kmer-size` | 31 | Longer k = more specific matching |
| `--min-kmer-count` | 2 | Noise filter threshold |
| `--early-stop-confidence` | 0.0 | Reserved; not read in the current release (no effect) |
| `--early-stop-min-kmers` | 1,000,000 | Reserved; not read in the current release (no effect) |
| `--strict-percentages` | true | Reserved; scoring always uses weighted query counts |

## Output

```tsv
Query_Files	Best_Match_Reference	Shared_Kmer_Fraction
reads_R1.fastq.gz,...	GCF_000253355.1	0.5440
```

Only the single best match is reported. The score represents the fraction of query sequencing depth attributable to k-mers shared with that reference.

## Use cases

1. **Pre-mapping reference selection**: Find the closest reference before read alignment to maximize mapping quality
2. **Lineage verification**: Cross-check ML-predicted lineage against known reference genomes
3. **Contamination detection**: Low scores may indicate mixed samples or wrong species

## See also

- [match](../match.md)
