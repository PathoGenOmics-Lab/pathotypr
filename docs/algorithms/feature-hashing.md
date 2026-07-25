# Feature hashing (the hashing trick)

**pathotypr converts each genome into a fixed-length numerical vector by hashing its k-mers into a fixed set of buckets — no explicit vocabulary required.**

This page explains why a vocabulary-based encoding does not scale, how the hashing trick replaces it, and how important buckets can be traced back to the k-mers and genomic coordinates they came from.

## Problem

Representing genomic sequences as numerical vectors for machine learning requires counting k-mers — short subsequences of length *k*. A naive vocabulary-based approach needs to enumerate all possible k-mers (4^k = 4 billion for k=16) or build a dictionary from the data (dataset-dependent, expensive to store).

## Solution: feature hashing

pathotypr uses the **hashing trick** (Weinberger et al., 2009): each k-mer is hashed into a fixed set of buckets using a bitmask, bypassing the need for an explicit vocabulary.

```text
Sequence: ACTGACTGCTA...
    ↓ extract k-mers (k=21)
    ACTGACTGCTA... → hash → bucket_idx = hash & 0xFFFFF
    ↓
Sparse vector: [(bucket_42, 3.0), (bucket_1087, 1.0), ...]
```

### Implementation

1. **K-mer extraction**: Uses [needletail](https://github.com/onecodex/needletail)'s `bit_kmers()` which encodes each k-mer as a 2-bit packed `u64`. This handles canonical k-mers (choosing the lexicographically smaller of forward and reverse complement) automatically.

2. **Hashing**: The bucket index is simply `bitkmer.0 & mask` where `mask = num_buckets - 1`.

3. **Bucket count**: 2^20 = 1,048,576 buckets.

4. **Sparse representation**: Each sequence produces a `Vec<(usize, f32)>` — only non-zero buckets are stored, sorted by bucket index. This is critical for efficiency: a typical 4.4 Mb genome at k=21 produces ~4.4M k-mers but only ~800K–1M unique bucket entries.

!!! info "Why no extra hash function is needed"
    Since needletail's 2-bit encoding already acts as a hash function (deterministic, uniform for biological sequences), no additional hash function is needed. The power-of-two constraint enables bitmask indexing (`& mask` instead of `% modulus`).

### Why this works

- **Stateless**: The same hash function is applied at train and predict time. No vocabulary to store or update.
- **Fixed dimensionality**: All genomes map to the same 1M-dimensional space, regardless of genome size or k-mer diversity.

!!! note "Collisions are acceptable"
    Multiple k-mers can map to the same bucket. With 1M buckets and ~4M k-mers per genome, collisions occur but empirically don't degrade classification accuracy because the Random Forest learns on aggregate count patterns, not individual k-mers.

### Reverse mapping

For feature importance analysis, pathotypr can reverse-map important buckets back to the original k-mers and their genomic coordinates:

1. **`reverse_map_buckets()`**: Re-scans all training sequences, collects k-mers that hash into important buckets. Used for the importance TSV.
2. **`reverse_map_with_coords()`**: Same but also records (sequence_index, 1-based position) for each hit. Used for the genomic coordinates TSV.

## Complexity

| Operation | Time | Space |
|---|---|---|
| Vectorize one sequence | O(L) where L = sequence length | O(min(4^k, B)) where B = buckets |
| Vectorize N sequences | O(N × L), parallelized with rayon | O(N × nnz) sparse |
| Reverse map | O(N × L × \|targets\|) | O(\|targets\| × hits) |

## Parameters

| Parameter | Value | Rationale |
|---|---|---|
| `num_buckets` | 2^20 = 1,048,576 | Balances collision rate vs. model size. Bucket count >> number of classes ensures discriminative power |
| Default k (train) | 21 | Long enough for specificity in bacterial genomes (4^21 ≈ 4.4 × 10^12 >> genome size), short enough that SNPs create distinguishable k-mer patterns |
| Default k (genotyping) | 31 | Longer k for higher specificity in marker-based genotyping where false positives are costly |

## See also

- [train](../train.md)
- [predict](../predict.md)
