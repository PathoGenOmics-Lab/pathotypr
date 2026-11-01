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

2. **Bucketing**: The bucket index is `bitkmer.0 & mask` where `mask = num_buckets - 1`. needletail packs the first base into the high bits, so masking to 20 bits keeps the **last ten bases** of the canonical k-mer. This is not a hash: see [What the buckets really are](#what-the-buckets-really-are) below.

3. **Bucket count**: 2^20 = 1,048,576 buckets.

4. **Sparse representation**: Each sequence produces a `Vec<(usize, f32)>` of its non-zero buckets, sorted by bucket index. H37Rv at k=21 produces 4.41M k-mers that fall into 652,008 buckets, 62% of the space.

## What the buckets really are

Masking the 2-bit code keeps its last ten bases, so the buckets are the 4^10
ten-mers, whatever k was asked for. k only decides which strand's ten bases are
kept, through the choice of canonical k-mer:

| k | buckets occupied by H37Rv (one position in seven) |
|---|---|
| 11 | 27.8% |
| 15 | 28.2% |
| 21 | 28.2% |
| 31 | 28.3% |

Two k-mers that share their last ten bases always share a bucket:
`AAAAAAAAAAACCCGGGTTTA` and `ACGTACGTAGACCCGGGTTTA` both land in 355068. Across
H37Rv at k=21, 96.0% of k-mer occurrences fall in a bucket shared with another
position, 6.77 of them per occupied bucket on average.

**It costs no accuracy in what was measured.** Synthetic genomes built on the
MTB ancestor with the lineage panel's own SNPs, six per lineage to train and
four held out to predict:

| lineages | own SNPs each | these buckets | a full 64-bit hash |
|---|---|---|---|
| L1 to L7 and A4 | 41 to 294 | 32/32 | 32/32 |
| L4.1 to L4.9 | 6 to 8 | 36/36 | 36/36 |
| A4.1 to A4.8 (*M. bovis*) | 2 to 5 | 32/32 | 32/32 |

The full hash gave slightly more decisive votes (a mean margin of 0.61 against
0.50 on the *M. bovis* sublineages) at twice the training time, and it changes
every feature index, so it would invalidate every model already trained. With no
call to gain, the buckets stay as they are; what changes is that they are
described as what they are.

!!! warning "Reading the importance reports"
    A bucket is not one k-mer. The `kmers` column of `*.importance.tsv` lists
    every k-mer of the training data that ends in the same ten bases, and the
    one a split actually used is among them without being singled out. Read the
    list as candidates, and confirm a candidate on the genomes before treating
    it as a marker.

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
