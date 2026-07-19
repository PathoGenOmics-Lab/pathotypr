# Marker-based genotyping from FASTQ reads

**How pathotypr calls lineage-defining SNP markers straight from FASTQ reads: it builds diagnostic k-mers for each marker, screens read k-mers through a Bloom filter and hash index, and aggregates reference/alternate evidence per marker before classifying lineages.**

The `split-fastq` command performs alignment-free genotyping directly from raw FASTQ reads. Given a set of known SNP markers and a reference genome, it builds diagnostic k-mers for each marker and scans reads for evidence of reference vs. alternate alleles.

## Architecture

```text
Marker TSV + Reference FASTA ──→ Build Diagnostic K-mers ──→ Build Index + Bloom Filter
                                                                       │
FASTQ Reads ──→ Stream in batches of 8192 ──→ 2-bit encode k-mers ──→ Bloom check ──→ Index lookup
                                                                                           │
                                                                                    Aggregate counts
                                                                                    [ref_count, alt_count]
                                                                                           │
                                                                                    Lineage Classification
```

## Step 1: Diagnostic k-mer construction

For each SNP marker at position P with REF allele and ALT allele:

1. Extract flanking sequence from the reference genome
2. Build two k-mers of length 31 (default):
   - **REF k-mer**: left_flank + REF + right_flank
   - **ALT k-mer**: left_flank + ALT + right_flank

```text
Reference:  ...ACTGACTG[C]TAGCTAGC...
                        ↑ SNP position

REF k-mer:  ACTGACTG[C]TAGCTAGCGATCGATCGATCG   (k=31)
ALT k-mer:  ACTGACTG[T]TAGCTAGCGATCGATCGATCG   (k=31)
```

The flanking region is centered on the variant, with separate right-flank lengths for ref and alt k-mers when allele lengths differ (MNVs).

!!! warning "Indels are intentionally skipped"

    Indels (insertions/deletions) are **intentionally skipped** in the FASTQ workflow:

    - Short reads spanning repetitive regions (e.g., PE/PPE genes in *M. tuberculosis*) produce k-mers that match indel signatures at incorrect loci
    - This generates false positive variant calls
    - Indels are reliably detected only via the FASTA assembly `classify` workflow, which uses full-length alignment context

## Step 2: Hash index + Bloom filter

### Hash index

All diagnostic k-mers are hashed (2-bit encoding, forward + reverse complement) into an `FxHashMap<u64, MarkerHitList>`:

```rust
enum MarkerHitList {
    Single((marker_id, is_alt)),        // Most common: one marker per k-mer
    Multi(Vec<(marker_id, is_alt)>),    // Rare: k-mer collision between markers
}
```

Each marker contributes 4 entries: REF forward, REF revcomp, ALT forward, ALT revcomp.

### Bloom filter

A 128 KB Bloom filter (2^20 = 1,048,576 bits) acts as a fast pre-check:

- **2 hash functions** using multiplicative hashing: `key × c₁ mod BITS` and `key × c₂ mod BITS`
- **False positive rate**: ~0.004% with ~32K entries and 2 hash functions
- **Fits in L2 cache** (128 KB vs 512 KB typical L2/core)

!!! note "Why a Bloom filter?"

    The Bloom filter rejects ~99.9% of non-marker k-mers in ~2 nanoseconds, before they reach the hash table lookup.

## Step 3: FASTQ scanning

### Batched parallel processing

```text
FASTQ file ──→ needletail parser ──→ ReadBatch (contiguous buffer)
                                          │
                                    ┌─────┴─────┐
                                    │  rayon     │
                                    │  chunks    │
                                    └─────┬─────┘
                                          │
                                    Thread-local counts
                                          │
                                    ┌─────┴─────┐
                                    │  reduce    │
                                    └─────┬─────┘
                                          │
                                    Global counts
```

1. **ReadBatch**: A contiguous byte buffer that stores reads back-to-back, avoiding per-read `Vec<u8>` allocations. Memory is reused across batches via `clear()`.

2. **Batch size**: 8,192 reads per batch. Chosen to amortize rayon scheduling overhead while keeping memory bounded.

3. **Chunk size**: Each rayon chunk processes `max(64, batch_size / (2 × num_threads))` reads with thread-local count arrays.

4. **Per-read scanning**:
   - Extract all k-mers (k=31) via `bit_kmers()` — **no canonical normalization** because the index already contains both forward and reverse complement entries
   - Bloom filter check: `bloom_may_contain(bloom, kmer.0)` — rejects ~99.9% of k-mers
   - Hash table lookup: if in Bloom, check exact match in `FxHashMap`
   - Increment `counts[marker_id][is_alt]`

5. **Reduction**: Thread-local counts are summed into global counts via rayon's `fold + reduce`.

### Cancellation

- Checked every 8,192 records (atomic flag, near-zero overhead)
- `AtomicBool` shared across threads for fast stop
- Parse errors captured via `Arc<Mutex<Option<String>>>` with poison recovery

## Step 4: Lineage classification

The `classify_split_fastq` module takes the raw `[ref_count, alt_count]` per marker and determines lineages:

1. For each marker, determine the called allele (REF or ALT) based on count ratios
2. Walk the hierarchical lineage columns (e.g., L4 → L4.9 → L4.9.1)
3. Report the deepest resolved lineage

### Marker TSV format

```tsv
#pos    ref    alt    level1    level2    level3    <empty>    annotation
761155  C      T      L4        L4.9      L4.9.1               gene_name
```

!!! info "Reading the lineage columns"

    Lineage columns are read until the first empty cell. Columns after the empty cell are treated as annotations (gene name, drug resistance, etc.).

## Performance

### Why pathotypr is ~2× faster than fastlin

| Factor | pathotypr | fastlin |
|---|---|---|
| **K-mer rejection** | Bloom filter (2 ns) | Direct hash lookup |
| **Read parallelism** | rayon batch 8192, chunks per thread | Per-file parallelism |
| **Gzip decompression** | zlib-ng (SIMD) | miniz_oxide (pure Rust) |
| **Memory allocations** | Contiguous ReadBatch buffer | Per-read Vec<u8> |

### Benchmarks (real MTB data)

| Sample | pathotypr | fastlin | Speedup |
|---|---|---|---|
| ERR551304 (L2, 158 MB PE) | 1.49 s | 3.23 s | 2.2× |
| ERR552797 (A4, 264 MB PE) | 2.30 s | 4.78 s | 2.1× |

### Memory trade-off

pathotypr uses ~28 MB vs fastlin's ~3 MB. The difference comes from:

- Full reference genome in memory (4.4 MB for MTB)
- Larger binary (needletail, rayon, indicatif dependencies)
- mimalloc thread arenas (~8 MB)

## See also

- [split-fastq](../split-fastq.md)
- [Marker format](../marker_format.md)
