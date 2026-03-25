# Assembly-Based Marker Classification

## Overview

The `classify` command calls known variant markers in assembled genomes (FASTA). Unlike the FASTQ `split-fastq` approach (which uses diagnostic k-mers on raw reads), this module works on full-length assemblies and supports SNPs, MNVs, **and indels**.

## Architecture

```
Marker TSV + Reference FASTA ──→ Generate marker k-mers ──→ MarkerIndex
                                                                │
Query FASTA(s) ──→ Scan each genome for markers ──────────────→ MarkerMatch results
                                                                │
                   Optional GFF ──→ Annotate variants ────────→ AA changes
                                                                │
                                                          Lineage classification
                                                          Masked FASTA output
                                                          Summary TSV + Excel
```

## Step 1: Marker K-mer Generation

For each marker in the TSV (position, REF, ALT):

1. Extract flanking sequence from the reference genome (minimum 10 bp each side, configurable via `--min-flank-bases`)
2. Build a diagnostic k-mer centered on the variant:
   - For SNPs/MNVs: k-mer contains the variant at the center with equal flanks
   - For indels: separate ref and alt k-mers with adjusted flank lengths to maintain k-mer size

3. Index k-mers for fast lookup:
   - **Encoded index** (`MarkerIndex::Encoded`): For k ≤ 32, k-mers are 2-bit encoded into `u64` keys using needletail's `BitNuclKmer`. Lookup is a single hash table query.
   - **Text index** (`MarkerIndex::Text`): For k > 32, k-mers are stored as strings. Fallback for unusually long k-mers.

### Index Structure

```rust
enum MarkerIndex {
    Encoded(FxHashMap<u64, MarkerKmerEntry>),   // k ≤ 32: 2-bit encoded
    Text(FxHashMap<String, MarkerKmerEntry>),    // k > 32: string keys
}
```

Each entry stores the marker metadata: position, alleles, lineage, gene name, and flanking information needed to distinguish REF from ALT matches.

## Step 2: Genome Scanning

For each query genome (parallelized across genomes with rayon):

1. Read the genome sequence with needletail
2. Extract all k-mers (canonical, same k as the marker k-mers)
3. Look up each k-mer in the MarkerIndex
4. For matches, determine whether the match supports REF or ALT allele by comparing the variant-spanning region

### Match Result

```rust
struct MarkerMatch {
    position: usize,       // Genomic position (1-based)
    ref_allele: String,
    alt_allele: String,
    found_allele: String,  // What was actually observed
    lineage: String,
    gene: Option<String>,
    mutation: Option<String>,
}
```

## Step 3: Annotation (Optional)

When GFF files are provided, each called variant is annotated:

1. **GFF parsing**: Extracts CDS features with strand, phase, and locus tags using an interval tree (`rust_lapper::Lapper`) for O(log n) position lookups
2. **Amino acid translation**: Determines the codon position, extracts the codon from the reference and query, translates both to amino acids
3. **Output**: `gene_name:REF_AA{position}ALT_AA` notation (e.g., `rpoB:S450L`)

Supports:
- Forward and reverse strand genes
- Multi-nucleotide variants spanning codon boundaries
- Synonymous vs. non-synonymous classification

## Step 4: Lineage Classification

Two modes:

### Simple (default)
- Each marker votes for its deepest lineage level
- The lineage with the most ALT-supporting markers wins
- Resolves nested lineages by walking the lineage hierarchy

### Nested (`--nested-classification`)
- Hierarchical classification: first resolves major lineage, then sub-lineage
- Only considers markers belonging to the parent lineage when resolving sub-lineages
- More robust for datasets with uneven marker distribution across lineages

## Step 5: Output

### Summary TSV (`{prefix}_summary.tsv`)
Per-genome: sample name, predicted lineage, marker counts (REF/ALT/missing), percentage resolved.

### Detailed TSV (`{prefix}_details.tsv`)
Per-marker-per-genome: position, alleles, found allele, lineage, gene, amino acid change.

### Masked FASTA (optional, `--output-masked-fasta`)
Query sequences with marker positions replaced by `N`. Useful for downstream phylogenetic analysis where marker sites should be excluded (e.g., removing sites under selection before building neutral phylogenies).

### Excel (optional, `--excel`)
Same data as TSV with conditional formatting for confidence metrics.

## Key Differences from FASTQ Mode (`split-fastq`)

| Aspect | `classify` (FASTA) | `split-fastq` (FASTQ) |
|---|---|---|
| **Input** | Assembled genomes | Raw sequencing reads |
| **Indels** | ✅ Supported | ❌ Skipped (false positives) |
| **Annotation** | ✅ AA changes via GFF | ❌ Not available |
| **Masked output** | ✅ Masked FASTA | ❌ Not applicable |
| **Speed** | Slower (full genome scan) | Faster (Bloom filter pre-check) |
| **Accuracy** | Higher (full alignment context) | May miss low-coverage markers |
| **Parallelism** | Per-genome (rayon) | Per-read batch (rayon) |

## Complexity

| Operation | Time | Space |
|---|---|---|
| Build marker index | O(M × k) where M = markers | O(M) |
| Scan one genome | O(L × log M) where L = genome length | O(L + M) |
| GFF annotation | O(V × log G) where V = variants, G = genes | O(G) |
| Total (N genomes) | O(N × L × log M) / num_threads | O(L + M + G) |
