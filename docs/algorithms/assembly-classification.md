# Assembly-based marker classification

**A walkthrough of how `classify` matches known variant markers in assembled genomes: diagnostic k-mer generation, genome scanning, optional GFF annotation, lineage classification, and output.**

The `classify` command calls known variant markers directly in assembled genomes (FASTA), supporting SNPs, MNVs, and indels.

!!! info "Assemblies, not reads"
    Unlike the FASTQ `split-fastq` path, which screens raw reads with diagnostic k-mers, this module scans full-length assemblies — so it can call indels and, with a GFF, report amino-acid changes.

## Architecture

```text
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

## Step 1: Marker k-mer generation

For each marker in the TSV (position, REF, ALT):

1. Extract flanking sequence from the reference genome (minimum 10 bp each side, configurable via `--min-flank-bases`)
2. Build a diagnostic k-mer centered on the variant:
   - For SNPs/MNVs: the k-mer holds the ALT allele at the center with equal flanks
   - For indels: the ALT allele is placed in the k-mer and the flank lengths are adjusted so the k-mer keeps its length of k

3. Index k-mers for fast lookup:
   - **Encoded index** (`MarkerIndex::Encoded`): For k ≤ 32, marker k-mers are 2-bit encoded into `u64` keys with a hand-rolled `encode_kmer_2bit` (needletail's `bit_kmers` is used to generate the *query* k-mers during genome scanning). Lookup is a single hash table query.
   - **Text index** (`MarkerIndex::Text`): Used when a marker k-mer contains a non-ACGT base that cannot be 2-bit encoded; k-mers are stored as strings. (The code also falls back for k > 32, but `validate_kmer_size` caps k at 31, so that branch is unreachable in practice.)

### Index structure

```rust
enum MarkerIndex {
    Encoded(FxHashMap<u64, MarkerKmerEntry>),   // k ≤ 32: 2-bit encoded
    Text(FxHashMap<String, MarkerKmerEntry>),    // k > 32 or non-ACGT: string keys
}
```

Each entry stores the marker metadata: position, alleles, lineage, gene name, and the flank lengths needed to locate the variant bases within a matching k-mer.

## Step 2: Genome scanning

For each query genome (parallelized across genomes with rayon):

1. Read every record (contig) of the genome with needletail
2. Scan each contig on **both strands** — the forward sequence and its reverse complement — because contig orientation in a draft assembly is arbitrary
3. Extract every k-mer of the same size k used for the marker k-mers and look it up in the `MarkerIndex`; a marker seen on several contigs, or on both strands, counts once
4. Because only the ALT (variant) k-mer is indexed, any hit means the genome carries the variant allele; the entry's flank lengths then locate the exact variant-spanning bases

### Match result

```rust
struct MarkerMatch {
    kmer: String,            // the diagnostic (ALT) k-mer that matched
    genome_position: usize,  // 0-based start of the match in the genome
    ref_position: usize,     // 1-based marker position in the reference
    lineage: String,
    ref_allele_len: usize,
    alt_allele_len: usize,
    left_flank_len: usize,   // flank bases before the variant in the k-mer
    ref_allele: String,
    alt_allele: String,
    gene: Option<String>,
    mutation: Option<String>,
}
```

## Step 3: Annotation (optional)

When GFF files are provided, each called variant is annotated:

1. **GFF parsing**: Reads `CDS` features and stores each one's strand and gene name (taken from the `gene`, `locus_tag`, `Name`, or `ID` attribute) in an interval tree (`rust_lapper::Lapper`) for O(log n) position lookups
2. **Amino acid translation**: Determines the codon position, extracts the codon from the reference and query, translates both to amino acids
3. **Output**: for each affected codon the change is written as `REF_AA{position}ALT_AA(REF_codon>ALT_codon)` using three-letter amino-acid codes (e.g. `Ser450Leu(TCG>TTG)`); synonymous changes get a trailing ` syn`, an MNV spanning several codons joins its per-codon changes with `;`, and a length-changing (indel) variant is reported as `frameshift_pos{N}`. The gene name and first affected codon are written to the separate `Gene` and `AA_Pos` columns.

Supports:

- Forward and reverse strand genes
- Multi-nucleotide variants spanning codon boundaries
- Synonymous vs. non-synonymous classification

## Step 4: Lineage classification

Two modes:

### Simple (default)
- Each matched marker increments the count of its full lineage path **and of every ancestor prefix** (`L4`, `L4;L4.9`, `L4;L4.9;L4.9.1`)
- The path with the highest count wins, which is normally the shallowest level
- Ties are broken deterministically by lineage name

### Nested (`--nested-classification`)
- Hierarchical classification: first resolves major lineage, then sub-lineage
- Only considers markers belonging to the parent lineage when resolving sub-lineages

!!! note "Nested vs. simple"
    Nested classification is more robust for datasets with an uneven marker distribution across lineages.

## Step 5: Output

### Summary TSV (`{prefix}_summary.tsv`)
Per genome, three columns: `genome`, `lineage:count` (marker tallies at **every** hierarchical level — each ancestor path prefix is counted separately, space-separated), and `major_lineage` (the final call).

### Detailed TSV (`{prefix}.tsv`)
One row per matched marker per genome, with columns: `genome`, `k-mer`, `k-merPOS`, `SNPgenome`, `SNPreference`, `REF`, `ALT`, `lineage`, `Gene`, `Gene_Start`, `Gene_End`, `AA_Pos`, `AA_Change`.

### Masked FASTA (optional, `--output-masked-fasta`)
Query sequences with marker positions replaced by `N`, written as `<input-stem>_masked.fasta`. Only produced for sequences in reference coordinates (records whose length matches the reference).

!!! tip "Masked output"
    Useful for downstream phylogenetic analysis where marker sites should be excluded — for example, removing sites under selection before building neutral phylogenies.

### Excel (optional, `--excel`)
Same data as the TSV outputs, with a frozen header row and an autofilter. No conditional formatting is applied — classify has no confidence column.

## Key differences from FASTQ mode (`split-fastq`)

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

## See also

- [classify](../classify.md)
- [Marker format](../marker_format.md)
- [Input formats](../input-formats.md)
