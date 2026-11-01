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
3. Build the matching reference k-mer from the same window: the same left flank, the REF allele, and as much right flank as keeps the length at k. This is how `split-fastq` builds its pair too.
4. Index both k-mers for fast lookup:
   - **Encoded index**: For k ≤ 32, marker k-mers are 2-bit encoded into `u64` keys with a hand-rolled `encode_kmer_2bit` (needletail's `bit_kmers` is used to generate the *query* k-mers during genome scanning). Lookup is a single hash table query.
   - **Text index**: Used when a marker k-mer contains a non-ACGT base that cannot be 2-bit encoded; k-mers are stored as strings. (The code also falls back for k > 32, but `validate_kmer_size` caps k at 31, so that branch is unreachable in practice.)

### Index structure

```rust
struct MarkerIndex {
    markers: Vec<MarkerKmerEntry>,     // sorted by position; a marker id is an index here
    alt: KmerLookup,                   // ALT k-mer -> marker id
    reference: KmerLookup,             // REF k-mer -> site id
    marker_site: Vec<Option<usize>>,   // the site of each marker, None without a REF k-mer
    site_count: usize,
}
```

Each entry stores the marker metadata: position, alleles, both k-mers, lineage, gene name, and the flank lengths needed to locate the variant bases within a matching k-mer.

Markers whose reference k-mers are identical form a **site**: same position, same reference allele, same window. Two markers with the same ALT k-mer cannot both be indexed; the later one in the TSV replaces the earlier, with a warning.

## Step 2: Genome scanning

For each query genome (parallelized across genomes with rayon):

1. Read every record (contig) of the genome with needletail
2. Scan each contig on **both strands** — the forward sequence and its reverse complement — because contig orientation in a draft assembly is arbitrary
3. Extract every k-mer of the same size k used for the marker k-mers and look it up in both lookups in a single pass; a marker seen on several contigs, or on both strands, counts once
4. An ALT hit means the genome carries the variant allele, and the entry's flank lengths then locate the exact variant-spanning bases. A REF hit means some copy of the window carries the reference allele
5. Once every contig has been scanned, each marker gets one of four states:

A k-mer says *what* it matched but not *where*, so a paralogous copy of either
allele reads as the marker's own locus. An anchor is only used when it occurs
once in the genome and once in the reference, so it does say where. Wherever an anchor can read the
locus, it decides:

| The anchors read | State |
|---|---|
| the marker's ALT | present |
| anything else | absent: the locus was read and does not carry it |
| nothing usable | the k-mers decide, below |

And when no anchor can read it:

| ALT k-mer | REF k-mer of its site | State |
|---|---|---|
| found, and the reference does not carry it | not found | present |
| found, but the reference carries it too | anything | undetermined, `alt_kmer_in_reference` |
| found | found | undetermined, `ambiguous_both_alleles` |
| not found | found as often as the reference has it | absent |
| not found | found, but the reference has more copies | undetermined, `ref_kmer_repeated` |
| not found | nothing of the site found | undetermined, `not_found` |
| anything | the marker has no reference k-mer at all | undetermined, `no_reference_kmer` |

Two anchors that both read the locus and disagree leave the marker undetermined
as `anchors_disagree`, whatever the k-mers found.

A k-mer is an exact match, so another change anywhere in the window (a
synonymous substitution, a lineage polymorphism) hides both alleles. Before the
anchors, that case was reported as absence.

### Match result

```rust
struct MarkerMatch {
    marker_id: usize,        // position of the marker in MarkerIndex::markers
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
2. **Amino acid translation**: Determines the codon position and translates the reference codon against the **sample's own codon**. A marker usually describes one base of three, so building the alternative codon from the reference and that one allele asserts the other two are reference. When the sample carries a second change in the same codon that assertion is false, and the amino acid reported is one nobody has: a sample reading TTC at rpoB 435 used to be reported as Tyr by the marker for its first base and as Val by the marker for its second. The bases around the variant are carried from wherever the marker was matched, so the codon is read rather than assumed
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
Per genome, five columns: `genome`, `lineage:count` (marker tallies at **every** hierarchical level, each ancestor path prefix counted separately, space-separated), `major_lineage` (the final call), `undetermined` (how many markers could be read neither way) and `undetermined:count` (those markers per first level of their path).

### Detailed TSV (`{prefix}.tsv`)
One row per matched marker per genome, with columns: `genome`, `k-mer`, `k-merPOS`, `SNPgenome`, `SNPreference`, `REF`, `ALT`, `lineage`, `Gene`, `Gene_Start`, `Gene_End`, `AA_Pos`, `AA_Change`, `alleles_found` (`ALT`, or `ALT,REF` when both k-mers turned up and an anchor settled it) and `evidence` (`kmer`, or the anchor that read the allele).

### Undetermined TSV (`{prefix}_undetermined.tsv`)
One row per undetermined marker per genome: `genome`, `SNPreference`, `REF`, `ALT`, `lineage`, `Gene`, `mutation`, `reason`.

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
