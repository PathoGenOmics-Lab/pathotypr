# `pathotypr classify`

Marker-based variant calling on assembled genomes.

## How it works

1. **Read** reference genome and marker definitions (TSV)
2. **Generate** diagnostic k-mers for each marker (ALT allele flanked by reference context)
3. **Scan** each query genome for marker k-mers using needletail's bit-packed k-mer iterator
4. **Annotate** matches with GFF gene information and amino acid changes (if GFF provided)
5. **Classify** each genome into a lineage using hierarchical path validation
6. **Export** detailed per-marker results and a summary with major lineage calls

## Usage

```bash
# Single genome
pathotypr classify \
  -m markers.tsv -r reference.fasta \
  -i sample.fasta --gff sample.gff3 \
  -o classify_run --excel

# Multiple genomes from a list
pathotypr classify \
  -m markers.tsv -r reference.fasta \
  -l genomes.tsv \
  -o classify_batch --nested-classification

# Multiple FASTA files (GUI batch mode)
pathotypr classify \
  -m markers.tsv -r reference.fasta \
  --input-files sample1.fasta sample2.fasta \
  -o batch_run
```

## Options

| Flag | Default | Description |
|---|---|---|
| `-m, --markers` | required | Marker definitions TSV |
| `-r, --reference` | required | Reference genome FASTA (single record) |
| `-i, --input` | — | Single query FASTA |
| `--input-files` | — | Multiple query FASTA files |
| `-l, --input-list` | — | TSV list of samples (name, FASTA path, optional GFF) |
| `--gff` | — | GFF3 annotation for `--input` mode |
| `--gff-files` | — | GFF3 files matched to `--input-files` by filename stem |
| `-o, --output-prefix` | required | Prefix for output files |
| `--kmer-size` | `31` | K-mer size for marker matching |
| `--min-flank-bases` | `10` | Minimum flanking bases on each side of allele |
| `--nested-classification` | off | Use hierarchical lineage path validation |
| `--output-masked-fasta` | off | Write masked FASTA (marker sites → N) |
| `-t, --threads` | all cores | Number of CPU threads |
| `--excel` | off | Also generate .xlsx files |

## Output files

| File | Content |
|---|---|
| `<prefix>.tsv` | Detailed: one row per marker match per genome |
| `<prefix>_summary.tsv` | Summary: one row per genome with major lineage |
| `<prefix>_*.xlsx` | Excel versions (if `--excel`) |
| `*_masked.fasta` | Masked genomes (if `--output-masked-fasta`) |

## Detailed output columns

`genome`, `k-mer`, `k-merPOS`, `SNPgenome`, `SNPreference`, `REF`, `ALT`, `lineage`, `Gene`, `Gene_Start`, `Gene_End`, `AA_Pos`, `AA_Change`

## Lineage classification logic

1. Find all lineages with a fully supported ancestor path (every parent node present)
2. Select the deepest valid candidate
3. If the most abundant lineage is from a different branch with more support, it wins
4. Ties broken by SNP count, then lexicographic order

Enable `--nested-classification` for hierarchical path validation; without it, the most abundant lineage wins directly.

## Technical details

- **Marker index**: u64-encoded k-mers for fast lookup (falls back to string index for non-ACGT bases)
- **GFF annotation**: CDS features parsed into an interval tree (rust-lapper) for O(log n) gene lookup
- **Amino acid translation**: handles both forward and reverse strand genes, multi-codon MNVs, and frameshift detection
- **Parallel processing**: genomes processed in chunks of 256 with rayon

## Algorithm Details

For in-depth documentation of the underlying algorithms:

- [Assembly Classification](algorithms/assembly-classification.md) — Marker k-mer generation, encoded vs text index, GFF annotation, indel support, lineage classification modes
