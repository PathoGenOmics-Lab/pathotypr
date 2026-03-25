# Marker File Format

pathotypr uses tab-separated (TSV) marker files to define lineage-specific SNPs for classification.

## Format

```
position	ref	alt	lineage
615938	A	G	L1
1799921	C	G	L2
4398141	A	G	L3
```

| Column | Description | Required |
|---|---|---|
| `position` | 1-based position in the reference genome | ✅ |
| `ref` | Reference allele (A, C, G, or T) | ✅ |
| `alt` | Alternative allele (A, C, G, or T) | ✅ |
| `lineage` | Lineage or group label | ✅ |

## Rules

- **Tab-separated**, no spaces around delimiters
- **Header row** is required (first line must contain `position`)
- **1-based coordinates** relative to the reference genome
- **Single nucleotide** only — no indels, no multi-allelic
- Lineage labels can use any naming convention (e.g., `L1`, `L1.1.2`, `ST258`, `clade_A`)
- Hierarchical lineages use semicolons: `L4;L4.1;L4.1.2`
- Duplicate positions with different lineages are allowed (shared markers)
- File encoding: UTF-8

## Hierarchical Lineages

For sub-lineage resolution, use semicolon-separated hierarchies:

```
position	ref	alt	lineage
615938	A	G	L1
789012	C	T	L1;L1.1
890123	G	A	L1;L1.1;L1.1.2
```

Each sample is classified at the deepest level where it matches sufficient markers.

## Drug Resistance Markers

DR markers use a separate file with the same format but include gene annotation:

```
position	ref	alt	lineage	gene	mutation
761155	C	T	RIF	rpoB	S450L
1673425	C	T	INH	katG	S315T
```

| Column | Description | Required |
|---|---|---|
| `gene` | Gene name | Optional |
| `mutation` | Amino acid change | Optional |

## Reference Genome

Your markers must match a reference genome in FASTA format. Both files are loaded together:

```bash
# Classify with custom markers
pathotypr classify \
  --input genomes/ \
  --markers my_markers.tsv \
  --reference my_reference.fasta \
  --output results/
```

## Validation

Before proposing a marker set, test it locally:

```bash
# Test on a small set of genomes with known lineages
pathotypr classify \
  --input test_genomes/ \
  --markers my_markers.tsv \
  --reference my_reference.fasta \
  --output validation/

# Check concordance with expected lineages
cat validation/*_summary.tsv
```

### Quality Criteria

| Metric | Minimum | Recommended |
|---|---|---|
| Markers per lineage | ≥ 5 | ≥ 20 |
| Validation genomes | ≥ 10 | ≥ 50 per lineage |
| Concordance | ≥ 95% | ≥ 99% |
| False positive rate | < 5% | < 1% |

## Example Marker Sets

| Organism | Lineages | Markers | Reference |
|---|---|---|---|
| *M. tuberculosis* | L1–L10, A1–A4 | 3,707 | H37Rv (GCF_000195955.2) |

## Proposing New Markers

Use the [Marker Proposals](https://github.com/mycolega/pathotypr/discussions/categories/ideas) discussion forum to propose new marker sets. See the discussion template for required information.

> **Note:** We currently only accept marker proposals for the *Mycobacterium tuberculosis* complex (MTBC). While pathotypr is organism-agnostic, community marker curation is limited to MTBC at this time.
