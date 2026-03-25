# Input Format Reference

## Training FASTA (`train`)

Standard FASTA format. The **first whitespace-delimited token** in each header is used as the class label.

```
>L4 sample_0001 extra_info
ACTGATCGATCG...
>L2 sample_0002
ACTGATCGATCG...
```

Labels: `L4`, `L2`.

Requirements:
- At least 2 distinct classes
- At least a few sequences per class (10+ recommended for reliable training)

## Marker TSV (`classify` and `split-fastq`)

Tab-separated, no header required (lines starting with `#` are skipped).

### Columns

| Column | Content | Required |
|---|---|---|
| 1 | Genomic position (1-based) | ✅ |
| 2 | REF allele | ✅ |
| 3 | ALT allele | ✅ |
| 4+ | Lineage hierarchy (one level per column) | ✅ (at least one) |
| after first empty | Annotation columns (gene, mutation) | Optional |

### Example

```
#pos    ref    alt    level1    level2           gene         mutation
761155  C      T      L4        L4.9             gyrA         Ser95Thr
2155168 G      A      L2        L2.2             katG         Ser315Thr
4247431 CC     C      L1        L1.1             ethA         frameshift
```

### Notes

- Lineage columns are read left-to-right until the first empty cell
- Hierarchical nesting: `L4` → `L4;L4.9` internally (semicolon-joined path)
- `classify` works best with SNP markers (single-base REF and ALT)
- `split-fastq` handles SNPs, MNVs, and small indels. However, indels in FASTQ mode are intentionally skipped because short reads produce unreliable k-mer matches across repetitive regions

## Sample List TSV (`classify --input-list`)

```
sample_name    /absolute/path/to/sample.fasta    /optional/path/to/sample.gff3
```

- Column 1: sample name (must be unique)
- Column 2: path to FASTA file
- Column 3: optional path to GFF3 annotation

## Sample List TSV (`split-fastq --input-list` and `match --input-list`)

```
sample_name    /path/to/reads_R1.fastq.gz    /path/to/reads_R2.fastq.gz
```

- Column 1: sample name (must be unique)
- Column 2+: paths to FASTQ files (one for single-end, two for paired-end)

## Reference FASTA

For `classify` and `split-fastq`: a **single-record** FASTA file containing the reference genome used to define marker positions.

For `match`: a **multi-record** FASTA file containing all candidate reference genomes.

## GFF3 Annotation (`classify`)

Standard GFF3 format. Only `CDS` features are used. Gene names are extracted from attributes in this priority order:

1. `gene=...`
2. `locus_tag=...`
3. `Name=...`
4. `ID=...`

Coordinates are 1-based (GFF3 standard) and converted internally to 0-based for interval tree queries.
