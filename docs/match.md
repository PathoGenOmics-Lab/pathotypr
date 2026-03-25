# `pathotypr match`

Find the best matching reference genome for a set of FASTQ reads.

## How it works

1. **Count** k-mers from FASTQ reads (parallel, with optional noise filtering)
2. **Compare** against references using streaming mode: processes references in batches, constant memory
3. **Score** each reference by weighted k-mer containment fraction
4. **Report** the best match

## Usage

```bash
# Match reads against reference genomes
pathotypr match \
  -i reads_R1.fastq.gz reads_R2.fastq.gz \
  -r references.fasta \
  -o match.tsv --excel

# From sample list
pathotypr match \
  -l samples.tsv \
  -r references.fasta \
  -o match.tsv
```

## Options

| Flag | Default | Description |
|---|---|---|
| `-i, --input` | — | One or more FASTQ files |
| `-l, --input-list` | — | TSV: sample_name → FASTQ path(s) |
| `-r, --references` | — | Multi-FASTA with reference genomes |
| `-k, --kmer-size` | `31` | K-mer size |
| `-o, --output` | stdout | Output TSV path |
| `-t, --threads` | all cores | Number of CPU threads |
| `--min-kmer-count` | `2` | Discard k-mers with fewer occurrences (noise filter) |
| `--excel` | off | Also generate .xlsx |
| `--strict-percentages` | on | Legacy-compatible weighted scoring |
| `--early-stop-confidence` | `0` (off) | Stop when confidence exceeds threshold |
| `--early-stop-min-kmers` | `1,000,000` | Minimum k-mers before early stop can trigger |

## Output columns

| Column | Description |
|---|---|
| `Query_Files` | Comma-separated input FASTQ paths |
| `Best_Match_Reference` | Header of the best-scoring reference |
| `Shared_Kmer_Fraction` | Weighted containment score (0–1) |

## Technical details

- **Streaming batch size**: `num_threads` references per batch, balancing parallelism vs memory
- **Constant memory**: regardless of reference count (1 to 500+), memory usage stays bounded
- **Noise filtering**: k-mers appearing only once are likely sequencing errors; filtered when total unique k-mers > 100K

## Algorithm Details

For in-depth documentation of the underlying algorithms:

- [Reference Matching](algorithms/reference-matching.md) — K-mer containment scoring, noise filtering, streaming batch processing, adaptive batch sizing, weighted score calculation
