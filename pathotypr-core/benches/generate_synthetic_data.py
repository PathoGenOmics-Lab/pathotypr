#!/usr/bin/env python3
"""
Generate synthetic test data for pathotypr benchmarks.

Creates a self-contained TEST_SYNTHETIC/ directory with:
  - reference FASTA (~4.4 MB, TB genome-sized)
  - sample FASTA (with mutations)
  - markers TSV (3700+ markers)
  - GFF3 annotation
  - compressed FASTQ (paired-end reads from sample)
  - trained model bundle (.bin.gz)

Usage:
    python generate_synthetic_data.py [--outdir TEST_SYNTHETIC]
"""

import argparse
import gzip
import hashlib
import os
import random
import struct
import sys
from pathlib import Path

GENOME_SIZE = 4_411_532  # ~TB genome
NUM_MARKERS = 3700
NUM_GENES = 4000
NUM_READS = 500_000  # 500k paired reads ≈ 30x coverage
READ_LEN = 150
KMER_SIZE = 21


def random_seq(length, rng):
    """Generate random DNA sequence."""
    return "".join(rng.choices("ACGT", k=length))


def mutate_base(base, rng):
    """Return a different base."""
    others = [b for b in "ACGT" if b != base]
    return rng.choice(others)


def generate_reference(outdir, rng):
    """Generate a random reference genome."""
    print("  Generating reference genome...")
    seq = random_seq(GENOME_SIZE, rng)
    path = outdir / "synthetic_reference.fasta"
    with open(path, "w") as f:
        f.write(">synthetic_reference\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
    return seq, path


def generate_markers(outdir, ref_seq, rng):
    """Generate marker positions with SNPs."""
    print("  Generating markers TSV...")
    # Pick random positions (avoiding edges for k-mer context)
    margin = 50
    positions = sorted(
        rng.sample(range(margin, len(ref_seq) - margin), NUM_MARKERS)
    )

    lineages = ["L1", "L2", "L3", "L4", "L5", "L6", "L7"]
    sub_lineages = {
        "L1": ["L1.1", "L1.2"],
        "L2": ["L2.1", "L2.2"],
        "L3": ["L3.1"],
        "L4": ["L4.1", "L4.2", "L4.3", "L4.9"],
        "L5": ["L5.1"],
        "L6": ["L6.1"],
        "L7": ["L7.1"],
    }

    path = outdir / "synthetic_markers.tsv"
    markers = []
    with open(path, "w") as f:
        f.write("position\tref_base\talt_base\tlineage1\tlineage2\tlineage3\tlineage4\tlineage5\n")
        for pos in positions:
            ref_base = ref_seq[pos]
            if ref_base not in "ACGT":
                ref_base = "A"
            alt_base = mutate_base(ref_base, rng)
            lin = rng.choice(lineages)
            sub = rng.choice(sub_lineages[lin])
            f.write(f"{pos}\t{ref_base}\t{alt_base}\t{lin}\t{sub}\t\t\t\n")
            markers.append((pos, ref_base, alt_base, lin))

    return markers, path


def generate_sample(outdir, ref_seq, markers, rng, lineage="L4"):
    """Generate a sample genome with mutations for a specific lineage."""
    print(f"  Generating sample genome (lineage {lineage})...")
    seq = list(ref_seq)

    # Apply mutations for matching lineage markers
    applied = 0
    for pos, ref_base, alt_base, lin in markers:
        if lin == lineage and rng.random() < 0.95:  # 95% of lineage markers mutated
            seq[pos] = alt_base
            applied += 1

    # Add some random noise mutations
    for _ in range(500):
        p = rng.randint(0, len(seq) - 1)
        seq[p] = mutate_base(seq[p], rng)

    sample_seq = "".join(seq)
    path = outdir / "synthetic_sample_L4.fasta"
    with open(path, "w") as f:
        f.write(">synthetic_sample_L4\n")
        for i in range(0, len(sample_seq), 80):
            f.write(sample_seq[i : i + 80] + "\n")

    return sample_seq, path


def generate_gff(outdir, ref_len, rng):
    """Generate synthetic GFF3 annotation."""
    print("  Generating GFF3 annotation...")
    path = outdir / "synthetic_annotation.gff3"

    genes = []
    pos = 100
    for i in range(NUM_GENES):
        gene_len = rng.randint(300, 3000)
        if pos + gene_len >= ref_len - 100:
            break
        strand = rng.choice(["+", "-"])
        genes.append((pos, pos + gene_len, strand, f"gene_{i:04d}"))
        pos += gene_len + rng.randint(50, 500)

    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region synthetic_reference 1 {ref_len}\n")
        for start, end, strand, name in genes:
            f.write(
                f"synthetic_reference\tsynth\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={name};Name={name}\n"
            )
            f.write(
                f"synthetic_reference\tsynth\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds_{name};Parent={name}\n"
            )

    return path


def generate_fastq(outdir, sample_seq, rng):
    """Generate paired-end FASTQ reads from sample."""
    print(f"  Generating {NUM_READS} paired-end reads...")
    r1_path = outdir / "synthetic_reads_R1.fastq.gz"
    r2_path = outdir / "synthetic_reads_R2.fastq.gz"
    seq_len = len(sample_seq)

    complement = str.maketrans("ACGTacgt", "TGCAtgca")

    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i in range(NUM_READS):
            frag_len = rng.randint(200, 500)
            start = rng.randint(0, seq_len - frag_len - 1)

            # R1: forward
            r1_seq = sample_seq[start : start + READ_LEN]
            # R2: reverse complement of end of fragment
            r2_start = start + frag_len - READ_LEN
            r2_seq = sample_seq[r2_start : r2_start + READ_LEN]
            r2_seq = r2_seq[::-1].translate(complement)

            # Add sequencing errors (~0.1%)
            r1_seq = list(r1_seq)
            r2_seq = list(r2_seq)
            for seq_list in [r1_seq, r2_seq]:
                for j in range(len(seq_list)):
                    if rng.random() < 0.001:
                        seq_list[j] = mutate_base(seq_list[j], rng)
            r1_seq = "".join(r1_seq)
            r2_seq = "".join(r2_seq)

            qual = "I" * READ_LEN  # Q40 quality

            f1.write(f"@read_{i}/1\n{r1_seq}\n+\n{qual}\n")
            f2.write(f"@read_{i}/2\n{r2_seq}\n+\n{qual}\n")

    return r1_path, r2_path


def generate_training_fasta(outdir, ref_seq, markers, rng):
    """Generate multi-sequence FASTA for train/predict benchmarks."""
    print("  Generating training FASTA (20 sequences)...")
    path = outdir / "synthetic_training.fasta"

    lineages = ["L1", "L2", "L3", "L4"]
    with open(path, "w") as f:
        for i in range(20):
            lin = lineages[i % len(lineages)]
            seq = list(ref_seq)
            for pos, ref_base, alt_base, m_lin in markers:
                if m_lin == lin and rng.random() < 0.9:
                    seq[pos] = alt_base
            # Random noise
            for _ in range(200):
                p = rng.randint(0, len(seq) - 1)
                seq[p] = mutate_base(seq[p], rng)
            sample_seq = "".join(seq)
            f.write(f">sample_{i:03d}|{lin}\n")
            for j in range(0, len(sample_seq), 80):
                f.write(sample_seq[j : j + 80] + "\n")

    return path


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic benchmark data")
    parser.add_argument(
        "--outdir",
        default=str(Path(__file__).resolve().parent.parent.parent / "TEST_SYNTHETIC"),
        help="Output directory",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"Generating synthetic data in {outdir}/")
    print(f"Seed: {args.seed}")
    print()

    ref_seq, ref_path = generate_reference(outdir, rng)
    markers, markers_path = generate_markers(outdir, ref_seq, rng)
    sample_seq, sample_path = generate_sample(outdir, ref_seq, markers, rng)
    gff_path = generate_gff(outdir, len(ref_seq), rng)
    r1_path, r2_path = generate_fastq(outdir, sample_seq, rng)
    training_path = generate_training_fasta(outdir, ref_seq, markers, rng)

    print()
    print("Generated files:")
    for p in sorted(outdir.iterdir()):
        size_mb = p.stat().st_size / (1024 * 1024)
        print(f"  {p.name:40s}  {size_mb:8.1f} MB")

    print()
    print("Done! Use these with:")
    print(f"  TEST_DIR={outdir} cargo bench --bench benchmarks")


if __name__ == "__main__":
    main()
