#!/usr/bin/env python3
"""
End-to-end scaling benchmark for pathotypr's 3 core modules.

Measures wall time AND peak RSS. Varies both input count and genome size.

Modules:
  1. classify      — FASTA assemblies → lineage
  2. split-fastq   — FASTQ reads → lineage
  3. match          — FASTQ → best reference

Scaling dimensions:
  A. Number of inputs (genomes / reads / references)
  B. Genome size (0.5 Mb .. 10 Mb)

Usage:
    python scaling_benchmark.py [--binary pathotypr] [--outdir results_scaling]
"""

import argparse
import gzip
import os
import platform
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Style — clean, readable
# ---------------------------------------------------------------------------
MM_TO_IN = 1 / 25.4
SINGLE_COL = 89 * MM_TO_IN
DOUBLE_COL = 183 * MM_TO_IN

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "DejaVu Sans", "Arial"],
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7,
    "legend.title_fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
DEFAULT_GENOME_SIZE = 4_411_532
READ_LEN = 150
SEED = 42
REPEATS = 3

# Scaling by count (fixed genome = 4.4 Mb)
CLASSIFY_N = [1, 2, 5, 10, 20, 50]
SPLITFASTQ_N = [50_000, 100_000, 250_000, 500_000, 1_000_000]
MATCH_N = [1, 2, 5, 10, 20]

# Scaling by genome size (fixed count: 5 genomes / 250k reads / 5 refs)
GENOME_SIZES = [500_000, 1_000_000, 2_000_000, 4_000_000, 6_000_000, 10_000_000]
GENOME_SIZE_LABELS = ["0.5", "1", "2", "4", "6", "10"]  # Mb

MODULE_COLORS = {
    "classify": "#E07A5F",
    "split-fastq": "#AF64D1",
    "match": "#0077B6",
}


# ---------------------------------------------------------------------------
# Time + RSS measurement
# ---------------------------------------------------------------------------
def detect_time_cmd():
    """Find a time command that reports peak RSS."""
    if platform.system() == "Darwin":
        # macOS: /usr/bin/time -l reports RSS in bytes
        if shutil.which("/usr/bin/time"):
            return "/usr/bin/time", "macos"
    if shutil.which("gtime"):
        return "gtime", "gnu"
    if shutil.which("/usr/bin/time"):
        return "/usr/bin/time", "gnu"
    return None, None


def parse_time_output(stderr_text, flavor):
    """Extract peak RSS in MB from time stderr output."""
    if not stderr_text:
        return None

    # macOS format: "  6930432  maximum resident set size" (number BEFORE text, bytes)
    # GNU format:   "Maximum resident set size (kbytes): 12345" (number AFTER text, KB)
    m = re.search(r"(\d+)\s+maximum resident set size", stderr_text, re.IGNORECASE)
    if not m:
        # GNU fallback
        m = re.search(r"Maximum resident set size[^0-9]*(\d+)", stderr_text)
    if not m:
        return None

    raw = int(m.group(1))
    if flavor == "macos":
        return raw / (1024 * 1024)  # bytes → MB
    else:
        return raw / 1024  # KB → MB


TIME_CMD, TIME_FLAVOR = detect_time_cmd()


def run_timed(cmd, timeout=600):
    """Run command with /usr/bin/time, return (wall_seconds, peak_rss_mb)."""
    import tempfile as _tf

    time_path = _tf.mktemp(suffix=".time")

    if TIME_CMD:
        if TIME_FLAVOR == "macos":
            time_flag = "-l"
        else:
            time_flag = "-v"
        # /usr/bin/time writes to stderr. We redirect ALL stderr (time + cmd)
        # to the time file, stdout goes to /dev/null
        escaped = " ".join(f'"{c}"' if " " in c else c for c in cmd)
        shell_cmd = f'{TIME_CMD} {time_flag} {escaped} >/dev/null 2>"{time_path}"'
    else:
        escaped = " ".join(f'"{c}"' if " " in c else c for c in cmd)
        shell_cmd = f'{escaped} >/dev/null 2>/dev/null'

    start = time.perf_counter()
    rc = os.system(f'{{ {shell_cmd} ; }} 2>>"{time_path}"')
    elapsed = time.perf_counter() - start

    exit_code = os.waitstatus_to_exitstatus(rc) if hasattr(os, 'waitstatus_to_exitstatus') else rc >> 8

    if exit_code != 0:
        print(f"    FAILED (exit {exit_code}): {' '.join(cmd[:4])}...")
        try:
            with open(time_path) as f:
                print(f"    stderr: {f.read()[:300]}")
        except Exception:
            pass
        if os.path.exists(time_path):
            os.unlink(time_path)
        return None, None

    rss_mb = None
    if TIME_CMD and os.path.exists(time_path):
        try:
            with open(time_path) as f:
                time_output = f.read()
            rss_mb = parse_time_output(time_output, TIME_FLAVOR)
        except Exception:
            pass

    if os.path.exists(time_path):
        os.unlink(time_path)
    return elapsed, rss_mb


# ---------------------------------------------------------------------------
# Synthetic data generators
# ---------------------------------------------------------------------------
def random_seq(length, rng):
    return "".join(rng.choices("ACGT", k=length))


def mutate_base(base, rng):
    others = [b for b in "ACGT" if b != base]
    return rng.choice(others)


def make_reference(path, genome_size, rng):
    seq = random_seq(genome_size, rng)
    with open(path, "w") as f:
        f.write(">synthetic_reference\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    return seq


def make_markers(path, ref_seq, rng, n_markers=3700):
    margin = 50
    max_pos = len(ref_seq) - margin
    actual_n = min(n_markers, max_pos - margin)
    positions = sorted(rng.sample(range(margin, max_pos), actual_n))
    lineages = ["L1", "L2", "L3", "L4"]
    sub = {"L1": "L1.1", "L2": "L2.1", "L3": "L3.1", "L4": "L4.1"}

    markers = []
    with open(path, "w") as f:
        f.write("position\tref_base\talt_base\tlineage1\tlineage2\tlineage3\tlineage4\tlineage5\n")
        for pos in positions:
            rb = ref_seq[pos] if ref_seq[pos] in "ACGT" else "A"
            ab = mutate_base(rb, rng)
            lin = rng.choice(lineages)
            f.write(f"{pos}\t{rb}\t{ab}\t{lin}\t{sub[lin]}\t\t\t\n")
            markers.append((pos, rb, ab, lin))
    return markers


def make_multifasta(path, ref_seq, markers, n_genomes, rng):
    lineages = ["L1", "L2", "L3", "L4"]
    with open(path, "w") as f:
        for i in range(n_genomes):
            lin = lineages[i % len(lineages)]
            seq = list(ref_seq)
            for pos, rb, ab, ml in markers:
                if ml == lin and rng.random() < 0.9:
                    seq[pos] = ab
            for _ in range(200):
                p = rng.randint(0, len(seq) - 1)
                seq[p] = mutate_base(seq[p], rng)
            f.write(f">g{i:04d}|{lin}\n")
            f.write("".join(seq) + "\n")


def make_fastq(r1_path, r2_path, ref_seq, n_reads, rng):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    seq_len = len(ref_seq)
    qual = "I" * READ_LEN
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i in range(n_reads):
            frag_len = rng.randint(200, 500)
            start = rng.randint(0, seq_len - frag_len - 1)
            r1 = ref_seq[start:start+READ_LEN]
            r2s = start + frag_len - READ_LEN
            r2 = ref_seq[r2s:r2s+READ_LEN][::-1].translate(complement)
            r1, r2 = list(r1), list(r2)
            for s in [r1, r2]:
                for j in range(len(s)):
                    if rng.random() < 0.001:
                        s[j] = mutate_base(s[j], rng)
            f1.write(f"@r{i}/1\n{''.join(r1)}\n+\n{qual}\n")
            f2.write(f"@r{i}/2\n{''.join(r2)}\n+\n{qual}\n")


def make_ref_db(path, ref_seq, n_refs, rng):
    with open(path, "w") as f:
        for i in range(n_refs):
            seq = list(ref_seq)
            for _ in range(int(len(seq) * 0.01)):
                p = rng.randint(0, len(seq) - 1)
                seq[p] = mutate_base(seq[p], rng)
            f.write(f">ref_{i:04d}\n")
            f.write("".join(seq) + "\n")


# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------
def bench_classify_by_count(binary, tmpdir, ref_seq, markers):
    """Vary number of genomes, fixed genome size."""
    print("\n=== classify: scaling by number of genomes ===")
    ref_path = str(tmpdir / "reference.fasta")
    markers_path = str(tmpdir / "markers.tsv")
    rows = []
    for n in CLASSIFY_N:
        rng = random.Random(SEED + n)
        fp = str(tmpdir / f"cls_n{n}.fasta")
        make_multifasta(fp, ref_seq, markers, n, rng)
        mb = os.path.getsize(fp) / (1024**2)
        times, rsses = [], []
        for rep in range(REPEATS):
            op = str(tmpdir / f"cls_n{n}_r{rep}")
            t, rss = run_timed([binary, "classify", "-m", markers_path, "-r", ref_path, "-i", fp, "-o", op])
            if t is not None:
                times.append(t)
                if rss is not None:
                    rsses.append(rss)
        if times:
            rows.append({"n_genomes": n, "input_mb": round(mb, 1),
                          "median_s": np.median(times), "std_s": np.std(times),
                          "median_rss_mb": np.median(rsses) if rsses else None})
            print(f"  {n:3d} genomes: {np.median(times):.3f}s, RSS {np.median(rsses):.0f} MB" if rsses else
                  f"  {n:3d} genomes: {np.median(times):.3f}s")
    return pd.DataFrame(rows)


def bench_classify_by_genome_size(binary, tmpdir):
    """Vary genome size, fixed 5 genomes."""
    print("\n=== classify: scaling by genome size ===")
    N_GENOMES = 5
    rows = []
    for gs, gs_label in zip(GENOME_SIZES, GENOME_SIZE_LABELS):
        rng = random.Random(SEED + gs)
        ref_path = str(tmpdir / f"ref_gs{gs}.fasta")
        mk_path = str(tmpdir / f"markers_gs{gs}.tsv")
        ref_seq = make_reference(ref_path, gs, rng)
        markers = make_markers(mk_path, ref_seq, rng, n_markers=min(3700, gs // 1200))
        fp = str(tmpdir / f"cls_gs{gs}.fasta")
        make_multifasta(fp, ref_seq, markers, N_GENOMES, rng)
        mb = os.path.getsize(fp) / (1024**2)
        times, rsses = [], []
        for rep in range(REPEATS):
            op = str(tmpdir / f"cls_gs{gs}_r{rep}")
            t, rss = run_timed([binary, "classify", "-m", mk_path, "-r", ref_path, "-i", fp, "-o", op])
            if t is not None:
                times.append(t)
                if rss is not None:
                    rsses.append(rss)
        if times:
            rows.append({"genome_mb": round(gs / 1e6, 1), "genome_label": gs_label,
                          "input_mb": round(mb, 1),
                          "median_s": np.median(times), "std_s": np.std(times),
                          "median_rss_mb": np.median(rsses) if rsses else None})
            rss_str = f", RSS {np.median(rsses):.0f} MB" if rsses else ""
            print(f"  {gs_label:>4s} Mb genome: {np.median(times):.3f}s{rss_str}")
    return pd.DataFrame(rows)


def bench_splitfastq_by_count(binary, tmpdir, ref_seq, markers):
    """Vary read count, fixed genome size."""
    print("\n=== split-fastq: scaling by number of reads ===")
    ref_path = str(tmpdir / "reference.fasta")
    markers_path = str(tmpdir / "markers.tsv")
    rows = []
    for nr in SPLITFASTQ_N:
        rng = random.Random(SEED + nr)
        r1 = str(tmpdir / f"fq_{nr}_R1.fastq.gz")
        r2 = str(tmpdir / f"fq_{nr}_R2.fastq.gz")
        print(f"  Generating {nr:,} reads...", end=" ", flush=True)
        make_fastq(r1, r2, ref_seq, nr, rng)
        mb = (os.path.getsize(r1) + os.path.getsize(r2)) / (1024**2)
        print(f"({mb:.1f} MB)")
        times, rsses = [], []
        for rep in range(REPEATS):
            op = str(tmpdir / f"fq_{nr}_r{rep}")
            t, rss = run_timed([binary, "split-fastq", "-m", markers_path, "-r", ref_path,
                                "-i", r1, "-i", r2, "--paired", "-o", op])
            if t is not None:
                times.append(t)
                if rss is not None:
                    rsses.append(rss)
        for f in [r1, r2]:
            if os.path.exists(f):
                os.remove(f)
        if times:
            rows.append({"n_reads": nr, "input_mb": round(mb, 1),
                          "median_s": np.median(times), "std_s": np.std(times),
                          "median_rss_mb": np.median(rsses) if rsses else None})
            rss_str = f", RSS {np.median(rsses):.0f} MB" if rsses else ""
            print(f"  {nr:>9,} reads: {np.median(times):.3f}s{rss_str}")
    return pd.DataFrame(rows)


def bench_splitfastq_by_genome_size(binary, tmpdir):
    """Vary genome size, fixed 250k reads."""
    print("\n=== split-fastq: scaling by genome size ===")
    N_READS = 250_000
    rows = []
    for gs, gs_label in zip(GENOME_SIZES, GENOME_SIZE_LABELS):
        rng = random.Random(SEED + gs)
        ref_path = str(tmpdir / f"ref_gs{gs}.fasta")
        mk_path = str(tmpdir / f"markers_gs{gs}.tsv")
        # Reuse if already generated by classify
        if not os.path.exists(ref_path):
            ref_seq = make_reference(ref_path, gs, rng)
            make_markers(mk_path, ref_seq, rng, n_markers=min(3700, gs // 1200))
        else:
            with open(ref_path) as f:
                ref_seq = "".join(line.strip() for line in f if not line.startswith(">"))
        r1 = str(tmpdir / f"fq_gs{gs}_R1.fastq.gz")
        r2 = str(tmpdir / f"fq_gs{gs}_R2.fastq.gz")
        print(f"  Generating reads for {gs_label} Mb genome...", end=" ", flush=True)
        make_fastq(r1, r2, ref_seq, N_READS, rng)
        mb = (os.path.getsize(r1) + os.path.getsize(r2)) / (1024**2)
        print(f"({mb:.1f} MB)")
        times, rsses = [], []
        for rep in range(REPEATS):
            op = str(tmpdir / f"fq_gs{gs}_r{rep}")
            t, rss = run_timed([binary, "split-fastq", "-m", mk_path, "-r", ref_path,
                                "-i", r1, "-i", r2, "--paired", "-o", op])
            if t is not None:
                times.append(t)
                if rss is not None:
                    rsses.append(rss)
        for f in [r1, r2]:
            if os.path.exists(f):
                os.remove(f)
        if times:
            rows.append({"genome_mb": round(gs / 1e6, 1), "genome_label": gs_label,
                          "input_mb": round(mb, 1),
                          "median_s": np.median(times), "std_s": np.std(times),
                          "median_rss_mb": np.median(rsses) if rsses else None})
            rss_str = f", RSS {np.median(rsses):.0f} MB" if rsses else ""
            print(f"  {gs_label:>4s} Mb genome: {np.median(times):.3f}s{rss_str}")
    return pd.DataFrame(rows)


def bench_match_by_count(binary, tmpdir, ref_seq):
    """Vary number of references, fixed genome size."""
    print("\n=== match: scaling by number of references ===")
    rng = random.Random(SEED)
    q1 = str(tmpdir / "query_R1.fastq.gz")
    q2 = str(tmpdir / "query_R2.fastq.gz")
    if not os.path.exists(q1):
        make_fastq(q1, q2, ref_seq, 100_000, rng)
    rows = []
    for nr in MATCH_N:
        rng = random.Random(SEED + nr)
        db = str(tmpdir / f"refdb_{nr}.fasta")
        make_ref_db(db, ref_seq, nr, rng)
        db_mb = os.path.getsize(db) / (1024**2)
        times, rsses = [], []
        for rep in range(REPEATS):
            for idx in Path(tmpdir).glob(f"refdb_{nr}*.ptidx"):
                idx.unlink()
            op = str(tmpdir / f"match_{nr}_r{rep}.tsv")
            t, rss = run_timed([binary, "match", "-r", db, "-i", q1, "-i", q2, "-o", op])
            if t is not None:
                times.append(t)
                if rss is not None:
                    rsses.append(rss)
        if times:
            rows.append({"n_refs": nr, "db_mb": round(db_mb, 1),
                          "median_s": np.median(times), "std_s": np.std(times),
                          "median_rss_mb": np.median(rsses) if rsses else None})
            rss_str = f", RSS {np.median(rsses):.0f} MB" if rsses else ""
            print(f"  {nr:3d} refs: {np.median(times):.3f}s{rss_str}")
    return pd.DataFrame(rows)


def bench_match_by_genome_size(binary, tmpdir):
    """Vary genome size, fixed 5 references."""
    print("\n=== match: scaling by genome/reference size ===")
    N_REFS = 5
    N_QUERY_READS = 100_000
    rows = []
    for gs, gs_label in zip(GENOME_SIZES, GENOME_SIZE_LABELS):
        rng = random.Random(SEED + gs)
        ref_path = str(tmpdir / f"ref_gs{gs}.fasta")
        if not os.path.exists(ref_path):
            ref_seq = make_reference(ref_path, gs, rng)
        else:
            with open(ref_path) as f:
                ref_seq = "".join(line.strip() for line in f if not line.startswith(">"))
        db = str(tmpdir / f"mdb_gs{gs}.fasta")
        make_ref_db(db, ref_seq, N_REFS, rng)
        q1 = str(tmpdir / f"mq_gs{gs}_R1.fastq.gz")
        q2 = str(tmpdir / f"mq_gs{gs}_R2.fastq.gz")
        make_fastq(q1, q2, ref_seq, N_QUERY_READS, rng)
        db_mb = os.path.getsize(db) / (1024**2)
        times, rsses = [], []
        for rep in range(REPEATS):
            for idx in Path(tmpdir).glob(f"mdb_gs{gs}*.ptidx"):
                idx.unlink()
            op = str(tmpdir / f"match_gs{gs}_r{rep}.tsv")
            t, rss = run_timed([binary, "match", "-r", db, "-i", q1, "-i", q2, "-o", op])
            if t is not None:
                times.append(t)
                if rss is not None:
                    rsses.append(rss)
        for f in [q1, q2]:
            if os.path.exists(f):
                os.remove(f)
        if times:
            rows.append({"genome_mb": round(gs / 1e6, 1), "genome_label": gs_label,
                          "db_mb": round(db_mb, 1),
                          "median_s": np.median(times), "std_s": np.std(times),
                          "median_rss_mb": np.median(rsses) if rsses else None})
            rss_str = f", RSS {np.median(rsses):.0f} MB" if rsses else ""
            print(f"  {gs_label:>4s} Mb genome: {np.median(times):.3f}s{rss_str}")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
def _save(fig, outdir, name, source_df=None):
    for ext in (".pdf", ".png"):
        fig.savefig(outdir / f"{name}{ext}")
    plt.close(fig)
    print(f"  -> {name}.pdf/.png")
    if source_df is not None:
        src = outdir / f"{name}_data.csv"
        source_df.to_csv(src, index=False)


def _fmt_count(n):
    """Format large counts for x-axis: 50000 → '50k', 1000000 → '1M'."""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.0f}M"
    elif n >= 1_000:
        return f"{n / 1_000:.0f}k"
    return str(int(n))


def _setup_axis(ax, color, title):
    """Common axis setup: grid, title. Y-limits set after data is plotted."""
    ax.set_title(title, fontweight="bold", fontsize=10)
    ax.grid(axis="y", linewidth=0.3, alpha=0.4, zorder=0)


def _finalize_ylim(ax, data_max):
    """Set y-axis from 0 to data_max * 1.15."""
    if data_max > 0:
        ax.set_ylim(0, data_max * 1.15)


def _annotate_last(ax, x, y, color, fmt="{:.2f}s"):
    """Annotate the last data point with its value, avoiding collisions."""
    if len(x) > 0 and len(y) > 0:
        ax.annotate(fmt.format(y.iloc[-1]),
                    xy=(x.iloc[-1], y.iloc[-1]),
                    xytext=(-5, 8), textcoords="offset points",
                    fontsize=6.5, color=color, fontweight="bold",
                    ha="right")


# ---------------------------------------------------------------------------
# Figure 1: Time scaling — 2 rows (by count, by genome size) x 3 cols
# ---------------------------------------------------------------------------
def plot_time_scaling(df_cls_n, df_fq_n, df_match_n,
                      df_cls_gs, df_fq_gs, df_match_gs, outdir):
    """2x3 grid: top row = scaling by count, bottom = scaling by genome size."""
    print("\nGenerating time scaling plots...")
    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE_COL * 1.1, DOUBLE_COL * 0.55))

    # --- Top row: scaling by input count ---
    panels_n = [
        (axes[0, 0], df_cls_n, "n_genomes", "Number of genomes", "classify"),
        (axes[0, 1], df_fq_n, "n_reads", "Number of reads", "split-fastq"),
        (axes[0, 2], df_match_n, "n_refs", "Number of references", "match"),
    ]

    for ax, df, xcol, xlabel, name in panels_n:
        color = MODULE_COLORS[name]
        _setup_axis(ax, color, name)
        if df.empty:
            continue

        ax.fill_between(df[xcol],
                        df["median_s"] - df["std_s"],
                        df["median_s"] + df["std_s"],
                        alpha=0.15, color=color, zorder=2)
        ax.plot(df[xcol], df["median_s"], "o-", color=color, markersize=4,
                linewidth=1.5, zorder=3)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("Time (s)", fontsize=8)
        _finalize_ylim(ax, df["median_s"].max())
        _annotate_last(ax, df[xcol], df["median_s"], color)

        # Explicit x-ticks for all panels to avoid auto-ticks with decimals
        ax.set_xticks(df[xcol].values)
        if name == "split-fastq":
            ax.set_xticklabels([_fmt_count(n) for n in df[xcol]], fontsize=7, rotation=40, ha="right")
        else:
            ax.set_xticklabels([_fmt_count(n) for n in df[xcol]], fontsize=7)

    # --- Bottom row: scaling by genome size ---
    panels_gs = [
        (axes[1, 0], df_cls_gs, "classify (5 genomes)"),
        (axes[1, 1], df_fq_gs, "split-fastq (250k reads)"),
        (axes[1, 2], df_match_gs, "match (5 refs)"),
    ]

    for ax, df, title in panels_gs:
        mod = title.split(" ")[0]
        color = MODULE_COLORS[mod]
        _setup_axis(ax, color, title)
        if df.empty:
            continue

        ax.fill_between(df["genome_mb"],
                        df["median_s"] - df["std_s"],
                        df["median_s"] + df["std_s"],
                        alpha=0.15, color=color, zorder=2)
        ax.plot(df["genome_mb"], df["median_s"], "o-", color=color, markersize=4,
                linewidth=1.5, zorder=3)
        ax.set_xlabel("Genome size (Mb)", fontsize=8)
        ax.set_ylabel("Time (s)", fontsize=8)
        ax.set_xticks(df["genome_mb"].values)
        ax.set_xticklabels(df["genome_label"].values if "genome_label" in df.columns
                           else [str(v) for v in df["genome_mb"]], fontsize=7)
        _finalize_ylim(ax, df["median_s"].max())
        _annotate_last(ax, df["genome_mb"], df["median_s"], color)

    fig.suptitle("Execution Time Scaling", fontweight="bold", fontsize=12, y=1.02)
    fig.tight_layout()
    _save(fig, outdir, "scaling_time")


# ---------------------------------------------------------------------------
# Figure 2: RSS scaling — 2 rows x 3 cols (same layout as time)
# ---------------------------------------------------------------------------
def plot_rss_scaling(df_cls_n, df_fq_n, df_match_n,
                     df_cls_gs, df_fq_gs, df_match_gs, outdir):
    """2x3 grid: top row = RSS by count, bottom = RSS by genome size."""
    print("Generating RSS scaling plots...")

    # Check if any RSS data exists
    all_dfs = [df_cls_n, df_fq_n, df_match_n, df_cls_gs, df_fq_gs, df_match_gs]
    has_rss = any("median_rss_mb" in df.columns and df["median_rss_mb"].notna().any()
                  for df in all_dfs if not df.empty)
    if not has_rss:
        print("  [skip] no RSS data available")
        return

    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE_COL * 1.1, DOUBLE_COL * 0.55))

    rss_color = "#666666"

    # --- Top row: RSS by input count ---
    panels_n = [
        (axes[0, 0], df_cls_n, "n_genomes", "Number of genomes", "classify"),
        (axes[0, 1], df_fq_n, "n_reads", "Number of reads", "split-fastq"),
        (axes[0, 2], df_match_n, "n_refs", "Number of references", "match"),
    ]

    for ax, df, xcol, xlabel, name in panels_n:
        color = MODULE_COLORS[name]
        _setup_axis(ax, color, name)
        if df.empty or "median_rss_mb" not in df.columns:
            continue

        rss = df["median_rss_mb"]
        if rss.isna().all():
            continue

        ax.fill_between(df[xcol], 0, rss, alpha=0.2, color=color, zorder=2)
        ax.plot(df[xcol], rss, "s-", color=color, markersize=4,
                linewidth=1.5, zorder=3)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("Peak RSS (MB)", fontsize=8)
        _finalize_ylim(ax, rss.max())

        # Annotate last point
        _annotate_last(ax, df[xcol], rss, color, fmt="{:.0f} MB")

        # Explicit x-ticks
        ax.set_xticks(df[xcol].values)
        if name == "split-fastq":
            ax.set_xticklabels([_fmt_count(n) for n in df[xcol]], fontsize=7, rotation=40, ha="right")
        else:
            ax.set_xticklabels([_fmt_count(n) for n in df[xcol]], fontsize=7)

    # --- Bottom row: RSS by genome size ---
    panels_gs = [
        (axes[1, 0], df_cls_gs, "classify (5 genomes)"),
        (axes[1, 1], df_fq_gs, "split-fastq (250k reads)"),
        (axes[1, 2], df_match_gs, "match (5 refs)"),
    ]

    for ax, df, title in panels_gs:
        mod = title.split(" ")[0]
        color = MODULE_COLORS[mod]
        _setup_axis(ax, color, title)
        if df.empty or "median_rss_mb" not in df.columns:
            continue

        rss = df["median_rss_mb"]
        if rss.isna().all():
            continue

        ax.fill_between(df["genome_mb"], 0, rss, alpha=0.2, color=color, zorder=2)
        ax.plot(df["genome_mb"], rss, "s-", color=color, markersize=4,
                linewidth=1.5, zorder=3)
        ax.set_xlabel("Genome size (Mb)", fontsize=8)
        ax.set_ylabel("Peak RSS (MB)", fontsize=8)
        ax.set_xticks(df["genome_mb"].values)
        ax.set_xticklabels(df["genome_label"].values if "genome_label" in df.columns
                           else [str(v) for v in df["genome_mb"]], fontsize=7)
        _finalize_ylim(ax, rss.max())
        _annotate_last(ax, df["genome_mb"], rss, color, fmt="{:.0f} MB")

    fig.suptitle("Peak Memory (RSS) Scaling", fontweight="bold", fontsize=12, y=1.02)
    fig.tight_layout()
    _save(fig, outdir, "scaling_rss")


# ---------------------------------------------------------------------------
# Figure 3: Combined summary — one panel per module, time + RSS overlaid
# ---------------------------------------------------------------------------
def plot_combined_summary(df_cls_n, df_fq_n, df_match_n, outdir):
    """3-panel summary: each module's time and RSS on the same chart (by count)."""
    print("Generating combined summary...")
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL * 1.1, SINGLE_COL * 1.0))

    panels = [
        (axes[0], df_cls_n, "n_genomes", "Number of genomes", "classify"),
        (axes[1], df_fq_n, "n_reads", "Number of reads", "split-fastq"),
        (axes[2], df_match_n, "n_refs", "Number of references", "match"),
    ]

    for ax, df, xcol, xlabel, name in panels:
        color = MODULE_COLORS[name]
        ax.set_title(name, fontweight="bold", fontsize=10)
        if df.empty:
            continue

        # Time on left axis
        ax.fill_between(df[xcol],
                        df["median_s"] - df["std_s"],
                        df["median_s"] + df["std_s"],
                        alpha=0.15, color=color, zorder=2)
        ln1 = ax.plot(df[xcol], df["median_s"], "o-", color=color, markersize=4,
                      linewidth=1.5, zorder=3, label="Time")
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("Time (s)", color=color, fontsize=8)
        ax.tick_params(axis="y", labelcolor=color)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", linewidth=0.3, alpha=0.3, zorder=0)

        # Explicit x-ticks
        ax.set_xticks(df[xcol].values)
        if name == "split-fastq":
            ax.set_xticklabels([_fmt_count(n) for n in df[xcol]], fontsize=7, rotation=40, ha="right")
        else:
            ax.set_xticklabels([_fmt_count(n) for n in df[xcol]], fontsize=7)

        # RSS on right axis — only if RSS varies meaningfully (>10% range)
        if "median_rss_mb" in df.columns and df["median_rss_mb"].notna().any():
            rss = df["median_rss_mb"]
            rss_range = rss.max() - rss.min()
            if rss_range > rss.mean() * 0.1:  # Only show if >10% variation
                ax2 = ax.twinx()
                ln2 = ax2.plot(df[xcol], rss, "s--", color="#888", markersize=3,
                               linewidth=1.0, alpha=0.7, label="Peak RSS")
                ax2.set_ylabel("Peak RSS (MB)", color="#888", fontsize=8)
                ax2.tick_params(axis="y", labelcolor="#888")
                ax2.spines["right"].set_visible(True)
                ax2.spines["right"].set_linewidth(0.5)
                ax2.set_ylim(bottom=0)

                # Combined legend
                lines = ln1 + ln2
                labels = [l.get_label() for l in lines]
                ax.legend(lines, labels, loc="upper left", fontsize=6.5, framealpha=0.9)
            else:
                # Show flat RSS as text annotation instead
                ax.text(0.98, 0.05, f"RSS: {rss.mean():.0f} MB (const.)",
                        transform=ax.transAxes, ha="right", fontsize=6.5,
                        color="#888", fontstyle="italic")

    fig.suptitle("Scaling by Input Count", fontweight="bold", fontsize=12, y=1.03)
    fig.tight_layout()
    _save(fig, outdir, "scaling_combined_summary")


# ---------------------------------------------------------------------------
# Re-plot from saved CSV data
# ---------------------------------------------------------------------------
def _replot_from_csv(outdir):
    """Re-generate plots from previously saved CSV data files."""
    count_csv = outdir / "scaling_by_count_data.csv"
    gs_csv = outdir / "scaling_by_genome_size_data.csv"

    if not count_csv.exists() and not gs_csv.exists():
        print("ERROR: No CSV data files found. Run benchmarks first (without --plot-only).")
        sys.exit(1)

    # Parse count-scaling data
    df_cls_n = df_fq_n = df_match_n = pd.DataFrame()
    if count_csv.exists():
        df = pd.read_csv(count_csv)
        if "module" in df.columns:
            for mod, sub in df.groupby("module"):
                sub = sub.copy()
                if "median_rss_mb" not in sub.columns:
                    sub["median_rss_mb"] = None
                if mod == "classify" and "input_param" in sub.columns:
                    sub["n_genomes"] = sub["input_param"].astype(int)
                    df_cls_n = sub
                elif mod == "split-fastq" and "input_param" in sub.columns:
                    sub["n_reads"] = sub["input_param"].astype(int)
                    df_fq_n = sub
                elif mod == "match" and "input_param" in sub.columns:
                    sub["n_refs"] = sub["input_param"].astype(int)
                    df_match_n = sub

    # Parse genome-size scaling data
    df_cls_gs = df_fq_gs = df_match_gs = pd.DataFrame()
    if gs_csv.exists():
        df = pd.read_csv(gs_csv)
        if "module" in df.columns:
            for mod, sub in df.groupby("module"):
                sub = sub.copy()
                if "median_rss_mb" not in sub.columns:
                    sub["median_rss_mb"] = None
                if mod == "classify":
                    df_cls_gs = sub
                elif mod == "split-fastq":
                    df_fq_gs = sub
                elif mod == "match":
                    df_match_gs = sub

    print()
    plot_time_scaling(df_cls_n, df_fq_n, df_match_n,
                      df_cls_gs, df_fq_gs, df_match_gs, outdir)
    plot_rss_scaling(df_cls_n, df_fq_n, df_match_n,
                     df_cls_gs, df_fq_gs, df_match_gs, outdir)
    plot_combined_summary(df_cls_n, df_fq_n, df_match_n, outdir)
    print("Done!")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Scaling benchmark for pathotypr")
    parser.add_argument("--binary", default="pathotypr")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--keep-tmp", action="store_true")
    parser.add_argument("--plot-only", action="store_true",
                        help="Re-generate plots from saved CSV data without running benchmarks")
    args = parser.parse_args()

    # Find binary
    binary = shutil.which(args.binary) or args.binary
    try:
        v = subprocess.run([binary, "--version"], capture_output=True, text=True)
        print(f"pathotypr: {v.stdout.strip()}")
    except FileNotFoundError:
        script_dir = Path(__file__).resolve().parent
        release_bin = script_dir.parent.parent / "target" / "release" / "pathotypr"
        if release_bin.exists():
            binary = str(release_bin)
            v = subprocess.run([binary, "--version"], capture_output=True, text=True)
            print(f"pathotypr: {v.stdout.strip()} (from target/release/)")
        else:
            print("ERROR: pathotypr not found. Build: cargo build --release -p pathotypr-core")
            sys.exit(1)

    if TIME_CMD:
        print(f"Time cmd:  {TIME_CMD} ({TIME_FLAVOR})")
    else:
        print("WARNING: no /usr/bin/time or gtime found — RSS will not be measured")

    outdir = Path(args.outdir) if args.outdir else \
        Path(__file__).resolve().parent.parent / "target" / "criterion" / "scaling"
    outdir.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        print(f"Plot-only mode. Reading CSVs from: {outdir}")
        _replot_from_csv(outdir)
        return

    tmpdir = Path(tempfile.mkdtemp(prefix="pathotypr_bench_"))
    print(f"Temp:      {tmpdir}")
    print(f"Output:    {outdir}")

    # --- Base data (4.4 Mb genome) ---
    rng = random.Random(SEED)
    print("\nGenerating base reference (4.4 Mb) and markers...")
    ref_seq = make_reference(str(tmpdir / "reference.fasta"), DEFAULT_GENOME_SIZE, rng)
    markers = make_markers(str(tmpdir / "markers.tsv"), ref_seq, rng)

    # --- Part 1: Scaling by count ---
    df_cls_n = bench_classify_by_count(binary, tmpdir, ref_seq, markers)
    df_fq_n = bench_splitfastq_by_count(binary, tmpdir, ref_seq, markers)
    df_match_n = bench_match_by_count(binary, tmpdir, ref_seq)

    # --- Part 2: Scaling by genome size ---
    df_cls_gs = bench_classify_by_genome_size(binary, tmpdir)
    df_fq_gs = bench_splitfastq_by_genome_size(binary, tmpdir)
    df_match_gs = bench_match_by_genome_size(binary, tmpdir)

    # --- Save raw CSVs for --plot-only reuse ---
    parts_n = []
    for df, xcol, name in [(df_cls_n, "n_genomes", "classify"),
                            (df_fq_n, "n_reads", "split-fastq"),
                            (df_match_n, "n_refs", "match")]:
        if not df.empty:
            d = df.copy()
            d["module"] = name
            d["input_param"] = d[xcol]
            parts_n.append(d)
    if parts_n:
        pd.concat(parts_n, ignore_index=True).to_csv(
            outdir / "scaling_by_count_data.csv", index=False)

    parts_gs = []
    for df, name in [(df_cls_gs, "classify"), (df_fq_gs, "split-fastq"), (df_match_gs, "match")]:
        if not df.empty:
            d = df.copy()
            d["module"] = name
            parts_gs.append(d)
    if parts_gs:
        pd.concat(parts_gs, ignore_index=True).to_csv(
            outdir / "scaling_by_genome_size_data.csv", index=False)

    # --- Plots ---
    print()
    plot_time_scaling(df_cls_n, df_fq_n, df_match_n,
                      df_cls_gs, df_fq_gs, df_match_gs, outdir)
    plot_rss_scaling(df_cls_n, df_fq_n, df_match_n,
                     df_cls_gs, df_fq_gs, df_match_gs, outdir)
    plot_combined_summary(df_cls_n, df_fq_n, df_match_n, outdir)

    if not args.keep_tmp:
        shutil.rmtree(tmpdir)
        print(f"\nCleaned up {tmpdir}")

    print(f"\nResults in: {outdir}/")
    print("Done!")


if __name__ == "__main__":
    main()
