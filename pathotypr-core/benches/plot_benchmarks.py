#!/usr/bin/env python3
"""
Parse Criterion.rs benchmark results and generate clean publication figures.

Reads from: target/criterion/<group>/<bench>/new/estimates.json
Outputs to: target/criterion/figures/

Requires: matplotlib, numpy, pandas
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Style: clean, readable, Nature-inspired but with larger fonts for clarity
# ---------------------------------------------------------------------------
MM_TO_IN = 1 / 25.4
SINGLE_COL = 89 * MM_TO_IN   # 3.50 in
DOUBLE_COL = 183 * MM_TO_IN  # 7.20 in

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "DejaVu Sans", "Arial"],
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7,
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
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Palette — muted, colorblind-friendly
PALETTE = {
    "io":                   "#4E79A7",
    "marker_generation":    "#59A14F",
    "find_markers":         "#E15759",
    "reverse_complement":   "#76B7B2",
    "analyze_genome":       "#F28E2B",
    "split_kmer":           "#B07AA1",
    "vectorizer":           "#4E4E4E",
    "train":                "#FF9DA7",
    "predict":              "#9C755F",
    "match":                "#EDC948",
}

# Friendly module labels
MODULE_LABELS = {
    "io":                   "I/O",
    "marker_generation":    "Marker Gen",
    "find_markers":         "Find Markers",
    "reverse_complement":   "Rev. Comp.",
    "analyze_genome":       "Genome Analysis",
    "split_kmer":           "Split k-mer",
    "vectorizer":           "Vectorizer",
    "train":                "Train",
    "predict":              "Predict",
    "match":                "Match",
}


def auto_unit(ns):
    """Convert nanoseconds to the best human unit, return (value, unit_str)."""
    if ns < 1_000:
        return ns, "ns"
    elif ns < 1_000_000:
        return ns / 1_000, "us"
    elif ns < 1_000_000_000:
        return ns / 1_000_000, "ms"
    else:
        return ns / 1_000_000_000, "s"


def format_time(ns):
    """Format nanoseconds to a short, readable string."""
    val, unit = auto_unit(ns)
    if val >= 100:
        return f"{val:.0f} {unit}"
    elif val >= 10:
        return f"{val:.1f} {unit}"
    else:
        return f"{val:.2f} {unit}"


def shorten_label(name):
    """Shorten benchmark names for cleaner labels."""
    # Remove __ artifacts from Criterion paths
    name = name.replace("__", "::")
    # Truncate if too long
    if len(name) > 45:
        name = name[:42] + "..."
    return name


def collect_results(criterion_dir):
    """Collect all benchmark results, keeping only the latest run per benchmark."""
    results = []
    criterion_path = Path(criterion_dir)

    if not criterion_path.exists():
        print(f"ERROR: Criterion directory not found: {criterion_path}")
        sys.exit(1)

    for group_dir in sorted(criterion_path.iterdir()):
        if not group_dir.is_dir() or group_dir.name in ("report", "figures"):
            continue

        group_name = group_dir.name
        for bench_dir in sorted(group_dir.iterdir()):
            if not bench_dir.is_dir():
                continue

            estimates_file = bench_dir / "new" / "estimates.json"
            if not estimates_file.exists():
                continue

            with open(estimates_file) as f:
                data = json.load(f)

            mean_ns = data["mean"]["point_estimate"]
            ci_lower = data["mean"]["confidence_interval"]["lower_bound"]
            ci_upper = data["mean"]["confidence_interval"]["upper_bound"]
            std_ns = data["std_dev"]["point_estimate"]

            results.append({
                "group": group_name,
                "benchmark": bench_dir.name,
                "mean_ns": mean_ns,
                "std_ns": std_ns,
                "ci_lower_ns": ci_lower,
                "ci_upper_ns": ci_upper,
            })

    df = pd.DataFrame(results)

    # Deduplicate: if there are old + new versions of same benchmark,
    # keep only the one with the most recent estimates.json mtime
    if not df.empty:
        # Group by (group, base_name) to detect duplicates like
        # "get_positions (markers)" vs "get_positions (3708 markers)"
        # For now, just keep all — they may be intentionally different params
        pass

    return df


def _save(fig, out_dir, name, source_df=None):
    """Save figure as PDF + PNG, optionally export source CSV."""
    for ext in (".pdf", ".png"):
        fig.savefig(out_dir / f"{name}{ext}")
    plt.close(fig)
    print(f"  -> {name}.pdf/.png")
    if source_df is not None:
        src = out_dir / f"{name}_data.csv"
        source_df.to_csv(src, index=False)


def _format_axis_time(ax, axis="x"):
    """Format an axis with human-readable time ticks."""
    def _fmt(x, _pos):
        return format_time(x)
    if axis == "x":
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(_fmt))
    else:
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(_fmt))


# =====================================================================
# Figure 1: Top Bottlenecks — the 12 slowest benchmarks
# =====================================================================
def plot_top_bottlenecks(df, out_dir):
    """Horizontal bar chart of top N slowest benchmarks (linear scale)."""
    print("Fig 1: Top bottlenecks...")

    top = df.nlargest(12, "mean_ns").sort_values("mean_ns", ascending=True)
    n = len(top)
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, max(2.0, 0.28 * n + 0.6)))

    y = np.arange(n)
    colors = [PALETTE.get(g, "#888") for g in top["group"]]

    bars = ax.barh(y, top["mean_ns"], color=colors, height=0.65, zorder=3)

    # CI error bars
    err_lower = top["mean_ns"].values - top["ci_lower_ns"].values
    err_upper = top["ci_upper_ns"].values - top["mean_ns"].values
    ax.errorbar(top["mean_ns"], y, xerr=[err_lower, err_upper],
                fmt="none", color="#333", capsize=2, linewidth=0.6, capthick=0.5, zorder=4)

    # Labels with module prefix
    labels = [f"{shorten_label(row['benchmark'])}" for _, row in top.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)

    # Time labels on bars
    x_max = top["mean_ns"].max()
    for i, (_, row) in enumerate(top.iterrows()):
        val = row["mean_ns"]
        # Place label inside or outside bar depending on size
        if val > x_max * 0.6:
            ax.text(val - x_max * 0.02, i, format_time(val),
                    va="center", ha="right", fontsize=7, fontweight="bold", color="white")
        else:
            ax.text(val + x_max * 0.015, i, format_time(val),
                    va="center", ha="left", fontsize=7, color="#333")

    ax.set_xlabel("Mean execution time")
    _format_axis_time(ax, "x")
    ax.set_xlim(0, x_max * 1.15)
    ax.grid(axis="x", linewidth=0.3, alpha=0.4, zorder=0)

    # Legend for module colors
    from matplotlib.patches import Patch
    seen = []
    handles = []
    for _, row in top.iterrows():
        g = row["group"]
        if g not in seen:
            seen.append(g)
            handles.append(Patch(facecolor=PALETTE.get(g, "#888"),
                                 label=MODULE_LABELS.get(g, g)))
    ax.legend(handles=handles, loc="lower right", framealpha=0.9, edgecolor="#ccc",
              fontsize=6.5, title="Module", title_fontsize=7)

    ax.set_title("Top Bottlenecks", fontweight="bold", fontsize=11)
    fig.tight_layout()

    export = top.copy()
    export["mean_formatted"] = export["mean_ns"].apply(format_time)
    _save(fig, out_dir, "fig1_top_bottlenecks", source_df=export)


# =====================================================================
# Figure 2: Pipeline Breakdowns (split-fastq, train, classify)
# =====================================================================
def plot_pipeline_breakdowns(df, out_dir):
    """Multi-panel pipeline breakdown with proportional (linear) bars."""
    print("Fig 2: Pipeline breakdowns...")

    pipelines = {
        "split-fastq": [
            ("io", "get_ref", "Load reference"),
            ("io", "get_positions", "Parse markers"),
            ("split_kmer", "build_markers", "Build markers"),
            ("split_kmer", "build_marker_index", "Build index"),
            ("split_kmer", "scan_fastq_with_index", "Scan FASTQ"),
        ],
        "train": [
            ("train", "read_fasta", "Read FASTA"),
            ("train", "parallel_kmer_discovery", "K-mer discovery"),
            ("train", "prepare_data", "Prepare data"),
        ],
        "classify (FASTA)": [
            ("io", "get_ref", "Load reference"),
            ("io", "get_positions", "Parse markers"),
            ("marker_generation", "generate_markerkmer", "Generate k-mers"),
            ("find_markers", "find_markers", "Find markers"),
        ],
    }

    valid_pipelines = {}
    for pipe_name, steps in pipelines.items():
        rows = []
        for group, bench_prefix, label in steps:
            match = df[(df["group"] == group) & (df["benchmark"].str.startswith(bench_prefix))]
            if not match.empty:
                # Take the first match (usually the most representative)
                row = match.iloc[0]
                rows.append({"step": label, "mean_ns": row["mean_ns"],
                             "group": group})
        if len(rows) >= 2:
            valid_pipelines[pipe_name] = pd.DataFrame(rows)

    if not valid_pipelines:
        print("  [skip] no valid pipelines found")
        return

    n_pipes = len(valid_pipelines)
    fig, axes = plt.subplots(1, n_pipes, figsize=(DOUBLE_COL, 2.4))
    if n_pipes == 1:
        axes = [axes]

    for ax, (pipe_name, pipe_df) in zip(axes, valid_pipelines.items()):
        pipe_df = pipe_df.sort_values("mean_ns", ascending=True)
        total = pipe_df["mean_ns"].sum()
        n = len(pipe_df)
        y = np.arange(n)

        colors = [PALETTE.get(g, "#888") for g in pipe_df["group"]]
        ax.barh(y, pipe_df["mean_ns"], color=colors, height=0.55, zorder=3)

        x_max = pipe_df["mean_ns"].max()
        for i, (_, row) in enumerate(pipe_df.iterrows()):
            pct = row["mean_ns"] / total * 100
            label = f"{format_time(row['mean_ns'])}"
            if pct >= 1.0:
                label += f" ({pct:.0f}%)"
            ax.text(row["mean_ns"] + x_max * 0.03, i, label,
                    va="center", fontsize=6.5, color="#333")

        ax.set_yticks(y)
        ax.set_yticklabels(pipe_df["step"], fontsize=7)
        ax.set_xlim(0, x_max * 1.45)
        ax.set_xlabel("")
        _format_axis_time(ax, "x")
        ax.grid(axis="x", linewidth=0.3, alpha=0.3, zorder=0)
        ax.set_title(pipe_name, fontweight="bold", fontsize=9)

    fig.suptitle("Pipeline Breakdown", fontweight="bold", fontsize=11, y=1.02)
    fig.tight_layout()
    _save(fig, out_dir, "fig2_pipeline_breakdowns")


# =====================================================================
# Figure 3: Module Facets — one small panel per module
# =====================================================================
def plot_module_facets(df, out_dir):
    """Grid of small panels, one per module, each with its own linear scale."""
    print("Fig 3: Module facets...")

    groups = sorted(df["group"].unique(), key=lambda g: -df[df["group"] == g]["mean_ns"].max())
    n_groups = len(groups)
    n_cols = 3
    n_rows = (n_groups + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(DOUBLE_COL, n_rows * 1.4 + 0.3))
    axes = axes.flatten() if n_groups > 1 else [axes]

    for idx, group in enumerate(groups):
        ax = axes[idx]
        gdf = df[df["group"] == group].sort_values("mean_ns", ascending=True)
        n = len(gdf)
        y = np.arange(n)
        color = PALETTE.get(group, "#888")

        ax.barh(y, gdf["mean_ns"], color=color, height=0.6, alpha=0.85, zorder=3)

        # CI error bars
        err_lo = gdf["mean_ns"].values - gdf["ci_lower_ns"].values
        err_hi = gdf["ci_upper_ns"].values - gdf["mean_ns"].values
        ax.errorbar(gdf["mean_ns"], y, xerr=[err_lo, err_hi],
                    fmt="none", color="#333", capsize=1.5, linewidth=0.4, capthick=0.4, zorder=4)

        short_names = [shorten_label(name) for name in gdf["benchmark"]]
        ax.set_yticks(y)
        ax.set_yticklabels(short_names, fontsize=5.5)

        x_max = gdf["mean_ns"].max()
        for i, (_, row) in enumerate(gdf.iterrows()):
            ax.text(row["mean_ns"] + x_max * 0.04, i, format_time(row["mean_ns"]),
                    va="center", fontsize=5.5, color="#555")

        ax.set_xlim(0, x_max * 1.4)
        _format_axis_time(ax, "x")
        ax.grid(axis="x", linewidth=0.2, alpha=0.3, zorder=0)
        ax.set_title(MODULE_LABELS.get(group, group), fontweight="bold", fontsize=8,
                     color=color)
        ax.tick_params(axis="x", labelsize=5.5)

    # Hide unused axes
    for idx in range(n_groups, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Benchmarks by Module", fontweight="bold", fontsize=11, y=1.01)
    fig.tight_layout()
    _save(fig, out_dir, "fig3_module_facets", source_df=df)


# =====================================================================
# Figure 4: Compact Overview — all benchmarks, log scale, well-organized
# =====================================================================
def plot_full_overview(df, out_dir):
    """All benchmarks in a single chart, grouped by module with separators."""
    print("Fig 4: Full overview...")

    # Sort: modules by max time (descending), benchmarks within by time
    module_order = (df.groupby("group")["mean_ns"].max()
                     .sort_values(ascending=False).index.tolist())

    rows = []
    for group in module_order:
        gdf = df[df["group"] == group].sort_values("mean_ns", ascending=False)
        for _, row in gdf.iterrows():
            rows.append(row)
    ordered = pd.DataFrame(rows).reset_index(drop=True)

    n = len(ordered)
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, max(3.5, 0.22 * n + 0.8)))

    y_pos = np.arange(n)
    colors = [PALETTE.get(g, "#888") for g in ordered["group"]]

    ax.barh(y_pos, ordered["mean_ns"], color=colors, height=0.7, zorder=3)

    # CI error bars
    err_lo = ordered["mean_ns"].values - ordered["ci_lower_ns"].values
    err_hi = ordered["ci_upper_ns"].values - ordered["mean_ns"].values
    ax.errorbar(ordered["mean_ns"], y_pos, xerr=[err_lo, err_hi],
                fmt="none", color="#333", capsize=1.5, linewidth=0.4, capthick=0.4, zorder=4)

    # Labels
    labels = [shorten_label(row["benchmark"]) for _, row in ordered.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=5.5)

    # Time annotations
    for i, (_, row) in enumerate(ordered.iterrows()):
        ax.text(row["mean_ns"] * 1.25, i, format_time(row["mean_ns"]),
                va="center", fontsize=5, color="#444")

    ax.set_xlabel("Mean execution time (log scale)")
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.grid(axis="x", linewidth=0.2, alpha=0.3, zorder=0)

    # Draw module group separators and labels
    prev_group = None
    group_start = 0
    for i, (_, row) in enumerate(ordered.iterrows()):
        if row["group"] != prev_group:
            if prev_group is not None:
                sep_y = i - 0.5
                ax.axhline(sep_y, color="#ccc", linewidth=0.5, zorder=1)
            prev_group = row["group"]
            group_start = i

    # Module color legend
    from matplotlib.patches import Patch
    seen = []
    handles = []
    for g in module_order:
        if g not in seen:
            seen.append(g)
            handles.append(Patch(facecolor=PALETTE.get(g, "#888"),
                                 label=MODULE_LABELS.get(g, g)))
    ax.legend(handles=handles, loc="lower right", framealpha=0.95, edgecolor="#ccc",
              fontsize=5.5, title="Module", title_fontsize=6.5, ncol=2)

    ax.set_title("All Benchmarks", fontweight="bold", fontsize=11)
    fig.tight_layout()

    export = ordered.copy()
    export["mean_formatted"] = export["mean_ns"].apply(format_time)
    _save(fig, out_dir, "fig4_full_overview", source_df=export)


# =====================================================================
# Figure 5: Module comparison — heaviest single benchmark per module
# =====================================================================
def plot_module_comparison(df, out_dir):
    """Compare modules by their heaviest benchmark (the true bottleneck)."""
    print("Fig 5: Module bottleneck comparison...")

    heaviest = df.loc[df.groupby("group")["mean_ns"].idxmax()]
    heaviest = heaviest.sort_values("mean_ns", ascending=True)

    n = len(heaviest)
    fig, ax = plt.subplots(figsize=(SINGLE_COL * 1.3, max(1.8, 0.3 * n + 0.4)))

    y = np.arange(n)
    colors = [PALETTE.get(g, "#888") for g in heaviest["group"]]

    ax.barh(y, heaviest["mean_ns"], color=colors, height=0.55, zorder=3)

    # CI error bars
    err_lo = heaviest["mean_ns"].values - heaviest["ci_lower_ns"].values
    err_hi = heaviest["ci_upper_ns"].values - heaviest["mean_ns"].values
    ax.errorbar(heaviest["mean_ns"], y, xerr=[err_lo, err_hi],
                fmt="none", color="#333", capsize=2, linewidth=0.5, capthick=0.5, zorder=4)

    # Module labels on y-axis, benchmark name as annotation
    labels = [MODULE_LABELS.get(row["group"], row["group"]) for _, row in heaviest.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontweight="bold")

    x_max = heaviest["mean_ns"].max()
    for i, (_, row) in enumerate(heaviest.iterrows()):
        val = row["mean_ns"]
        bench = shorten_label(row["benchmark"])
        if val > x_max * 0.5:
            ax.text(val - x_max * 0.02, i, format_time(val),
                    va="center", ha="right", fontsize=7, fontweight="bold", color="white")
            ax.text(val - x_max * 0.02, i + 0.32, bench,
                    va="center", ha="right", fontsize=5, color="white", alpha=0.85)
        else:
            ax.text(val + x_max * 0.015, i, format_time(val),
                    va="center", fontsize=7, color="#333")
            ax.text(val + x_max * 0.015, i + 0.32, bench,
                    va="center", fontsize=5, color="#888")

    ax.set_xlim(0, x_max * 1.12)
    ax.set_xlabel("Heaviest benchmark per module")
    _format_axis_time(ax, "x")
    ax.grid(axis="x", linewidth=0.3, alpha=0.3, zorder=0)
    ax.set_title("Module Bottlenecks", fontweight="bold", fontsize=10)
    fig.tight_layout()

    export = heaviest.copy()
    export["mean_formatted"] = export["mean_ns"].apply(format_time)
    _save(fig, out_dir, "fig5_module_bottlenecks", source_df=export)


def find_criterion_dir():
    """Locate the criterion output directory."""
    script_dir = Path(__file__).resolve().parent
    # script_dir = .../pathotypr-core/benches/
    # workspace  = .../pathotypr-final/   (two levels up)
    candidates = [
        script_dir.parent.parent / "target" / "criterion",      # from benches/ → workspace/target/
        script_dir.parent / "target" / "criterion",             # pathotypr-core/target/ (unlikely)
        Path.cwd() / "target" / "criterion",                    # from workspace root (cargo bench cwd)
        Path.cwd() / "pathotypr-final" / "target" / "criterion", # from project root
    ]
    for c in candidates:
        if c.exists() and any(c.iterdir()):
            return c
    return None


def main():
    criterion_dir = find_criterion_dir()
    if criterion_dir is None:
        print("ERROR: Could not find target/criterion/")
        print("Run `cargo bench --bench benchmarks` first.")
        sys.exit(1)

    out_dir = criterion_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {criterion_dir}")
    print(f"Output:  {out_dir}\n")

    df = collect_results(criterion_dir)
    if df.empty:
        print("No benchmark results found!")
        sys.exit(1)

    print(f"Found {len(df)} benchmarks in {df['group'].nunique()} modules\n")

    # Generate all figures
    plot_top_bottlenecks(df, out_dir)
    plot_pipeline_breakdowns(df, out_dir)
    plot_module_facets(df, out_dir)
    plot_full_overview(df, out_dir)
    plot_module_comparison(df, out_dir)

    # Export raw data
    raw = df.copy()
    raw["mean_formatted"] = raw["mean_ns"].apply(format_time)
    raw.to_csv(out_dir / "all_benchmarks_raw.csv", index=False)
    print(f"\nRaw data: {out_dir / 'all_benchmarks_raw.csv'}")
    print("Done!")


if __name__ == "__main__":
    main()
