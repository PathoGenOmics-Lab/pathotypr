#!/usr/bin/env python3
"""Generate benchmark visualizations for pathotypr docs."""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
RESULTS_FILE = SCRIPT_DIR / "results.json"
OUTPUT_DIR = SCRIPT_DIR / "figures"
OUTPUT_DIR.mkdir(exist_ok=True)

# Theme support (light / dark)
import sys
DARK = '--dark' in sys.argv

LIGHT_STYLE = {
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': '#333333',
    'axes.labelcolor': '#333333',
    'text.color': '#333333',
    'xtick.color': '#333333',
    'ytick.color': '#333333',
}

DARK_STYLE = {
    'figure.facecolor': '#1a1a2e',
    'axes.facecolor': '#1a1a2e',
    'axes.edgecolor': '#cccccc',
    'axes.labelcolor': '#e0e0e0',
    'text.color': '#e0e0e0',
    'xtick.color': '#cccccc',
    'ytick.color': '#cccccc',
}

# Style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 12,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.2,
})
plt.rcParams.update(DARK_STYLE if DARK else LIGHT_STYLE)

SUFFIX = '-dark' if DARK else ''

# Colors
TRAIN_COLOR = '#6366f1'   # indigo
PREDICT_COLOR = '#f43f5e'  # rose
MODEL_COLOR = '#10b981'    # emerald
MEM_TRAIN_COLOR = '#8b5cf6'  # violet
MEM_PREDICT_COLOR = '#f97316'  # orange

with open(RESULTS_FILE) as f:
    data = json.load(f)

genomes = [d['n_genomes'] for d in data]
train_times = [d['train_time_s'] for d in data]
predict_times = [d['predict_time_s'] for d in data]
model_sizes_kb = [d['model_size_bytes'] / 1024 for d in data]
train_mem_mb = [d['train_peak_rss_bytes'] / (1024 * 1024) if d['train_peak_rss_bytes'] > 0 else None for d in data]
predict_mem_mb = [d['predict_peak_rss_bytes'] / (1024 * 1024) if d['predict_peak_rss_bytes'] > 0 else None for d in data]

# --- Figure 1: Training & Prediction Time ---
fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(genomes, train_times, 'o-', color=TRAIN_COLOR, linewidth=2.5, markersize=8, label='Train (100 trees)', zorder=5)
ax.plot(genomes, predict_times, 's-', color=PREDICT_COLOR, linewidth=2.5, markersize=8, label='Predict', zorder=5)

# Fill area
ax.fill_between(genomes, train_times, alpha=0.12, color=TRAIN_COLOR)
ax.fill_between(genomes, predict_times, alpha=0.12, color=PREDICT_COLOR)

# Annotate last points
ax.annotate(f'{train_times[-1]:.2f}s', (genomes[-1], train_times[-1]),
            textcoords="offset points", xytext=(10, 5), fontsize=11, color=TRAIN_COLOR, fontweight='bold')
ax.annotate(f'{predict_times[-1]:.3f}s', (genomes[-1], predict_times[-1]),
            textcoords="offset points", xytext=(10, -15), fontsize=11, color=PREDICT_COLOR, fontweight='bold')

ax.set_xlabel('Number of genomes', fontsize=13, fontweight='bold')
ax.set_ylabel('Time (seconds)', fontsize=13, fontweight='bold')
ax.set_title('Training & Prediction Speed', fontsize=15, fontweight='bold', pad=15)
ax.legend(fontsize=12, frameon=False, )
ax.set_xlim(0, max(genomes) * 1.15)
ax.set_ylim(0, max(train_times) * 1.2)
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))

fig.savefig(OUTPUT_DIR / f'speed{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'speed{SUFFIX}.png'}")

# --- Figure 2: Throughput (genomes/second) ---
fig, ax = plt.subplots(figsize=(9, 5))

train_throughput = [g / t for g, t in zip(genomes, train_times)]
predict_throughput = [g / t for g, t in zip(genomes, predict_times)]

bars_x = np.arange(len(genomes))
width = 0.35

bars1 = ax.bar(bars_x - width/2, train_throughput, width, color=TRAIN_COLOR, alpha=0.85, label='Train', zorder=5)
bars2 = ax.bar(bars_x + width/2, predict_throughput, width, color=PREDICT_COLOR, alpha=0.85, label='Predict', zorder=5)

# Value labels on bars
for bar, val in zip(bars1, train_throughput):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(predict_throughput)*0.02,
            f'{val:.0f}', ha='center', va='bottom', fontsize=9, color=TRAIN_COLOR, fontweight='bold')
for bar, val in zip(bars2, predict_throughput):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(predict_throughput)*0.02,
            f'{val:,.0f}', ha='center', va='bottom', fontsize=9, color=PREDICT_COLOR, fontweight='bold')

ax.set_xlabel('Dataset size', fontsize=13, fontweight='bold')
ax.set_ylabel('Genomes / second', fontsize=13, fontweight='bold')
ax.set_title('Throughput', fontsize=15, fontweight='bold', pad=15)
ax.set_xticks(bars_x)
ax.set_xticklabels([f'{g:,}' for g in genomes], rotation=30, ha='right')
ax.legend(fontsize=12, frameon=False, )
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))

fig.savefig(OUTPUT_DIR / f'throughput{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'throughput{SUFFIX}.png'}")

# --- Figure 3: Model Size ---
fig, ax = plt.subplots(figsize=(9, 4))

ax.bar(range(len(genomes)), model_sizes_kb, color=MODEL_COLOR, alpha=0.85, zorder=5)
for i, (g, s) in enumerate(zip(genomes, model_sizes_kb)):
    ax.text(i, s + max(model_sizes_kb)*0.03, f'{s:.1f} KB', ha='center', va='bottom',
            fontsize=10, color=MODEL_COLOR, fontweight='bold')

ax.set_xlabel('Training genomes', fontsize=13, fontweight='bold')
ax.set_ylabel('Model size (KB)', fontsize=13, fontweight='bold')
ax.set_title('Compressed Model Size (zstd)', fontsize=15, fontweight='bold', pad=15)
ax.set_xticks(range(len(genomes)))
ax.set_xticklabels([f'{g:,}' for g in genomes], rotation=30, ha='right')
ax.set_ylim(0, max(model_sizes_kb) * 1.25)
ax.grid(axis='y', alpha=0.3, linestyle='--')

fig.savefig(OUTPUT_DIR / f'model_size{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'model_size{SUFFIX}.png'}")

# --- Figure 4: Peak Memory (if available) ---
if any(m is not None for m in train_mem_mb):
    fig, ax = plt.subplots(figsize=(9, 5))
    
    valid_genomes = [g for g, m in zip(genomes, train_mem_mb) if m is not None]
    valid_train = [m for m in train_mem_mb if m is not None]
    valid_predict = [m for m, t in zip(predict_mem_mb, train_mem_mb) if t is not None]
    
    ax.plot(valid_genomes, valid_train, 'o-', color=MEM_TRAIN_COLOR, linewidth=2.5, markersize=8, label='Train', zorder=5)
    ax.plot(valid_genomes, valid_predict, 's-', color=MEM_PREDICT_COLOR, linewidth=2.5, markersize=8, label='Predict', zorder=5)
    ax.fill_between(valid_genomes, valid_train, alpha=0.12, color=MEM_TRAIN_COLOR)
    ax.fill_between(valid_genomes, valid_predict, alpha=0.12, color=MEM_PREDICT_COLOR)
    
    ax.annotate(f'{valid_train[-1]:.0f} MB', (valid_genomes[-1], valid_train[-1]),
                textcoords="offset points", xytext=(10, 5), fontsize=11, color=MEM_TRAIN_COLOR, fontweight='bold')
    ax.annotate(f'{valid_predict[-1]:.0f} MB', (valid_genomes[-1], valid_predict[-1]),
                textcoords="offset points", xytext=(10, -15), fontsize=11, color=MEM_PREDICT_COLOR, fontweight='bold')
    
    ax.set_xlabel('Number of genomes', fontsize=13, fontweight='bold')
    ax.set_ylabel('Peak RSS (MB)', fontsize=13, fontweight='bold')
    ax.set_title('Peak Memory Usage', fontsize=15, fontweight='bold', pad=15)
    ax.legend(fontsize=12, frameon=False, )
    ax.set_xlim(0, max(valid_genomes) * 1.15)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    fig.savefig(OUTPUT_DIR / f'memory{SUFFIX}.png')
    plt.close()
    print(f"  ✅ {OUTPUT_DIR / f'memory{SUFFIX}.png'}")

# --- Figure 5: Combined dashboard ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Speed
ax = axes[0]
ax.plot(genomes, train_times, 'o-', color=TRAIN_COLOR, linewidth=2, markersize=6, label='Train')
ax.plot(genomes, predict_times, 's-', color=PREDICT_COLOR, linewidth=2, markersize=6, label='Predict')
ax.fill_between(genomes, train_times, alpha=0.1, color=TRAIN_COLOR)
ax.fill_between(genomes, predict_times, alpha=0.1, color=PREDICT_COLOR)
ax.set_xlabel('Genomes')
ax.set_ylabel('Time (s)')
ax.set_title('Speed', fontweight='bold')
ax.legend(fontsize=10, frameon=False, )
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Throughput
ax = axes[1]
ax.bar(range(len(genomes)), predict_throughput, color=PREDICT_COLOR, alpha=0.8, label='Predict')
ax.set_xlabel('Dataset')
ax.set_ylabel('Genomes/s')
ax.set_title('Predict Throughput', fontweight='bold')
ax.set_xticks(range(len(genomes)))
ax.set_xticklabels([f'{g:,}' for g in genomes], rotation=45, ha='right', fontsize=9)
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))

# Model size
ax = axes[2]
ax.bar(range(len(genomes)), model_sizes_kb, color=MODEL_COLOR, alpha=0.8)
ax.set_xlabel('Training genomes')
ax.set_ylabel('Size (KB)')
ax.set_title('Model Size', fontweight='bold')
ax.set_xticks(range(len(genomes)))
ax.set_xticklabels([f'{g:,}' for g in genomes], rotation=45, ha='right', fontsize=9)
ax.grid(axis='y', alpha=0.3, linestyle='--')

fig.suptitle('pathotypr Performance Benchmarks', fontsize=16, fontweight='bold', y=1.02)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / f'dashboard{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'dashboard{SUFFIX}.png'}")

# --- Figure 6: All Modules — Real MTB (2 panels: time + memory) ---
modules = ['train\n(10)', 'train\n(50)', 'predict\n(5)', 'classify\n(5)', 'split-fastq\n(500K)', 'split-fastq\n(full)', 'match\n(20 refs)']
times_real = [0.61, 55.03, 0.26, 0.10, 1.45, 10.46, 78.11]
mem_real = [302, 1381, 198, 92, 26, 26, 4606]
time_labels = ['0.61s', '55.0s', '0.26s', '0.10s', '1.4s', '10.5s', '78.1s']
mem_labels = ['302 MB', '1.4 GB', '198 MB', '92 MB', '26 MB', '26 MB', '4.6 GB']
bar_colors = ['#6366f1', '#8b5cf6', '#10b981', '#f43f5e', '#f59e0b', '#3b82f6', '#f97316']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('pathotypr — Real M. tuberculosis Benchmarks (Mac mini M4, 4 threads)', fontsize=15, fontweight='bold', y=1.02)

for ax, vals, labels, ylabel, title in [
    (ax1, times_real, time_labels, 'Time (seconds)', 'Execution Time'),
    (ax2, mem_real, mem_labels, 'Peak RAM (MB)', 'Peak Memory Usage')]:
    bars = ax.bar(range(len(modules)), vals, color=bar_colors, alpha=0.85, zorder=5)
    ax.set_yscale('log')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xticks(range(len(modules)))
    ax.set_xticklabels(modules, fontsize=9)
    ax.grid(axis='y', alpha=0.3, linestyle='--', zorder=0)
    for bar, lbl in zip(bars, labels):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.15, lbl,
                ha='center', va='bottom', fontsize=9, fontweight='bold')

fig.tight_layout()
fig.savefig(OUTPUT_DIR / f'all_modules{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'all_modules{SUFFIX}.png'}")

# --- Figure 7: Train scaling real MTB ---
train_genomes_real = [10, 15, 20, 25, 30, 35, 40, 45, 50]
train_times_real = [0.6, 3.8, 8.2, 13.5, 20.5, 27.8, 35.2, 44.0, 55.0]
train_ram_real = [0.30, 0.38, 0.50, 0.62, 0.75, 0.88, 1.00, 1.18, 1.38]

fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()

ln1 = ax1.plot(train_genomes_real, train_times_real, 'o-', color='#6366f1', linewidth=2.5, markersize=8, label='Time (s)', zorder=5)
ln2 = ax2.plot(train_genomes_real, train_ram_real, 's-', color='#f43f5e', linewidth=2.5, markersize=8, label='Peak RAM (GB)', zorder=5)

ax1.set_xlabel('Number of genomes', fontsize=13, fontweight='bold')
ax1.set_ylabel('Time (seconds)', fontsize=13, fontweight='bold', color='#6366f1')
ax2.set_ylabel('Peak RAM (GB)', fontsize=13, fontweight='bold', color='#f43f5e')
ax1.set_title('train: Scaling with Dataset Size (real MTB, k=21)', fontsize=14, fontweight='bold', pad=15)
ax1.grid(axis='y', alpha=0.3, linestyle='--')

ax1.annotate('0.6s', (10, 0.6), textcoords="offset points", xytext=(-15, 10), fontsize=11, color='#6366f1', fontweight='bold')
ax1.annotate('55.0s', (50, 55.0), textcoords="offset points", xytext=(5, 5), fontsize=11, color='#6366f1', fontweight='bold')

lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, fontsize=11, loc='upper left')

fig.tight_layout()
fig.savefig(OUTPUT_DIR / f'train_scaling_real{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'train_scaling_real{SUFFIX}.png'}")

# --- Figure 8: Split-FASTQ scaling ---
sfq_reads = [1.0, 4.25, 8.5]
sfq_times = [1.45, 5.1, 10.46]
sfq_ram = [26, 26, 26]

fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()

ln1 = ax1.plot(sfq_reads, sfq_times, 'o-', color='#3b82f6', linewidth=2.5, markersize=8, label='Time (s)', zorder=5)
ln2 = ax2.plot(sfq_reads, sfq_ram, 's--', color='#f43f5e', linewidth=3, markersize=8, label='Peak RAM (MB)', zorder=5)

ax1.set_xlabel('Reads (millions)', fontsize=13, fontweight='bold')
ax1.set_ylabel('Time (seconds)', fontsize=13, fontweight='bold', color='#3b82f6')
ax2.set_ylabel('Peak RAM (MB)', fontsize=13, fontweight='bold', color='#f43f5e')
ax2.set_ylim(0, 58)
ax1.set_title('split-fastq: Linear Time, Constant Memory', fontsize=14, fontweight='bold', pad=15)
ax1.grid(axis='y', alpha=0.3, linestyle='--')

lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, fontsize=11, loc='upper left')

fig.tight_layout()
fig.savefig(OUTPUT_DIR / f'split_fastq_scaling{SUFFIX}.png')
plt.close()
print(f"  ✅ {OUTPUT_DIR / f'split_fastq_scaling{SUFFIX}.png'}")

# --- Figure 9: pathotypr vs fastlin comparison ---
COMP_FILE = SCRIPT_DIR / "results_comparison.json"
if COMP_FILE.exists():
    with open(COMP_FILE) as f:
        comp = json.load(f)

    samples = [c['sample'] for c in comp]
    sample_labels = [f"{c['sample']}\n({c['fastq_size_mb']:.0f} MB)" for c in comp]
    pt_times = [c['pathotypr_time_s'] for c in comp]
    fl_times = [c['fastlin_time_s'] for c in comp]
    pt_mem = [c['pathotypr_rss_mb'] for c in comp]
    fl_mem = [c['fastlin_rss_mb'] for c in comp]

    PT_COLOR = '#6366f1'
    FL_COLOR = '#f97316'

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('pathotypr vs fastlin — Real TB FASTQ Samples', fontsize=15, fontweight='bold', y=1.02)

    x = np.arange(len(samples))
    w = 0.35

    for ax, pt_vals, fl_vals, ylabel, title, fmt in [
        (ax1, pt_times, fl_times, 'Time (seconds)', 'Speed Comparison', '{:.1f}s'),
        (ax2, pt_mem, fl_mem, 'Peak RSS (MB)', 'Memory Usage', '{:.0f} MB')]:
        bars1 = ax.bar(x - w/2, pt_vals, w, color=PT_COLOR, alpha=0.85, label='pathotypr', zorder=5)
        bars2 = ax.bar(x + w/2, fl_vals, w, color=FL_COLOR, alpha=0.85, label='fastlin', zorder=5)
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(sample_labels, fontsize=10)
        ax.legend(fontsize=11, frameon=False, )
        ax.grid(axis='y', alpha=0.3, linestyle='--', zorder=0)
        for bar, val in zip(bars1, pt_vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(pt_vals + fl_vals)*0.03,
                    fmt.format(val), ha='center', va='bottom', fontsize=10, color=PT_COLOR, fontweight='bold')
        for bar, val in zip(bars2, fl_vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(pt_vals + fl_vals)*0.03,
                    fmt.format(val), ha='center', va='bottom', fontsize=10, color=FL_COLOR, fontweight='bold')

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f'comparison{SUFFIX}.png')
    plt.close()
    print(f"  ✅ {OUTPUT_DIR / f'comparison{SUFFIX}.png'}")

print("\n📊 All figures generated!")
