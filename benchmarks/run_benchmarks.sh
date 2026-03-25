#!/usr/bin/env bash
# Run pathotypr benchmarks at various data sizes.
# Outputs JSON results to benchmarks/results.json
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
BIN="$REPO_DIR/target/release/pathotypr"
RESULTS="$SCRIPT_DIR/results.json"
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Build release binary
echo "🔨 Building release binary..."
cd "$REPO_DIR"
cargo build --release -p pathotypr-core --bin pathotypr 2>/dev/null

# Generate synthetic FASTA at various sizes
generate_fasta() {
    local path=$1 n_classes=$2 n_per_class=$3
    python3 -c "
import random
random.seed(42)
motifs = [
    'ATCGATCGATCGATCGATCGATCGATCG',
    'GCTAGCTAGCTAGCTAGCTAGCTAGCTA',
    'TTAACCGGTTAACCGGTTAACCGGTTAA',
    'CCGGAATTCCGGAATTCCGGAATTCCGG',
    'AACCTTGGAACCTTGGAACCTTGGAACC',
    'GGTTCCAAGGTTCCAAGGTTCCAAGGTT',
    'TCGATCGATCGATCGATCGATCGATCGA',
    'AGCTAGCTAGCTAGCTAGCTAGCTAGCT',
]
with open('$path', 'w') as f:
    for ci in range($n_classes):
        motif = motifs[ci % len(motifs)]
        for i in range($n_per_class):
            f.write(f'>Class{ci} sample_{ci}_{i}\n')
            seq = ''
            for j in range(40):
                seq += motif
                base = 'ACGT'[(i + j) % 4]
                if j % 3 == 0:
                    seq += base
            f.write(seq[:1000] + '\n')
"
}

echo "[]" > "$RESULTS"

# Benchmark configurations: (n_classes, n_per_class, total_genomes)
SIZES=(
    "3 10 30"
    "3 30 90"
    "3 100 300"
    "5 100 500"
    "5 200 1000"
    "5 400 2000"
    "8 500 4000"
)

echo "📊 Running benchmarks..."
echo "[" > "$RESULTS"
first=true

for config in "${SIZES[@]}"; do
    read -r nc npc total <<< "$config"
    
    FASTA="$TMPDIR/bench_${total}.fasta"
    MODEL="$TMPDIR/model_${total}.pathotypr.zst"
    PRED_OUT="$TMPDIR/pred_${total}.tsv"
    
    echo "  → $total genomes ($nc classes × $npc each)..."
    generate_fasta "$FASTA" "$nc" "$npc"
    
    # Measure file size
    fasta_size=$(wc -c < "$FASTA")
    
    # TRAIN benchmark (3 runs, take median)
    train_times=()
    for run in 1 2 3; do
        rm -f "$MODEL" "$MODEL".*
        t_start=$(python3 -c "import time; print(time.time())")
        "$BIN" train -i "$FASTA" -o "$MODEL" -k 11 --max-depth 10 --min-samples-leaf 2 -t 4 2>/dev/null
        t_end=$(python3 -c "import time; print(time.time())")
        elapsed=$(python3 -c "print(round($t_end - $t_start, 4))")
        train_times+=("$elapsed")
    done
    train_median=$(python3 -c "ts=sorted([${train_times[0]},${train_times[1]},${train_times[2]}]); print(ts[1])")
    
    # Model size
    model_size=$(wc -c < "$MODEL")
    
    # PREDICT benchmark (3 runs, take median)
    predict_times=()
    for run in 1 2 3; do
        rm -f "$PRED_OUT"
        t_start=$(python3 -c "import time; print(time.time())")
        "$BIN" predict -i "$FASTA" -m "$MODEL" -o "$PRED_OUT" -t 4 2>/dev/null
        t_end=$(python3 -c "import time; print(time.time())")
        elapsed=$(python3 -c "print(round($t_end - $t_start, 4))")
        predict_times+=("$elapsed")
    done
    predict_median=$(python3 -c "ts=sorted([${predict_times[0]},${predict_times[1]},${predict_times[2]}]); print(ts[1])")
    
    # Memory (peak RSS via /usr/bin/time if available)
    if command -v /usr/bin/time &>/dev/null; then
        train_mem=$( { /usr/bin/time -l "$BIN" train -i "$FASTA" -o "$MODEL" -k 11 --max-depth 10 --min-samples-leaf 2 -t 4 2>&1 >/dev/null; } 2>&1 | grep 'maximum resident' | awk '{print $1}' || echo "0")
        predict_mem=$( { /usr/bin/time -l "$BIN" predict -i "$FASTA" -m "$MODEL" -o "$PRED_OUT" -t 4 2>&1 >/dev/null; } 2>&1 | grep 'maximum resident' | awk '{print $1}' || echo "0")
    else
        train_mem="0"
        predict_mem="0"
    fi
    
    if [ "$first" = true ]; then first=false; else echo "," >> "$RESULTS"; fi
    cat >> "$RESULTS" <<EOF
  {
    "n_genomes": $total,
    "n_classes": $nc,
    "n_per_class": $npc,
    "fasta_size_bytes": $fasta_size,
    "train_time_s": $train_median,
    "predict_time_s": $predict_median,
    "model_size_bytes": $model_size,
    "train_peak_rss_bytes": $train_mem,
    "predict_peak_rss_bytes": $predict_mem,
    "threads": 4,
    "kmer_size": 11,
    "seq_length": 1000
  }
EOF
    echo "    train=${train_median}s predict=${predict_median}s model=$(echo "scale=1; $model_size/1024" | bc)KB"
done

echo "]" >> "$RESULTS"
echo "✅ Results written to $RESULTS"
