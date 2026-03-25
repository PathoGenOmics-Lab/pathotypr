#!/usr/bin/env bash
# =============================================================================
# Run pathotypr benchmarks and generate publication figures
#
# Usage:
#   bash run_benchmarks.sh              # Run all with real data (TEST_PATHOTYPR)
#   bash run_benchmarks.sh --synthetic  # Generate synthetic data first
#   bash run_benchmarks.sh --plot-only  # Only regenerate plots from existing results
#   bash run_benchmarks.sh --filter "split_kmer"  # Run specific group only
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CORE_DIR="$(dirname "$SCRIPT_DIR")"
WORKSPACE_DIR="$(dirname "$CORE_DIR")"

USE_SYNTHETIC=false
PLOT_ONLY=false
FILTER=""

for arg in "$@"; do
    case "$arg" in
        --synthetic) USE_SYNTHETIC=true ;;
        --plot-only) PLOT_ONLY=true ;;
        --filter) shift; FILTER="$1" ;;
        --filter=*) FILTER="${arg#--filter=}" ;;
        *) ;;
    esac
done

echo "========================================"
echo "  pathotypr benchmark suite"
echo "========================================"
echo ""

# --- Step 1: Generate synthetic data if requested ---
if $USE_SYNTHETIC; then
    SYNTH_DIR="${WORKSPACE_DIR}/TEST_SYNTHETIC"
    echo "Step 1: Generating synthetic test data..."
    python3 "${SCRIPT_DIR}/generate_synthetic_data.py" --outdir "$SYNTH_DIR"
    echo ""

    # Point benchmarks to synthetic data
    export PATHOTYPR_TEST_DIR="$SYNTH_DIR"
    echo "Using synthetic data from: $SYNTH_DIR"
else
    REAL_DIR="${WORKSPACE_DIR}/TEST_PATHOTYPR"
    if [ ! -d "$REAL_DIR" ]; then
        echo "WARNING: TEST_PATHOTYPR not found at $REAL_DIR"
        echo "Use --synthetic to generate test data."
        exit 1
    fi
    echo "Using real data from: $REAL_DIR"
fi

echo ""

# --- Step 2: Run benchmarks ---
if ! $PLOT_ONLY; then
    echo "Step 2: Running Criterion benchmarks..."
    cd "$WORKSPACE_DIR"

    if [ -n "$FILTER" ]; then
        echo "  Filter: $FILTER"
        cargo bench --bench benchmarks -p pathotypr-core -- "$FILTER"
    else
        cargo bench --bench benchmarks -p pathotypr-core
    fi

    echo ""
fi

# --- Step 3: Generate plots ---
echo "Step 3: Generating publication figures..."
cd "$CORE_DIR"
python3 "${SCRIPT_DIR}/plot_benchmarks.py"

echo ""
echo "========================================"
echo "  Results"
echo "========================================"
echo ""
echo "Criterion HTML reports:"
echo "  ${WORKSPACE_DIR}/target/criterion/report/index.html"
echo ""
echo "Publication figures (Nature style):"
ls -1 "${WORKSPACE_DIR}/target/criterion/figures/"*.pdf 2>/dev/null | while read f; do
    echo "  $f"
done
echo ""
echo "Source data (CSV):"
ls -1 "${WORKSPACE_DIR}/target/criterion/figures/"*.csv 2>/dev/null | while read f; do
    echo "  $f"
done
