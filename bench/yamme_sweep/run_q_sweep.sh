#!/usr/bin/env bash
# Run yamme with current defaults (no -q filter) on every fixture so that
# score_q_sweep.py can decide a reasonable default for the -q flag.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
OUT_DIR="$HERE/results/q_sweep"
mkdir -p "$OUT_DIR"

# Make sure yamtk is built at current defaults.
(cd "$REPO" && make -s yamtk)

DATASETS=(
  D1 D2 D3 D4 D5 D6 D7 D8
  D9   D9_v2  D9_v3  D9_v4  D9_v5
  D10  D10_v2 D10_v3 D10_v4 D10_v5
  D11  D11_v2 D11_v3 D11_v4 D11_v5
  D12 D13 D14
  D15 D16
)

for d in "${DATASETS[@]}"; do
  "$REPO/yamtk" me -i "$FIX_DIR/${d}_pos.fa" \
    -k 6 -K 15 -N 8 -s 42 \
    -o "$OUT_DIR/${d}.tsv" -M '' > /dev/null 2>&1
  echo "  [done] $d"
done

echo "TSVs written to: $OUT_DIR"
