#!/usr/bin/env bash
# Sweep -t (per-motif stop p-value) at current defaults (incl. -q 1e-3) over
# all fixtures, to pick a sensible default for -t given the -q filter.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
OUT_DIR="$HERE/results/t_sweep"
TIMING="$HERE/results/timings_t_sweep.tsv"

mkdir -p "$OUT_DIR"
(cd "$REPO" && make -s yamtk)
: > "$TIMING"
echo -e "dataset\tt\twall_s" >> "$TIMING"

# Disable the q-filter so we can study -t in isolation; q-filter will be
# applied at scoring time so we can also see how -t affects post-q-filter
# recall.
DISABLE_Q="-q 1"

TS=(1.0 0.5 0.05 0.01 0.001 0.0001)

DATASETS=(
  D1 D2 D3 D4 D5 D6 D7 D8
  D9   D9_v2  D9_v3  D9_v4  D9_v5
  D10  D10_v2 D10_v3 D10_v4 D10_v5
  D11  D11_v2 D11_v3 D11_v4 D11_v5
  D12 D13 D14
  D15 D16
)

for t in "${TS[@]}"; do
  echo "[t=${t}]"
  for d in "${DATASETS[@]}"; do
    tsv="$OUT_DIR/${d}_t${t}.tsv"
    TIMEFORMAT='%R'
    wall=$({ time "$REPO/yamtk" me -i "$FIX_DIR/${d}_pos.fa" \
      -k 6 -K 15 -N 8 -s 42 -t "$t" $DISABLE_Q \
      -o "$tsv" -M '' > /dev/null 2>&1; } 2>&1)
    echo -e "${d}\t${t}\t${wall}" >> "$TIMING"
  done
done
echo "Done.  TSVs: $OUT_DIR    timings: $TIMING"
