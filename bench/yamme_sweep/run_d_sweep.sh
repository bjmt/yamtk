#!/usr/bin/env bash
# Sweep -D (cross-width dedup overlap threshold) at current defaults
# (-t 1e-3, -q 1 to study dedup in isolation; q-filter is applied at scoring).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
OUT_DIR="$HERE/results/d_sweep"
TIMING="$HERE/results/timings_d_sweep.tsv"

mkdir -p "$OUT_DIR"
(cd "$REPO" && make -s yamtk)
: > "$TIMING"
echo -e "dataset\td\twall_s" >> "$TIMING"

DS=(0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

DATASETS=(
  D1 D2 D3 D4 D5 D6 D7 D8
  D9   D9_v2  D9_v3  D9_v4  D9_v5
  D10  D10_v2 D10_v3 D10_v4 D10_v5
  D11  D11_v2 D11_v3 D11_v4 D11_v5
  D12 D13 D14
  D15 D16
)

for d in "${DS[@]}"; do
  echo "[D=${d}]"
  for ds in "${DATASETS[@]}"; do
    tsv="$OUT_DIR/${ds}_D${d}.tsv"
    TIMEFORMAT='%R'
    wall=$({ time "$REPO/yamtk" me -i "$FIX_DIR/${ds}_pos.fa" \
      -k 6 -K 15 -N 8 -s 42 -D "$d" -q 1 \
      -o "$tsv" -M '' > /dev/null 2>&1; } 2>&1)
    echo -e "${ds}\t${d}\t${wall}" >> "$TIMING"
  done
done
echo "Done.  TSVs: $OUT_DIR    timings: $TIMING"
