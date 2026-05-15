#!/usr/bin/env bash
# Follow-up one-axis sweeps at the main-grid winner (k=4, r=2):
#   HAMMING_MISMATCH=2
#   MIN_REFINE_HITS=5, MIN_REFINE_HITS=20
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
BUILD_DIR="$HERE/results/builds"
TSV_DIR="$HERE/results/tsv"
TIMING="$HERE/results/timings_followup.tsv"

mkdir -p "$BUILD_DIR" "$TSV_DIR"
: > "$TIMING"
echo -e "dataset\tvariant\twall_s" >> "$TIMING"

# Build objects for non-yamme are already cached from the main sweep.
NON_YAMME_OBJS=()
for src in yamdedup yamenr yamscan yamshuf yamtk; do
  NON_YAMME_OBJS+=("$BUILD_DIR/${src}.o")
done

DATASETS=(D1 D2 D3 D4 D5 D6)

build_variant() {
  local label="$1"; shift
  local me_obj="$BUILD_DIR/yamme_${label}.o"
  local bin="$BUILD_DIR/yamtk_${label}"
  echo "[build] ${label}"
  cc -std=gnu99 -O2 -DTOP_K_SEEDS=4 -DREFINE_PASSES=2 "$@" \
     -c "$REPO/src/yamme.c" -o "$me_obj"
  cc -std=gnu99 -O2 "${NON_YAMME_OBJS[@]}" "$me_obj" \
     -o "$bin" -lz -lm -pthread
  echo "$bin"
}

run_variant() {
  local label="$1"
  local bin="$BUILD_DIR/yamtk_${label}"
  for d in "${DATASETS[@]}"; do
    tsv="$TSV_DIR/${d}_${label}.tsv"
    TIMEFORMAT='%R'
    wall=$({ time "$bin" me -i "$FIX_DIR/${d}_pos.fa" \
      -k 6 -K 15 -N 5 -s 42 -o "$tsv" -M '' > /dev/null 2>&1; } 2>&1)
    echo -e "${d}\t${label}\t${wall}" >> "$TIMING"
    echo "  [run] $d  ${label}  ${wall}s"
  done
}

build_variant hm2_k4_r2  -DHAMMING_MISMATCH=2
run_variant   hm2_k4_r2

build_variant mrh5_k4_r2  -DMIN_REFINE_HITS=5
run_variant   mrh5_k4_r2

build_variant mrh20_k4_r2  -DMIN_REFINE_HITS=20
run_variant   mrh20_k4_r2

echo "Done.  Follow-up timings: $TIMING"
