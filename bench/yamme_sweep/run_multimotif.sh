#!/usr/bin/env bash
# Focused multi-motif comparison on D7 (3 motifs of decreasing strength)
# and D8 (4 moderate motifs). Compares:
#   current defaults: k=20, r=2, MRH=10
#   recommended:      k=4,  r=2, MRH=20
#   also: k=10, r=2, MRH=20 (a hedge in case k=4 is too tight here)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
BUILD_DIR="$HERE/results/builds"
TSV_DIR="$HERE/results/tsv"
TIMING="$HERE/results/timings_multimotif.tsv"

mkdir -p "$BUILD_DIR" "$TSV_DIR"
: > "$TIMING"
echo -e "dataset\tvariant\twall_s" >> "$TIMING"

NON_YAMME_OBJS=()
for src in yamdedup yamenr yamscan yamshuf yamtk; do
  NON_YAMME_OBJS+=("$BUILD_DIR/${src}.o")
done

DATASETS=(D7 D8 D9 D10 D11)

build_variant() {
  local label="$1"; shift
  local me_obj="$BUILD_DIR/yamme_${label}.o"
  local bin="$BUILD_DIR/yamtk_${label}"
  echo "[build] ${label}"
  cc -std=gnu99 -O2 "$@" -c "$REPO/src/yamme.c" -o "$me_obj"
  cc -std=gnu99 -O2 "${NON_YAMME_OBJS[@]}" "$me_obj" \
     -o "$bin" -lz -lm -pthread
}

run_variant() {
  local label="$1"
  local bin="$BUILD_DIR/yamtk_${label}"
  for d in "${DATASETS[@]}"; do
    tsv="$TSV_DIR/${d}_${label}.tsv"
    meme="$TSV_DIR/${d}_${label}.meme"
    TIMEFORMAT='%R'
    wall=$({ time "$bin" me -i "$FIX_DIR/${d}_pos.fa" \
      -k 6 -K 15 -N 8 -s 42 -o "$tsv" -M "$meme" > /dev/null 2>&1; } 2>&1)
    echo -e "${d}\t${label}\t${wall}" >> "$TIMING"
    echo "  [run] $d  ${label}  ${wall}s"
  done
}

# Current default
build_variant default      -DTOP_K_SEEDS=20 -DREFINE_PASSES=2  -DMIN_REFINE_HITS=10
run_variant   default

# Recommended (k=4, r=2, MRH=20)
build_variant recommended  -DTOP_K_SEEDS=4  -DREFINE_PASSES=2  -DMIN_REFINE_HITS=20
run_variant   recommended

# Hedge: k=10
build_variant hedge_k10    -DTOP_K_SEEDS=10 -DREFINE_PASSES=2  -DMIN_REFINE_HITS=20
run_variant   hedge_k10

echo "Done.  Timings: $TIMING"
