#!/usr/bin/env bash
# Verify which (TOP_K_SEEDS, MIN_REFINE_HITS) preset gives best PWM fidelity
# on degenerate fixtures (D9–D11), to choose the "ambiguous mode" defaults.
# REFINE_PASSES stays at 2 (the sweep showed it has marginal effect on PWM).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
BUILD_DIR="$HERE/results/builds"
TSV_DIR="$HERE/results/tsv"
TIMING="$HERE/results/timings_degen_verify.tsv"

mkdir -p "$BUILD_DIR" "$TSV_DIR"
: > "$TIMING"
echo -e "dataset\tvariant\twall_s" >> "$TIMING"

NON_YAMME_OBJS=()
for src in yamdedup yamenr yamscan yamshuf yamtk; do
  NON_YAMME_OBJS+=("$BUILD_DIR/${src}.o")
done

DATASETS=(D9 D10 D11)

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

# Recommended baseline (for reference)
build_variant degen_baseline   -DTOP_K_SEEDS=4  -DREFINE_PASSES=2 -DMIN_REFINE_HITS=20
run_variant   degen_baseline

# Candidate 1: mirror of old default
build_variant degen_k20_mrh10  -DTOP_K_SEEDS=20 -DREFINE_PASSES=2 -DMIN_REFINE_HITS=10
run_variant   degen_k20_mrh10

# Candidate 2: even looser hits floor
build_variant degen_k20_mrh5   -DTOP_K_SEEDS=20 -DREFINE_PASSES=2 -DMIN_REFINE_HITS=5
run_variant   degen_k20_mrh5

# Candidate 3: more seeds
build_variant degen_k40_mrh10  -DTOP_K_SEEDS=40 -DREFINE_PASSES=2 -DMIN_REFINE_HITS=10
run_variant   degen_k40_mrh10

# Candidate 4: more seeds + lower floor + a bit more refinement
build_variant degen_k40_r5_mrh5  -DTOP_K_SEEDS=40 -DREFINE_PASSES=5 -DMIN_REFINE_HITS=5
run_variant   degen_k40_r5_mrh5

echo "Done.  Timings: $TIMING"
