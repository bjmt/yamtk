#!/usr/bin/env bash
# Full verification of yamme presets against STREME on all degenerate fixtures.
#
# yamme presets:
#   default:     k=4,  r=2, MRH=20
#   proposed_A:  k=20, r=2, MRH=10   (proposed -A "ambiguous" mode)
#   hedge:       k=10, r=2, MRH=10   (intermediate)
#
# STREME: default settings, --minw 6 --maxw 15 --nmotifs 3 --seed 42
#
# Fixtures: D9..D11 (5 replicates each) + D12..D14 (big) + D15..D16 (stress).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
BUILD_DIR="$HERE/results/builds"
TSV_DIR="$HERE/results/tsv"
STREME_DIR="$HERE/results/streme"
TIMING="$HERE/results/timings_full_verify.tsv"

STREME_BIN="$HOME/meme/bin/streme"
if [ ! -x "$STREME_BIN" ]; then
  echo "ERROR: $STREME_BIN not found or not executable"
  exit 1
fi

mkdir -p "$BUILD_DIR" "$TSV_DIR" "$STREME_DIR"
: > "$TIMING"
echo -e "dataset\tvariant\twall_s" >> "$TIMING"

NON_YAMME_OBJS=()
for src in yamdedup yamenr yamscan yamshuf yamtk; do
  NON_YAMME_OBJS+=("$BUILD_DIR/${src}.o")
done

DATASETS=(
  D9   D9_v2  D9_v3  D9_v4  D9_v5
  D10  D10_v2 D10_v3 D10_v4 D10_v5
  D11  D11_v2 D11_v3 D11_v4 D11_v5
  D12 D13 D14
  D15 D16
)

build_variant() {
  local label="$1"; shift
  local me_obj="$BUILD_DIR/yamme_${label}.o"
  local bin="$BUILD_DIR/yamtk_${label}"
  echo "[build] ${label}"
  cc -std=gnu99 -O2 "$@" -c "$REPO/src/yamme.c" -o "$me_obj"
  cc -std=gnu99 -O2 "${NON_YAMME_OBJS[@]}" "$me_obj" \
     -o "$bin" -lz -lm -pthread
}

run_yamme() {
  local label="$1"
  local bin="$BUILD_DIR/yamtk_${label}"
  for d in "${DATASETS[@]}"; do
    tsv="$TSV_DIR/${d}_${label}.tsv"
    meme="$TSV_DIR/${d}_${label}.meme"
    TIMEFORMAT='%R'
    wall=$({ time "$bin" me -i "$FIX_DIR/${d}_pos.fa" \
      -k 6 -K 15 -N 8 -s 42 -o "$tsv" -M "$meme" > /dev/null 2>&1; } 2>&1)
    echo -e "${d}\t${label}\t${wall}" >> "$TIMING"
  done
  echo "[done] yamme ${label}"
}

run_streme() {
  for d in "${DATASETS[@]}"; do
    outdir="$STREME_DIR/${d}"
    rm -rf "$outdir"
    TIMEFORMAT='%R'
    wall=$({ time "$STREME_BIN" --p "$FIX_DIR/${d}_pos.fa" \
      --oc "$outdir" --minw 6 --maxw 15 --nmotifs 3 --seed 42 \
      --verbosity 1 > /dev/null 2>&1; } 2>&1)
    echo -e "${d}\tstreme\t${wall}" >> "$TIMING"
  done
  echo "[done] streme"
}

build_variant default      -DTOP_K_SEEDS=4  -DREFINE_PASSES=2 -DMIN_REFINE_HITS=20
run_yamme     default

build_variant proposed_A   -DTOP_K_SEEDS=20 -DREFINE_PASSES=2 -DMIN_REFINE_HITS=10
run_yamme     proposed_A

build_variant hedge        -DTOP_K_SEEDS=10 -DREFINE_PASSES=2 -DMIN_REFINE_HITS=10
run_yamme     hedge

run_streme

echo "Done.  Timings: $TIMING"
