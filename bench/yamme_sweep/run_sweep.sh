#!/usr/bin/env bash
# Sweep yamme over a (TOP_K_SEEDS, REFINE_PASSES) grid across all D{1..6}
# fixtures. Builds one yamtk binary per grid point with -D overrides,
# then runs it across all fixtures and saves TSVs + timing.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
FIX_DIR="$HERE/fixtures"
BUILD_DIR="$HERE/results/builds"
TSV_DIR="$HERE/results/tsv"
TIMING="$HERE/results/timings.tsv"

mkdir -p "$BUILD_DIR" "$TSV_DIR"
: > "$TIMING"
echo -e "dataset\ttop_k\trefine\twall_s" >> "$TIMING"

# Compile the non-yamme objects once, reuse for all grid points.
NON_YAMME_OBJS=()
for src in yamdedup yamenr yamscan yamshuf yamtk; do
  obj="$BUILD_DIR/${src}.o"
  if [ ! -f "$obj" ] || [ "$REPO/src/${src}.c" -nt "$obj" ]; then
    cc -std=gnu99 -O2 -c "$REPO/src/${src}.c" -o "$obj"
  fi
  NON_YAMME_OBJS+=("$obj")
done

KS=(4 10 20 40)
RS=(2 5 10 20)
DATASETS=(D1 D2 D3 D4 D5 D6)

for k in "${KS[@]}"; do
  for r in "${RS[@]}"; do
    bin="$BUILD_DIR/yamtk_k${k}_r${r}"
    me_obj="$BUILD_DIR/yamme_k${k}_r${r}.o"
    echo "[build] k=$k r=$r"
    cc -std=gnu99 -O2 -DTOP_K_SEEDS=$k -DREFINE_PASSES=$r \
       -c "$REPO/src/yamme.c" -o "$me_obj"
    cc -std=gnu99 -O2 "${NON_YAMME_OBJS[@]}" "$me_obj" \
       -o "$bin" -lz -lm -pthread

    for d in "${DATASETS[@]}"; do
      tsv="$TSV_DIR/${d}_k${k}_r${r}.tsv"
      # Use the shell's builtin `time` via `TIMEFORMAT` to get plain seconds.
      TIMEFORMAT='%R'
      wall=$({ time "$bin" me -i "$FIX_DIR/${d}_pos.fa" \
        -k 6 -K 15 -N 5 -s 42 -o "$tsv" -M '' > /dev/null 2>&1; } 2>&1)
      echo -e "${d}\t${k}\t${r}\t${wall}" >> "$TIMING"
      echo "  [run] $d  k=$k  r=$r  ${wall}s"
    done
  done
done

echo "Done.  Timings: $TIMING"
echo "TSVs:    $TSV_DIR"
