#!/usr/bin/env bash
#
# yamref_param_sweep.sh — sweep yamref refinement parameters over a grid
# and report recovered consensus + IC + hits + runtime per cell.
#
# Designed to inform whether yamref's compile-time defaults
# (DEFAULT_HIT_PVAL, DEFAULT_R_PASSES, MIN_IC_BITS) are well-chosen.
#
# Benchmark: refine a 6-bp narrow seed (CGTACG) against an implanted
# 8-bp truth (ACGTACGT) in test/fixtures/me_implanted.fa. Always uses
# auto-extend (-E) so all three swept parameters are exercised.
#
# Usage:
#   bash scripts/yamref_param_sweep.sh [yamtk-binary]
#
# Defaults to ./yamtk if no argument given. Run from repo root.

set -euo pipefail

YAMTK="${1:-./yamtk}"
SEED="test/fixtures/ref_seed_narrow.meme"
INPUT="test/fixtures/me_implanted.fa"
TRUTH="ACGTACGT"

if [ ! -x "$YAMTK" ]; then
    echo "error: yamtk binary not found at '$YAMTK'" >&2
    echo "       build first with 'make yamtk' or pass an explicit path." >&2
    exit 1
fi
if [ ! -f "$SEED" ] || [ ! -f "$INPUT" ]; then
    echo "error: fixtures not found; run from the repo root." >&2
    exit 1
fi

# Grid
T_VALUES=(1e-2 1e-3 5e-4 1e-4 5e-5 1e-5)
N_VALUES=(1 2 3 5)
I_VALUES=(0.25 0.5 0.75 1.0)

ms_now() {
    if date +%s%3N 2>/dev/null | grep -qE '^[0-9]+$'; then
        date +%s%3N
    else
        # macOS / BSD date lacks %3N; fall back to seconds and python
        python3 -c 'import time; print(int(time.time()*1000))'
    fi
}

# Parse stderr summary line:
#   [name] w: W0->W | IC: IC0->IC bits | hits: H0->H | consensus=CONS
parse_summary() {
    local stderr="$1"
    local field="$2"
    case "$field" in
        width)     echo "$stderr" | sed -n 's/.* w: [0-9]*->\([0-9]*\) .*/\1/p' | head -1 ;;
        ic)        echo "$stderr" | sed -n 's/.* IC: [0-9.]*->\([0-9.]*\) bits .*/\1/p' | head -1 ;;
        hits)      echo "$stderr" | sed -n 's/.* hits: [0-9]*->\([0-9]*\) .*/\1/p' | head -1 ;;
        consensus) echo "$stderr" | sed -n 's/.*consensus=\([A-Z]*\).*/\1/p' | head -1 ;;
    esac
}

printf 't\tn\tT\twidth\tconsensus\tic_bits\thits\truntime_ms\trecovered\n'

for t in "${T_VALUES[@]}"; do
    for n in "${N_VALUES[@]}"; do
        for i in "${I_VALUES[@]}"; do
            t0=$(ms_now)
            stderr=$("$YAMTK" ref -m "$SEED" -i "$INPUT" -E \
                              -t "$t" -n "$n" -T "$i" -o /dev/null 2>&1 || true)
            t1=$(ms_now)
            dt=$((t1 - t0))
            width=$(parse_summary "$stderr" width)
            ic=$(parse_summary "$stderr" ic)
            hits=$(parse_summary "$stderr" hits)
            consensus=$(parse_summary "$stderr" consensus)
            [ -z "$width" ]     && width="-"
            [ -z "$ic" ]        && ic="-"
            [ -z "$hits" ]      && hits="-"
            [ -z "$consensus" ] && consensus="-"
            if [ "$consensus" = "$TRUTH" ]; then recovered=1; else recovered=0; fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$t" "$n" "$i" "$width" "$consensus" "$ic" "$hits" "$dt" "$recovered"
        done
    done
done
