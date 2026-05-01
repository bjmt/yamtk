#!/usr/bin/env bash
# Bug #4: score_seq_in_bed used '<' instead of '<=' for the final window,
# silently skipping the hit when bed_size == motif_size.
# BED [2,7] (+) is exactly 5 bp wide = motif width; the sole valid hit must appear.
source "$TESTDIR/lib.sh"

out=$("$YAMTK" scan -m "$TESTDIR/motif.meme" \
    -s "$TESTDIR/dna.fa" \
    -x "$TESTDIR/fixtures/exact_width.bed" \
    -t 0.9 -j 1 2>/dev/null)

# Expect the hit at position 3-7 on seq "1" strand "+"
# BED scan output has extra columns (bed_range, bed_name) before seq fields
if ! echo "$out" | grep -qP '\t3\t7\t\+' 2>/dev/null \
        && ! echo "$out" | awk -F'\t' '$4==3&&$5==7&&$6=="+"' | grep -q .; then
    echo "not ok - missing hit at 1:3-7(+) in exact-width BED scan" >&2
    echo "Output was:" >&2; echo "$out" >&2; exit 1
fi
echo "ok - exact-width BED interval produces final-position hit"
