#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf -x k=2 (BED ranges, Euler)" \
    "$TESTDIR/expected/shuf_bed_euler.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -x "$TESTDIR/dna.bed" -k 2 -s 42
