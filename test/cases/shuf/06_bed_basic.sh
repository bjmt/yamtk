#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf -x k=1 (BED ranges, Fisher-Yates)" \
    "$TESTDIR/expected/shuf_bed_k1.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -x "$TESTDIR/dna.bed" -k 1 -s 42
