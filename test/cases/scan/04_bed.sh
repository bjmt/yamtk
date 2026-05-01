#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "scan with BED restriction" \
    "$TESTDIR/expected/scan_bed.txt" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" \
        -x "$TESTDIR/dna.bed" -t 0.9 -j 1
