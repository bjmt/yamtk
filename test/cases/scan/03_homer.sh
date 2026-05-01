#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "scan HOMER motif" \
    "$TESTDIR/expected/scan_homer.txt" \
    "$YAMTK" scan -m "$TESTDIR/motif.homer" -s "$TESTDIR/dna.fa" -t 0.05 -j 1
