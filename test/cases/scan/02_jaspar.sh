#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "scan JASPAR motif" \
    "$TESTDIR/expected/scan_jaspar.txt" \
    "$YAMTK" scan -m "$TESTDIR/motif.jaspar" -s "$TESTDIR/dna.fa" -t 0.05 -j 1
