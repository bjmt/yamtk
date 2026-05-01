#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "scan MEME motif" \
    "$TESTDIR/expected/scan_meme.txt" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 1
