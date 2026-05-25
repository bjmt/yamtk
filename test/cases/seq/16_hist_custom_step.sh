#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seq hist -b 0.1 TSV matches golden" \
    "$TESTDIR/expected/seq_hist_step0p1.txt" \
    "$YAMTK" seq -a hist -b 0.1 -i "$TESTDIR/dna.fa"
