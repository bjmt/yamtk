#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seq hist default-step TSV matches golden" \
    "$TESTDIR/expected/seq_hist_default.txt" \
    "$YAMTK" seq -a hist -i "$TESTDIR/dna.fa"
