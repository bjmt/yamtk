#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seq stats TSV matches golden" \
    "$TESTDIR/expected/seq_stats.txt" \
    "$YAMTK" seq -a stats -i "$TESTDIR/dna.fa"
