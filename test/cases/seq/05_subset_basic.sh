#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seq subset basic FASTA matches golden" \
    "$TESTDIR/expected/seq_subset_basic.fa" \
    "$YAMTK" seq -a subset -x "$TESTDIR/dna.bed" -i "$TESTDIR/dna.fa"
