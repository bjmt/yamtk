#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seq mask soft (default lowercase) matches golden" \
    "$TESTDIR/expected/seq_mask_soft.fa" \
    "$YAMTK" seq -a mask -x "$TESTDIR/dna.bed" -i "$TESTDIR/dna.fa"
