#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seq mask hard (-N) matches golden" \
    "$TESTDIR/expected/seq_mask_hard.fa" \
    "$YAMTK" seq -a mask -N -x "$TESTDIR/dna.bed" -i "$TESTDIR/dna.fa"
