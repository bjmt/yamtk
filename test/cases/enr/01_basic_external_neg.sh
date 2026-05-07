#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_enr_golden "enr basic external neg" \
    "$TESTDIR/expected/enr_basic.txt" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4
