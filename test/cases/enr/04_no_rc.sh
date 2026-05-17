#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_enr_golden "enr forward-only -R" \
    "$TESTDIR/expected/enr_fwd_only.txt" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 \
        -R
