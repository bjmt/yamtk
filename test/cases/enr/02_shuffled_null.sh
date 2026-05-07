#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_enr_golden "enr shuffled null k=2 seed=42" \
    "$TESTDIR/expected/enr_shuffled.txt" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 \
        -s 42 \
        -k 2
