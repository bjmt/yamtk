#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -x and -X together: both sets are BED-restricted.
assert_enr_golden "enr -x -X" \
    "$TESTDIR/expected/enr_bed_pos_neg.txt" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 \
        -x "$TESTDIR/fixtures/enr_pos.bed" \
        -X "$TESTDIR/fixtures/enr_neg.bed"
