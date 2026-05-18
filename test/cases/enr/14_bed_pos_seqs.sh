#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -x with -T seqs: each BED region is one unit; pos_n must reflect the
# region count, and the implant motif's seq_hits should equal it (50/50).
assert_enr_golden "enr -x positives (seqs)" \
    "$TESTDIR/expected/enr_bed_pos_seqs.txt" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 \
        -x "$TESTDIR/fixtures/enr_pos.bed"

# pos_n column (index 4 after consensus inserted at 3) must be the BED region count.
pos_n=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -x "$TESTDIR/fixtures/enr_pos.bed" 2>/dev/null \
    | grep -v '^#' | awk '$1=="ebox"{print $4; exit}')
if [ "$pos_n" -ne 50 ]; then
    FAIL "pos_n under -x equals BED region count (got $pos_n, expected 50)"
fi
PASS "pos_n under -x equals BED region count"
