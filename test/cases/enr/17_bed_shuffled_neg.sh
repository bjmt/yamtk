#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -x without -n: shuffled negatives are built from BED slices of positives,
# one shuffled seq per region. neg_n should match the BED region count.
assert_enr_golden "enr -x shuffled (seed=42 k=2)" \
    "$TESTDIR/expected/enr_bed_shuffled.txt" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 -s 42 -k 2 \
        -x "$TESTDIR/fixtures/enr_pos.bed"

# neg_n (col 7 after consensus inserted at 3) must equal BED region count.
neg_n=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -s 42 -k 2 -x "$TESTDIR/fixtures/enr_pos.bed" 2>/dev/null \
    | grep -v '^#' | awk '$1=="ebox"{print $7; exit}')
if [ "$neg_n" -ne 50 ]; then
    FAIL "neg_n under shuffled+-x equals BED region count (got $neg_n, expected 50)"
fi
PASS "neg_n under shuffled+-x equals BED region count"
