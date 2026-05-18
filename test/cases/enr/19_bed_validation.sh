#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
TMP_BAD=$(mktemp -t yamenr_bed_XXXXXX)
TMP_MISS=$(mktemp -t yamenr_bed_XXXXXX)
TMP_OOB=$(mktemp -t yamenr_bed_XXXXXX)
trap 'rm -f "$TMP_BAD" "$TMP_MISS" "$TMP_OOB"' EXIT

# -X without -n
assert_exit_nonzero "enr -X without -n fails" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -X "$TESTDIR/fixtures/enr_neg.bed"

# -X without -x
assert_exit_nonzero "enr -X without -x fails" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -X "$TESTDIR/fixtures/enr_neg.bed"

# BED with start >= end
printf "pos1\t50\t50\n" > "$TMP_BAD"
assert_exit_nonzero "enr -x with bad BED (start>=end) fails" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -x "$TMP_BAD"

# BED references unknown sequence
printf "unknown_seq\t10\t20\n" > "$TMP_MISS"
assert_exit_nonzero "enr -x with unknown BED seq fails" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -x "$TMP_MISS"

# BED end past seq length
printf "pos1\t0\t10000\n" > "$TMP_OOB"
assert_exit_nonzero "enr -x with out-of-bounds BED fails" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -x "$TMP_OOB"

# -l + -x without -n (shuffled-from-BED is not streamable)
assert_exit_nonzero "enr -l + -x without -n fails" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -x "$TESTDIR/fixtures/enr_pos.bed" -l
