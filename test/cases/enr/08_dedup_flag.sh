#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

# Default (no -d): duplicate motif name should abort with an error message.
assert_exit_nonzero "enr aborts on duplicate motif name by default" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_dup_motif.meme" \
        -t 5e-4

assert_stderr_contains "enr stderr mentions 'duplicate motif name'" \
    "duplicate motif name" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_dup_motif.meme" \
        -t 5e-4

# With -d: exits 0 and emits two rows; second name gets __N2 suffix.
assert_exit0 "enr -d deduplicates and exits 0" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_dup_motif.meme" \
        -t 5e-4 \
        -d

out=$("$YAMTK" enr \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_dup_motif.meme" \
    -t 5e-4 -d 2>/dev/null | grep -v '^#')

n2_row=$(echo "$out" | grep '__N2')
if [ -z "$n2_row" ]; then
    FAIL "enr -d: second duplicate row should have __N2 suffix"
fi
PASS "enr -d: duplicate row has __N2 suffix"

n_rows=$(echo "$out" | grep -c '^ebox')
if [ "$n_rows" -ne 2 ]; then
    echo "Expected 2 ebox rows, got $n_rows" >&2
    FAIL "enr -d: two output rows for two deduped motifs"
fi
PASS "enr -d: two output rows for two deduped motifs"
