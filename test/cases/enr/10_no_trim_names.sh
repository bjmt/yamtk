#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

# Default: JASPAR motif name "ebox extra_info" trimmed to "ebox" at first whitespace.
trimmed=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_jaspar_spacedname.jaspar" \
    -t 5e-4 2>/dev/null | grep -v '^#' | awk -F'\t' '{print $1}')

if [ "$trimmed" != "ebox" ]; then
    echo "Default motif name: expected 'ebox', got '$trimmed'" >&2
    FAIL "enr default trims name to first word"
fi
PASS "enr default trims name to first word"

# With -r: full name including space is preserved.
full=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_jaspar_spacedname.jaspar" \
    -t 5e-4 -r 2>/dev/null | grep -v '^#' | awk -F'\t' '{print $1}')

if [ "$full" != "ebox extra_info" ]; then
    echo "With -r motif name: expected 'ebox extra_info', got '$full'" >&2
    FAIL "enr -r preserves full name with space"
fi
PASS "enr -r preserves full name with space"

# Verify -r -r (duplicate) exits non-zero.
assert_exit_nonzero "enr -r twice exits non-zero" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_jaspar_spacedname.jaspar" \
        -t 5e-4 -r -r
