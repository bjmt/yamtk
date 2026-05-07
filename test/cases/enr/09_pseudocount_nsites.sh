#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

# Verify -p and -N are accepted and produce valid output.
assert_exit0 "enr -p 10 exits 0" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 -p 10

assert_exit0 "enr -N 100 exits 0" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 -N 100

# Parity check: pos_site_hits from enr -p 10 must match yamtk scan -p 10.
enr_hits=$("$YAMTK" enr \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -p 10 2>/dev/null | awk '$1=="ebox"{print $5}')

scan_hits=$("$YAMTK" scan \
    -s "$TESTDIR/fixtures/enr_pos.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -p 10 -j 1 2>/dev/null | grep -v '^#' | awk '$5=="ebox"' | wc -l | tr -d ' ')

if [ "$enr_hits" -ne "$scan_hits" ]; then
    echo "enr -p 10 pos_site_hits=$enr_hits, scan -p 10 site_hits=$scan_hits" >&2
    FAIL "enr -p parity with yamtk scan"
fi
PASS "enr -p parity with yamtk scan"

# Reject invalid values.
assert_exit_nonzero "enr -p 0 exits non-zero" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 -p 0

assert_exit_nonzero "enr -N 0 exits non-zero" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 -N 0
