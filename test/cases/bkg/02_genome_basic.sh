#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

bedfile=$(mktemp -t yamtkbkgbasicbed.XXXXXX)
trap 'rm -f "$bedfile"' EXIT

assert_golden "bkg genome basic FASTA matches golden" \
    "$TESTDIR/expected/bkg_genome_basic.fa" \
    "$YAMTK" bkg -i "$TESTDIR/fixtures/enr_pos.fa" -g "$TESTDIR/fixtures/enr_neg.fa" \
        -s 1 -G 0.1 -B "$bedfile"

if ! diff -u "$TESTDIR/expected/bkg_genome_basic.bed" "$bedfile" >/dev/null 2>&1; then
    diff -u "$TESTDIR/expected/bkg_genome_basic.bed" "$bedfile" >&2 || true
    FAIL "bkg genome basic BED differs from golden"
fi
PASS "bkg genome basic BED matches golden"
