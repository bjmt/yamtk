#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

bedfile=$(mktemp -t yamtkbkgstrandbed.XXXXXX)
trap 'rm -f "$bedfile"' EXIT

assert_golden "bkg genome strand FASTA matches golden" \
    "$TESTDIR/expected/bkg_genome_strand.fa" \
    "$YAMTK" bkg -i "$TESTDIR/fixtures/enr_pos.fa" -g "$TESTDIR/fixtures/enr_neg.fa" \
        -s 1 -G 0.1 -r -B "$bedfile"

if ! diff -u "$TESTDIR/expected/bkg_genome_strand.bed" "$bedfile" >/dev/null 2>&1; then
    diff -u "$TESTDIR/expected/bkg_genome_strand.bed" "$bedfile" >&2 || true
    FAIL "bkg genome strand BED differs from golden"
fi
PASS "bkg genome strand BED matches golden"

# Sanity: both '+' and '-' strands must appear.
plus=$(awk '$6 == "+"' "$bedfile" | wc -l | tr -d ' ')
minus=$(awk '$6 == "-"' "$bedfile" | wc -l | tr -d ' ')
if [ "$plus" -eq 0 ] || [ "$minus" -eq 0 ]; then
    FAIL "bkg genome strand: expected both strands; got +=$plus -=$minus"
fi
PASS "bkg genome strand: both '+' and '-' present (+=$plus -=$minus)"
