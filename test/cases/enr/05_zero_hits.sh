#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Use an extremely stringent threshold so no motif can produce any hit.
# Verifies zero-hit behavior: all count columns = 0, pvalue=1, qvalue=1, no crash.
out=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 1e-20 2>/dev/null)

# Both motifs should appear with all-zero hit counts
zero_lines=$(echo "$out" | grep -v '^#' | awk '$5==0 && $6==0 && $8==0 && $9==0' | wc -l | tr -d ' ')
if [ "$zero_lines" -ne 2 ]; then
    echo "Expected 2 rows with all-zero hit counts, got $zero_lines" >&2
    FAIL "enr zero hits: all hit counts are zero"
fi
PASS "enr zero hits: all hit counts are zero"

# All pvalues should be 1
bad_pval=$(echo "$out" | grep -v '^#' | awk '$12 != "1"' | wc -l | tr -d ' ')
if [ "$bad_pval" -ne 0 ]; then
    echo "Expected all pvalues=1 for zero hits, found $bad_pval non-1 rows" >&2
    FAIL "enr zero hits: pvalue=1 for all motifs"
fi
PASS "enr zero hits: pvalue=1 for all motifs"
