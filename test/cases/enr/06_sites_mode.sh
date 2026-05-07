#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_enr_golden "enr -T sites mode" \
    "$TESTDIR/expected/enr_sites.txt" \
    "$YAMTK" enr \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 \
        -T sites

# Verify seq_hits and site_hits columns are identical to seqs-mode output
# (counts don't depend on test mode, only effect/pvalue/qvalue differ).
seqs_out=$("$YAMTK" enr \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 2>/dev/null | grep -v '^#')

sites_out=$("$YAMTK" enr \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -T sites 2>/dev/null | grep -v '^#')

# Columns 3-8 (counts) should be identical between modes
seqs_counts=$(echo "$seqs_out"  | awk '{print $3,$4,$5,$6,$7,$8}')
sites_counts=$(echo "$sites_out" | awk '{print $3,$4,$5,$6,$7,$8}')
if [ "$seqs_counts" != "$sites_counts" ]; then
    echo "Count columns differ between seqs and sites modes:" >&2
    diff <(echo "$seqs_counts") <(echo "$sites_counts") >&2 || true
    FAIL "enr sites: count columns identical to seqs mode"
fi
PASS "enr sites: count columns identical to seqs mode"
