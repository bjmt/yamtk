#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_enr_golden "enr -T ranksum mode" \
    "$TESTDIR/expected/enr_ranksum.txt" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 \
        -T ranksum

# Verify count columns are identical to seqs-mode output (counts are mode-independent).
seqs_out=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 2>/dev/null | grep -v '^#')

ranksum_out=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -T ranksum 2>/dev/null | grep -v '^#')

seqs_counts=$(echo "$seqs_out"    | awk '{print $4,$5,$6,$7,$8,$9}')
ranks_counts=$(echo "$ranksum_out" | awk '{print $4,$5,$6,$7,$8,$9}')
if [ "$seqs_counts" != "$ranks_counts" ]; then
    echo "Count columns differ between seqs and ranksum modes:" >&2
    diff <(echo "$seqs_counts") <(echo "$ranks_counts") >&2 || true
    FAIL "enr ranksum: count columns identical to seqs mode"
fi
PASS "enr ranksum: count columns identical to seqs mode"

# Spiked motif (ebox) should have AUC > 0.7
ebox_auc=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -T ranksum 2>/dev/null \
    | awk '$1=="ebox"{print $10}')
ebox_ok=$(echo "$ebox_auc" | awk '{print ($1+0 > 0.7) ? "yes" : "no"}')
if [ "$ebox_ok" != "yes" ]; then
    echo "ebox AUC=$ebox_auc, expected > 0.7" >&2
    FAIL "enr ranksum: spiked motif AUC > 0.7"
fi
PASS "enr ranksum: spiked motif AUC > 0.7"
