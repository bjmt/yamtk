#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

# Fixed-count mode: deterministic golden FASTA + truth BED at -n 3 -s 1.
assert_golden "seed fixed (-n 3 -s 1) FASTA matches golden" \
    "$TESTDIR/expected/seed_fixed_n3_s1.fa" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
                  -i "$TESTDIR/fixtures/me_implanted.fa" \
                  -n 3 -s 1

tmpdir=$(mktemp -d /tmp/yamseed_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
              -i "$TESTDIR/fixtures/me_implanted.fa" \
              -n 3 -s 1 -O "$tmpdir/truth.bed" >/dev/null 2>&1
if ! diff -u "$TESTDIR/expected/seed_fixed_n3_s1_truth.bed" "$tmpdir/truth.bed" >/dev/null 2>&1; then
    diff -u "$TESTDIR/expected/seed_fixed_n3_s1_truth.bed" "$tmpdir/truth.bed" >&2 || true
    FAIL "seed fixed truth BED differs from golden"
fi
PASS "seed fixed truth BED matches golden"

# Every sequence should have exactly 3 insertions (all fixture seqs are L=100,
# motifs are 6 bp -> ample capacity).
uniq_counts=$(awk '{print $1}' "$tmpdir/truth.bed" | sort | uniq -c | awk '{print $1}' | sort -u)
if [ "$uniq_counts" != "3" ]; then
    FAIL "seed fixed: expected exactly 3 insertions per sequence, got counts: $uniq_counts"
fi
PASS "seed fixed produces exactly N insertions per sequence"

# -n + -f mutual exclusion must fail.
assert_exit_nonzero "seed -n with -f rejected" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
                  -i "$TESTDIR/fixtures/me_implanted.fa" \
                  -n 3 -f 0.01 -s 1
