#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmpdir=$(mktemp -d /tmp/yamseed_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
              -i "$TESTDIR/fixtures/me_implanted.fa" \
              -f 0.01 -s 1 -O "$tmpdir/truth.bed" >/dev/null 2>&1
if ! diff -u "$TESTDIR/expected/seed_random_s1_truth.bed" "$tmpdir/truth.bed" >/dev/null 2>&1; then
    diff -u "$TESTDIR/expected/seed_random_s1_truth.bed" "$tmpdir/truth.bed" >&2 || true
    FAIL "seed random truth BED differs from golden"
fi
PASS "seed random truth BED matches golden"
