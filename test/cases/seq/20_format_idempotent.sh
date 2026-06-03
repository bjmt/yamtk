#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# format is idempotent: running it twice gives the same output as once.
tmpdir=$(mktemp -d /tmp/yamseq_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seq -a format -i "$TESTDIR/dna.fa" 2>/dev/null > "$tmpdir/once.fa"
"$YAMTK" seq -a format -i "$tmpdir/once.fa" 2>/dev/null > "$tmpdir/twice.fa"
if ! diff -u "$tmpdir/once.fa" "$tmpdir/twice.fa" >/dev/null 2>&1; then
    diff -u "$tmpdir/once.fa" "$tmpdir/twice.fa" >&2 || true
    FAIL "format format differs from format"
fi
PASS "format is idempotent"
