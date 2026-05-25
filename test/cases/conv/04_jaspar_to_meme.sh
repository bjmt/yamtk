#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmp=$(mktemp /tmp/yamconv_test.XXXXXX)
trap 'rm -f "$tmp"' EXIT
"$YAMTK" conv -m "$TESTDIR/fixtures/ref_seed.jaspar" -t meme > "$tmp" 2>/dev/null
if ! diff -u "$TESTDIR/expected/conv_jaspar_to_meme.out" "$tmp" >/dev/null; then
    diff -u "$TESTDIR/expected/conv_jaspar_to_meme.out" "$tmp" >&2 || true
    FAIL "JASPAR -> MEME output differs from golden"
fi
PASS "JASPAR -> MEME matches golden"
