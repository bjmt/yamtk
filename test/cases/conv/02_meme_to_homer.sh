#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmp=$(mktemp /tmp/yamconv_test.XXXXXX)
trap 'rm -f "$tmp"' EXIT
"$YAMTK" conv -m "$TESTDIR/fixtures/enr_motifs.meme" -t homer > "$tmp" 2>/dev/null
if ! diff -u "$TESTDIR/expected/conv_meme_to_homer.out" "$tmp" >/dev/null; then
    diff -u "$TESTDIR/expected/conv_meme_to_homer.out" "$tmp" >&2 || true
    FAIL "MEME -> HOMER output differs from golden"
fi
PASS "MEME -> HOMER matches golden"
