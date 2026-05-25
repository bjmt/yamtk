#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmp=$(mktemp /tmp/yamconv_test.XXXXXX)
trap 'rm -f "$tmp"' EXIT
"$YAMTK" conv -m "$TESTDIR/fixtures/enr_motifs.meme" -t hocomoco > "$tmp" 2>/dev/null
if ! diff -u "$TESTDIR/expected/conv_meme_to_hocomoco.out" "$tmp" >/dev/null; then
    diff -u "$TESTDIR/expected/conv_meme_to_hocomoco.out" "$tmp" >&2 || true
    FAIL "MEME -> HOCOMOCO output differs from golden"
fi
PASS "MEME -> HOCOMOCO matches golden"
