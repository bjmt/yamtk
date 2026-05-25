#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmp=$(mktemp /tmp/yamconv_test.XXXXXX)
trap 'rm -f "$tmp"' EXIT
"$YAMTK" conv -m "$TESTDIR/fixtures/enr_motifs.meme" -t meme > "$tmp" 2>/dev/null
if ! diff -u "$TESTDIR/expected/conv_meme_to_meme.out" "$tmp" >/dev/null; then
    diff -u "$TESTDIR/expected/conv_meme_to_meme.out" "$tmp" >&2 || true
    FAIL "MEME -> MEME pass-through differs from golden"
fi
PASS "MEME -> MEME pass-through matches golden"
