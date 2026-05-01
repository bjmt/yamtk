#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

actual=$("$YAMTK" dedup -i "$TESTDIR/fixtures/scan_output_for_dedup.txt" 2>/dev/null | grep -v '^##yamscan ')

if ! diff -u "$TESTDIR/expected/dedup_basic.txt" <(echo "$actual") >/dev/null 2>&1; then
    diff -u "$TESTDIR/expected/dedup_basic.txt" <(echo "$actual") >&2 || true
    echo "not ok - dedup output differs from golden" >&2; exit 1
fi
echo "ok - dedup basic matches golden"
