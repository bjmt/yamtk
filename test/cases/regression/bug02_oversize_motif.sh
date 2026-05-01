#!/usr/bin/env bash
# Bug #2: oversize motif (> 50 cols) previously continued past the error log
# and wrote OOB into the PWM array. Now each parser calls badexit.
source "$TESTDIR/lib.sh"

for fmt in meme jaspar homer; do
    fixture="$TESTDIR/fixtures/oversize_$fmt.$fmt"
    if "$YAMTK" scan -m "$fixture" -s "$TESTDIR/dna.fa" -j 1 >/dev/null 2>&1; then
        echo "not ok - $fmt oversize motif: scan exited 0, expected non-zero" >&2; exit 1
    fi
    err=$("$YAMTK" scan -m "$fixture" -s "$TESTDIR/dna.fa" -j 1 2>&1 >/dev/null || true)
    if ! echo "$err" | grep -qi 'too large\|too wide\|too big\|exceeds\|maximum\|error'; then
        echo "not ok - $fmt oversize motif: no diagnostic in stderr" >&2
        echo "  stderr was: $err" >&2; exit 1
    fi
    echo "ok - $fmt oversize motif exits 1 with diagnostic"
done
