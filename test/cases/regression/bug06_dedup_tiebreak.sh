#!/usr/bin/env bash
# Bug #6: compare_scores returned 0 on ties, making qsort non-deterministic.
# Same input run twice must produce byte-identical output.
source "$TESTDIR/lib.sh"

run1=$("$YAMTK" dedup -i "$TESTDIR/fixtures/tied_scores.txt" 2>/dev/null | grep -v '^##yamscan ')
run2=$("$YAMTK" dedup -i "$TESTDIR/fixtures/tied_scores.txt" 2>/dev/null | grep -v '^##yamscan ')

if [ "$run1" != "$run2" ]; then
    diff <(echo "$run1") <(echo "$run2") >&2
    echo "not ok - dedup tiebreak is non-deterministic across runs" >&2; exit 1
fi
echo "ok - dedup tiebreak is deterministic (identical output across two runs)"
