#!/usr/bin/env bash
# Bug #5: dedup eviction used last.end instead of max(end) in the bucket.
# Input has hits [1,1000],[2,3],[4,5],[6,500] on the same (seq,motif,strand).
# [1,1000] has the lowest p-value (highest priority) and overlaps all others;
# with the fix only [1,1000] survives.  Before the fix, eviction happened
# prematurely at [4,5] because last.end=3 < start=4, and [6,500] escaped.
source "$TESTDIR/lib.sh"

out=$("$YAMTK" dedup -i "$TESTDIR/fixtures/overlapping_bucket.txt" 2>/dev/null)
count=$(echo "$out" | grep -c '^s1' || true)

if [ "$count" -ne 1 ]; then
    echo "not ok - expected 1 survivor (the [1,1000] hit), got $count" >&2
    echo "Output:" >&2; echo "$out" >&2; exit 1
fi
winner_end=$(echo "$out" | grep '^s1' | awk '{print $3}')
if [ "$winner_end" != "1000" ]; then
    echo "not ok - expected winner end=1000, got $winner_end" >&2; exit 1
fi
echo "ok - dedup eviction uses max_end; [1,1000] correctly dominates [6,500]"
