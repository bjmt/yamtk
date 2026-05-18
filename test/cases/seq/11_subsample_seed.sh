#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Same seed -> identical output, every time.
a=$("$YAMTK" seq -a subsample -n 2 -s 42 -i "$TESTDIR/dna.fa" 2>/dev/null)
b=$("$YAMTK" seq -a subsample -n 2 -s 42 -i "$TESTDIR/dna.fa" 2>/dev/null)
if [ "$a" != "$b" ]; then
    FAIL "subsample with same seed is deterministic"
fi
PASS "subsample with same seed is deterministic"

# Two seeds should at minimum sometimes disagree. With only 3 input records and
# 2 kept, every seed produces one of 3 outcomes ({1,2}, {1,3}, {2,3}). Across
# a handful of seeds we should see at least 2 distinct outputs.
distinct=$(for s in 1 2 3 4 5 6 7 8 9 10; do
    "$YAMTK" seq -a subsample -n 2 -s "$s" -i "$TESTDIR/dna.fa" 2>/dev/null \
        | awk '/^>/ {printf "%s,", $1} END {print ""}'
done | sort -u | wc -l | tr -d ' ')
if [ "$distinct" -lt 2 ]; then
    FAIL "subsample varies with seed (saw $distinct distinct outputs across 10 seeds)"
fi
PASS "subsample varies with seed"
