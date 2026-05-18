#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Reservoir size larger than input -> output all input sequences, in order.
out=$("$YAMTK" seq -a subsample -n 100 -s 42 -i "$TESTDIR/dna.fa" 2>/dev/null)
n=$(echo "$out" | grep -c '^>')
if [ "$n" -ne 3 ]; then
    FAIL "subsample -n 100 on 3-seq input keeps all 3 (got $n)"
fi
PASS "subsample -n > input keeps all input"

names=$(echo "$out" | awk '/^>/ {print $1}' | tr -d '>' | tr '\n' ' ')
if [ "$names" != "1 2 3 " ]; then
    FAIL "subsample -n > input preserves order (got '$names')"
fi
PASS "subsample -n > input preserves order"
