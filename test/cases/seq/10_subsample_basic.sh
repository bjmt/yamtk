#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Subsample 2 of 3 records.
out=$("$YAMTK" seq -a subsample -n 2 -s 42 -i "$TESTDIR/dna.fa" 2>/dev/null)
n=$(echo "$out" | grep -c '^>')
if [ "$n" -ne 2 ]; then
    FAIL "subsample -n 2 emits 2 records (got $n)"
fi
PASS "subsample -n 2 emits exactly 2 records"

# Output must be in input order: names are integers, must be strictly increasing.
ordered=$(echo "$out" | awk '/^>/ {print $1}' | tr -d '>' \
    | awk 'NR==1 || $1>prev {prev=$1; ok=1; next} {ok=0; exit} END {print ok}')
if [ "$ordered" != "1" ]; then
    names=$(echo "$out" | awk '/^>/ {printf "%s ", $1}')
    FAIL "subsample preserves input order (got '$names')"
fi
PASS "subsample output preserves input order"
