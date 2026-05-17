#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Input has 3 seqs; with -n 4 we expect 12 output records.
n_in=$("$YAMTK" seq -a stats -i "$TESTDIR/dna.fa" 2>/dev/null | grep -v '^#' | wc -l | tr -d ' ')
n_out=$("$YAMTK" seq -a dup -n 4 -i "$TESTDIR/dna.fa" 2>/dev/null | grep -c '^>')
expected=$((n_in * 4))
if [ "$n_out" -ne "$expected" ]; then
    FAIL "dup output record count: expected $expected, got $n_out"
fi
PASS "dup -n 4 produces N*input records"

# Suffix check: first input ("1") should yield 1-1 ... 1-4
suffixes=$("$YAMTK" seq -a dup -n 4 -i "$TESTDIR/dna.fa" 2>/dev/null \
  | awk '/^>1-/{print}' | sort -u | tr '\n' ' ')
expected_suffix=">1-1 >1-2 >1-3 >1-4 "
if [ "$suffixes" != "$expected_suffix" ]; then
    FAIL "dup name suffix mismatch: got '$suffixes' expected '$expected_suffix'"
fi
PASS "dup name suffixes are name-1 .. name-N"
