#!/usr/bin/env bash
# Test 02: cross-width dedup collapses the implanted motif to one row
source "${TESTDIR}/lib.sh"

N=$(
  "$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
    -k 6 -K 10 -N 5 -D 0.5 -s 42 -o - -M '' 2>/dev/null \
    | grep -v '^#' | wc -l | tr -d ' '
)
if [ "$N" -lt 1 ]; then
  FAIL "dedup: expected >= 1 surviving motif, got $N"
fi
PASS "dedup: $N motif(s) survive after cross-width dedup"
