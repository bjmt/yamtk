#!/usr/bin/env bash
# Test 03: two runs with the same seed produce identical output
source "${TESTDIR}/lib.sh"

RUN1=$("$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 8 -K 8 -N 2 -s 42 -o - -M '' 2>/dev/null | grep -v '^##yamme ')
RUN2=$("$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 8 -K 8 -N 2 -s 42 -o - -M '' 2>/dev/null | grep -v '^##yamme ')

if [ "$RUN1" != "$RUN2" ]; then
  FAIL "reproducibility: two runs with -s 42 differ"
fi
PASS "reproducibility: two runs with -s 42 are identical"
