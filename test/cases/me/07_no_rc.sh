#!/usr/bin/env bash
# Test 07: -R (no RC) is deterministic and differs from RC-enabled run
source "${TESTDIR}/lib.sh"

RUN1=$("$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 8 -K 8 -N 1 -R -s 42 -o - -M '' 2>/dev/null | grep -v '^##yamme ')
RUN2=$("$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 8 -K 8 -N 1 -R -s 42 -o - -M '' 2>/dev/null | grep -v '^##yamme ')

if [ "$RUN1" != "$RUN2" ]; then
  FAIL "no-RC reproducibility: two -R runs with -s 42 differ"
fi
PASS "no-RC: two -R runs with -s 42 are identical"
