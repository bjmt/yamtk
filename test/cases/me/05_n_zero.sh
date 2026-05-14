#!/usr/bin/env bash
# Test 05: -N 0 exits 0 and produces no data rows
source "${TESTDIR}/lib.sh"

assert_exit0 "-N 0 exits 0" \
  "$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" -N 0 -o - -M ''

N=$("$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" -N 0 -o - -M '' 2>/dev/null \
  | { grep -v '^#' || true; } | wc -l | tr -d ' ')
if [ "$N" -ne 0 ]; then
  FAIL "-N 0: expected 0 data rows, got $N"
fi
PASS "-N 0: 0 data rows"
