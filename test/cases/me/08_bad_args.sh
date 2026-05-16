#!/usr/bin/env bash
# Test 08: invalid arguments exit non-zero with informative stderr
source "${TESTDIR}/lib.sh"

assert_exit_nonzero "no -i flag" \
  "$YAMTK" me -k 8 -K 8 -N 1 -s 42 -o - -O ''

assert_stderr_contains "no -i stderr" "Missing -i" \
  "$YAMTK" me -k 8 -K 8 -N 1 -s 42 -o - -O ''

assert_exit_nonzero "-k greater than -K" \
  "$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" -k 10 -K 5 -o - -O ''

assert_stderr_contains "-k > -K stderr" "must be <= -K" \
  "$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" -k 10 -K 5 -o - -O ''
