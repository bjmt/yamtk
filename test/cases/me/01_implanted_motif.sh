#!/usr/bin/env bash
# Test 01: implanted motif recovered with fixed width and seed
source "${TESTDIR}/lib.sh"

assert_me_golden "implanted ACGTACGT w=8 s=42" \
  "${TESTDIR}/expected/me_implanted_w8.txt" \
  "$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 8 -K 8 -N 1 -s 42 -o - -M ''
