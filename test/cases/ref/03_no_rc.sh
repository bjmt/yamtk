#!/usr/bin/env bash
# Test 03: -R disables reverse-strand scoring
source "${TESTDIR}/lib.sh"

assert_ref_golden "forward-only refinement" \
  "${TESTDIR}/expected/ref_norc.txt" \
  "$YAMTK" ref -R -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_stdout_contains "strands header is forward only" "strands: +$" \
  "$YAMTK" ref -R -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
