#!/usr/bin/env bash
# Test 03: -R disables reverse-strand scoring.
# Explicit -t 5e-4 because the default 1e-4 is unreachable for a 6-7bp
# motif on a small forward-only set (would fall back to p=0.01).
source "${TESTDIR}/lib.sh"

assert_ref_golden "forward-only refinement" \
  "${TESTDIR}/expected/ref_norc.txt" \
  "$YAMTK" ref -R -t 5e-4 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_stdout_contains "strands header is forward only" "strands: +$" \
  "$YAMTK" ref -R -t 5e-4 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
