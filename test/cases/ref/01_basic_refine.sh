#!/usr/bin/env bash
# Test 01: refine a MEME seed against implanted positives
source "${TESTDIR}/lib.sh"

assert_ref_golden "basic MEME refine" \
  "${TESTDIR}/expected/ref_basic.txt" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
