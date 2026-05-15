#!/usr/bin/env bash
# Test 09: -1 builds a seed motif from an IUPAC consensus string
source "${TESTDIR}/lib.sh"

assert_ref_golden "exact consensus seed refines to implant" \
  "${TESTDIR}/expected/ref_consensus.txt" \
  "$YAMTK" ref -1 ACGTACGT -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

# Narrow consensus + auto-extend should recover the full implant
assert_stdout_contains "narrow consensus + -E recovers ACGTACGT" "MOTIF CGTACG ACGTACGT" \
  "$YAMTK" ref -1 CGTACG -E -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_exit_nonzero "-m and -1 conflict" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" -1 ACGT \
  -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_stderr_contains "-m + -1 error message" "cannot both be used" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" -1 ACGT \
  -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_exit_nonzero "unknown letter in consensus" \
  "$YAMTK" ref -1 ACGTQ -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_stderr_contains "unknown letter error message" "Unknown letter" \
  "$YAMTK" ref -1 ACGTQ -i "${TESTDIR}/fixtures/me_implanted.fa"
