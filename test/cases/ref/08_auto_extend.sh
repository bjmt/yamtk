#!/usr/bin/env bash
# Test 08: -E auto-extend grows flanks until IC drops, then reverts last step
source "${TESTDIR}/lib.sh"

assert_ref_golden "auto-extend recovers full implant from narrow seed" \
  "${TESTDIR}/expected/ref_auto_extend.txt" \
  "$YAMTK" ref -E -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_stdout_contains "final consensus is full 8bp implant" "MOTIF narrow_seed CGTACG ACGTACGT" \
  "$YAMTK" ref -E -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_exit_nonzero "-E with -e is rejected" \
  "$YAMTK" ref -E -e 1 -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_stderr_contains "-E + -e error message" "both -e and -E" \
  "$YAMTK" ref -E -e 1 -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa"
