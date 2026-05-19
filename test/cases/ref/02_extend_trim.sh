#!/usr/bin/env bash
# Test 02: -e 1 -T 0.5 extends a narrowed seed and trims back to the full 8bp implant
source "${TESTDIR}/lib.sh"

assert_ref_golden "extend + IC trim recovers full implant" \
  "${TESTDIR}/expected/ref_extend_trim.txt" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o -

assert_stdout_contains "consensus = full 8bp implant" "ACGTACGT" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o -
