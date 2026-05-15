#!/usr/bin/env bash
# Test 04: -r 0 returns a motif (no refinement, no extend) and exits 0
source "${TESTDIR}/lib.sh"

assert_exit0 "yamtk ref -r 0 with no extend" \
  "$YAMTK" ref -r 0 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

# -r 0 with no extend and no trim: output preserves the original 7bp width and consensus.
assert_stdout_contains "preserves original width" "w= 7" \
  "$YAMTK" ref -r 0 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_stdout_contains "preserves original consensus" "MOTIF motif_1 ACGTACG" \
  "$YAMTK" ref -r 0 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
