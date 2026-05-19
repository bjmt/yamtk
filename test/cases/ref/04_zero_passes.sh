#!/usr/bin/env bash
# Test 04: -n 0 returns a motif (no refinement, no extend) and exits 0
source "${TESTDIR}/lib.sh"

assert_exit0 "yamtk ref -n 0 with no extend" \
  "$YAMTK" ref -n 0 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

# -n 0 with no extend and no trim: output preserves the original 7bp width and consensus.
assert_stdout_contains "preserves original width" "w= 7" \
  "$YAMTK" ref -n 0 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

assert_stdout_contains "preserves original consensus" "MOTIF motif_1 ACGTACG" \
  "$YAMTK" ref -n 0 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
