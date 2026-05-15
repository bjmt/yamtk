#!/usr/bin/env bash
# Test 05: JASPAR-format seed input
source "${TESTDIR}/lib.sh"

assert_ref_golden "JASPAR seed refinement" \
  "${TESTDIR}/expected/ref_jaspar.txt" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.jaspar" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
