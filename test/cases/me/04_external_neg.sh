#!/usr/bin/env bash
# Test 04: external negative FASTA, golden match
source "${TESTDIR}/lib.sh"

assert_me_golden "external neg k=6 s=42" \
  "${TESTDIR}/expected/me_external_neg.txt" \
  "$YAMTK" me \
  -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -n "${TESTDIR}/fixtures/me_external_neg.fa" \
  -N 1 -k 6 -K 6 -s 42 -o - -M ''
