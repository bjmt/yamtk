#!/usr/bin/env bash
# Test 06: refined motif can be fed back through yamtk scan
source "${TESTDIR}/lib.sh"

tmp="$(mktemp -d -t yamref-roundtrip-XXXXXX)"
trap 'rm -rf "$tmp"' EXIT

assert_exit0 "ref produces MEME output" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o "$tmp/refined.meme"

assert_exit0 "scan accepts refined MEME" \
  "$YAMTK" scan -m "$tmp/refined.meme" -s "${TESTDIR}/fixtures/me_implanted.fa" -j 1

assert_stdout_contains "scan finds hits in implants" "ACGTACG" \
  "$YAMTK" scan -m "$tmp/refined.meme" -s "${TESTDIR}/fixtures/me_implanted.fa" -j 1
