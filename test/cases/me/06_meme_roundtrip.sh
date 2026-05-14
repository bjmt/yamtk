#!/usr/bin/env bash
# Test 06: MEME output from me can be scanned by yamtk scan
source "${TESTDIR}/lib.sh"

TMPDIR=$(mktemp -d)

"$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 8 -K 8 -N 1 -s 42 -o /dev/null -M "${TMPDIR}/out.meme" 2>/dev/null

HITS=$("$YAMTK" scan -m "${TMPDIR}/out.meme" \
  -s "${TESTDIR}/fixtures/me_implanted.fa" -j 1 2>/dev/null \
  | grep -v '^#' | wc -l | tr -d ' ')

rm -rf "$TMPDIR"

if [ "$HITS" -lt 1 ]; then
  FAIL "MEME round-trip: scan found 0 hits, expected >= 1"
fi
PASS "MEME round-trip: scan found $HITS hit(s)"
