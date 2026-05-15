#!/usr/bin/env bash
# Test 09: -j 1 and -j 4 produce byte-identical TSV output at the same seed
source "${TESTDIR}/lib.sh"

T1=$(mktemp -t yamme_j1.XXXXXX)
T4=$(mktemp -t yamme_j4.XXXXXX)

"$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 6 -K 10 -N 3 -s 42 -j 1 -o "$T1" -M '' 2>/dev/null
"$YAMTK" me -i "${TESTDIR}/fixtures/me_implanted.fa" \
  -k 6 -K 10 -N 3 -s 42 -j 4 -o "$T4" -M '' 2>/dev/null

if ! diff -u <(grep -v '^##yamme ' "$T1") <(grep -v '^##yamme ' "$T4") >/dev/null 2>&1; then
  diff -u <(grep -v '^##yamme ' "$T1") <(grep -v '^##yamme ' "$T4") >&2 || true
  rm -f "$T1" "$T4"
  FAIL "-j 1 and -j 4 with -s 42 produce different output"
fi

rm -f "$T1" "$T4"
PASS "-j 1 and -j 4 with -s 42 are byte-identical"
