#!/usr/bin/env bash
# Test 10: random sequences with nothing planted should give no motifs.
# (me picks seeds as the best of thousands of k-mers, so without correcting
# for that search it used to report motifs here, e.g. CAGCAG at -s 12.)
source "${TESTDIR}/lib.sh"

for seed in 1 12 42; do
  n=$("$YAMTK" me -i "${TESTDIR}/fixtures/me_null.fa" -s "$seed" -o - -M '' 2>/dev/null \
      | awk '!/^#/ && NF { n++ } END { print n + 0 }')
  if [ "$n" -ne 0 ]; then
    FAIL "null positives give no motifs (seed $seed: got $n)"
  fi
done
PASS "null positives give no motifs (3 seeds)"
