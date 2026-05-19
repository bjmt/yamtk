#!/usr/bin/env bash
# Test 10: per-motif stderr summary always prints; -Q drops motifs when refined IC < seed IC
source "${TESTDIR}/lib.sh"

# Summary line on stderr (no -v needed)
assert_stderr_contains "summary reports width transition" "w: 6->8" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o /dev/null

assert_stderr_contains "summary reports IC delta" "IC: " \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o /dev/null

assert_stderr_contains "summary reports hit count" "hits: " \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o /dev/null

# Refinement that improves IC passes the gate
assert_exit0 "-Q passes when IC goes up" \
  "$YAMTK" ref -Q -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 1 -T 0.5 -o /dev/null

# Refinement that doesn't improve IC is dropped under -Q (the already-fit seed
# has nsites=38 baked in; rebuilding from positives lands at a slightly lower IC).
assert_exit_nonzero "-Q drops when IC goes down" \
  "$YAMTK" ref -Q -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o /dev/null

assert_stderr_contains "-Q reports the drop" "dropping under -Q" \
  "$YAMTK" ref -Q -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -o /dev/null
