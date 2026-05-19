#!/usr/bin/env bash
# Test 13: -r preserves the identifier+altname pair on output (matches sibling
# tools' "do not trim names" convention). Without -r, the altname is dropped
# and the refined consensus serves as the sole altname.
source "${TESTDIR}/lib.sh"

# Default: drop the seed altname; refined consensus is the altname.
# ref_seed.meme is "MOTIF motif_1 ACGTACG" — under refinement the consensus
# stays ACGTACG, so the output line is also "MOTIF motif_1 ACGTACG"
# (id + consensus, not id + seed-altname).
assert_stdout_contains "default drops seed altname (id + consensus)" \
  "^MOTIF motif_1 ACGTACG$" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
    -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

# With -r: preserve the full identifier + altname; no consensus appended.
# Output should be exactly "MOTIF motif_1 ACGTACG" (the original line).
# Hard to distinguish from above for this fixture because consensus == altname,
# so use ref_seed_narrow.meme where the seed altname (CGTACG) differs from the
# refined consensus (~ACGTACGT after auto-extend).
assert_stdout_contains "-r preserves the seed's identifier+altname pair" \
  "^MOTIF narrow_seed CGTACG$" \
  "$YAMTK" ref -r -E -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
    -i "${TESTDIR}/fixtures/me_implanted.fa" -o -

# And confirm the default for the same input does NOT keep the seed's altname:
# auto-extend recovers ACGTACGT, so the output line is "MOTIF narrow_seed ACGTACGT".
assert_stdout_contains "default replaces seed altname with refined consensus" \
  "^MOTIF narrow_seed ACGTACGT$" \
  "$YAMTK" ref -E -m "${TESTDIR}/fixtures/ref_seed_narrow.meme" \
    -i "${TESTDIR}/fixtures/me_implanted.fa" -o -
