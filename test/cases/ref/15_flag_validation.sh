#!/usr/bin/env bash
# Test 15: legacy / malformed flags fail cleanly under the new flag set.
source "${TESTDIR}/lib.sh"

# -T now requires a numeric arg; bare -T should fail (or eat the next arg).
# Using -T at end of args triggers the "option requires an argument" path.
assert_exit_nonzero "-T without argument fails" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
    -i "${TESTDIR}/fixtures/me_implanted.fa" -T

# -I has been removed; unknown option.
assert_exit_nonzero "-I no longer accepted" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
    -i "${TESTDIR}/fixtures/me_implanted.fa" -E -I 0.5

# Old-style -r N (integer arg) — N is mistaken for a positional arg and yamref
# rejects extras. In practice this manifests as an unknown / extra argument.
# Use -r 2 -e 0 -m ... ; if -r is now a no-arg flag, "2" becomes a stray arg.
# getopt + main_ref together should not produce successful refinement output
# of the form we'd expect.
assert_exit_nonzero "old-style -r 2 (numeric arg) fails" \
  "$YAMTK" ref -r 2 -m "${TESTDIR}/fixtures/ref_seed.meme" \
    -i "${TESTDIR}/fixtures/me_implanted.fa"
