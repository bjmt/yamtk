#!/usr/bin/env bash
# Audit: ref and cmp adopted the user's LC_NUMERIC, so a comma-decimal locale
# broke parsing of their own motif files.
source "$TESTDIR/lib.sh"
# (-t 0.05 so the seed has hits on the small fixture; the failure mode being
# tested is "Failed to parse background value: 0.249493".)
assert_exit0 "ref under de_DE" env LC_ALL=de_DE.UTF-8 "$YAMTK" ref -m "$TESTDIR/motif.meme" -i "$TESTDIR/dna.fa" -t 0.05 -n 0
assert_exit0 "cmp under de_DE" env LC_ALL=de_DE.UTF-8 "$YAMTK" cmp -m "$TESTDIR/motif.meme" -t "$TESTDIR/motif.meme" -q 1
