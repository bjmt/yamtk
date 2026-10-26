#!/usr/bin/env bash
# Audit: negative counts wrapped through strtoull (seq dup -n -1 ran forever),
# and -b nan / -b 1e-300 slipped past the hist range check.
source "$TESTDIR/lib.sh"
assert_stderr_contains "seq dup -n -1" "parse" "$YAMTK" seq -a dup -i "$TESTDIR/dna.fa" -n -1
assert_stderr_contains "seq subsample -n -1" "parse" "$YAMTK" seq -a subsample -i "$TESTDIR/dna.fa" -n -1
assert_exit_nonzero "seq hist -b nan" "$YAMTK" seq -a hist -i "$TESTDIR/dna.fa" -b nan
assert_exit_nonzero "seq hist -b 1e-300" "$YAMTK" seq -a hist -i "$TESTDIR/dna.fa" -b 1e-300
assert_exit0 "seq hist -b 0.5" "$YAMTK" seq -a hist -i "$TESTDIR/dna.fa" -b 0.5
assert_exit_nonzero "shuf -r -1" "$YAMTK" shuf -i "$TESTDIR/dna.fa" -k 2 -r -1
assert_exit_nonzero "seed -n -1" "$YAMTK" seed -m "$TESTDIR/motif.meme" -i "$TESTDIR/dna.fa" -n -1
