#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# subsample needs either -n or -f
assert_exit_nonzero "subsample without -n/-f fails" \
    "$YAMTK" seq -a subsample -i "$TESTDIR/dna.fa"
# -n and -f are mutually exclusive
assert_exit_nonzero "subsample with both -n and -f fails" \
    "$YAMTK" seq -a subsample -n 2 -f 0.5 -i "$TESTDIR/dna.fa"
# -f must be in (0,1)
assert_exit_nonzero "subsample -f 0 fails" \
    "$YAMTK" seq -a subsample -f 0 -i "$TESTDIR/dna.fa"
assert_exit_nonzero "subsample -f 1 fails" \
    "$YAMTK" seq -a subsample -f 1 -i "$TESTDIR/dna.fa"
assert_exit_nonzero "subsample -f >1 fails" \
    "$YAMTK" seq -a subsample -f 1.5 -i "$TESTDIR/dna.fa"
# -f and -s outside subsample fail
assert_exit_nonzero "-f outside subsample fails" \
    "$YAMTK" seq -a stats -f 0.5 -i "$TESTDIR/dna.fa"
assert_exit_nonzero "-s outside subsample fails" \
    "$YAMTK" seq -a stats -s 7 -i "$TESTDIR/dna.fa"
