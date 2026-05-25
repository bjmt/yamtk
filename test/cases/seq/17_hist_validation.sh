#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -b only valid with -a hist
assert_exit_nonzero "-b outside hist fails" \
    "$YAMTK" seq -a stats -b 0.1 -i "$TESTDIR/dna.fa"
# -b must be > 0
assert_exit_nonzero "hist -b 0 fails" \
    "$YAMTK" seq -a hist -b 0 -i "$TESTDIR/dna.fa"
# -b must be <= 1
assert_exit_nonzero "hist -b 1.5 fails" \
    "$YAMTK" seq -a hist -b 1.5 -i "$TESTDIR/dna.fa"
# unparseable -b fails
assert_exit_nonzero "hist -b abc fails" \
    "$YAMTK" seq -a hist -b abc -i "$TESTDIR/dna.fa"
# -b 1.0 is allowed (single bin)
assert_exit0 "hist -b 1.0 succeeds" \
    "$YAMTK" seq -a hist -b 1.0 -i "$TESTDIR/dna.fa"
