#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -a is required
assert_exit_nonzero "seq without -a fails" \
    "$YAMTK" seq -i "$TESTDIR/dna.fa"
# -x not valid with non-BED actions
assert_exit_nonzero "seq -a rc with -x fails" \
    "$YAMTK" seq -a rc -x "$TESTDIR/dna.bed" -i "$TESTDIR/dna.fa"
# -a subset without -x fails
assert_exit_nonzero "seq -a subset without -x fails" \
    "$YAMTK" seq -a subset -i "$TESTDIR/dna.fa"
# -a dup without -n fails
assert_exit_nonzero "seq -a dup without -n fails" \
    "$YAMTK" seq -a dup -i "$TESTDIR/dna.fa"
# -N only with mask
assert_exit_nonzero "seq -N without -a mask fails" \
    "$YAMTK" seq -a rc -N -i "$TESTDIR/dna.fa"
# -l only with format
assert_exit_nonzero "seq -l without -a format fails" \
    "$YAMTK" seq -a rc -l 80 -i "$TESTDIR/dna.fa"
# unknown action fails
assert_exit_nonzero "seq -a foo (unknown action) fails" \
    "$YAMTK" seq -a foo -i "$TESTDIR/dna.fa"
