#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# format rewraps lines to the default 60-bp width and trims names to first word.
assert_golden "seq format rewrap matches golden" \
    "$TESTDIR/expected/seq_format_rewrap.fa" \
    "$YAMTK" seq -a format -i "$TESTDIR/dna.fa"
# -r retains the full header.
assert_golden "seq format -r keeps full names" \
    "$TESTDIR/expected/seq_format_keepname.fa" \
    "$YAMTK" seq -a format -r -i "$TESTDIR/dna.fa"
