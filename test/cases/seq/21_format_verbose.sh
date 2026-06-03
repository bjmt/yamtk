#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Verbose mode prints a file-wide composition report to stderr.
assert_stderr_contains "seq format -v reports composition" \
    "Composition" \
    "$YAMTK" seq -a format -v -i "$TESTDIR/dna.fa"
assert_stderr_contains "seq format -v reports soft-mask count" \
    "soft-masked" \
    "$YAMTK" seq -a format -v -i "$TESTDIR/dna.fa"
# Very verbose adds a per-sequence line (named with the first sequence's name).
assert_stderr_contains "seq format -w reports per-sequence line" \
    "\[1\]" \
    "$YAMTK" seq -a format -w -i "$TESTDIR/dna.fa"
