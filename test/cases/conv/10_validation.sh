#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_exit_nonzero "conv rejects unknown target format" \
    "$YAMTK" conv -m "$TESTDIR/fixtures/enr_motifs.meme" -t bogus
assert_exit_nonzero "conv requires -t" \
    "$YAMTK" conv -m "$TESTDIR/fixtures/enr_motifs.meme"
assert_exit_nonzero "conv requires -m" \
    "$YAMTK" conv -t meme
