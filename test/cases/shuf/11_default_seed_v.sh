#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# With no -s flag and -v, yamshuf should print the chosen (time-based)
# RNG seed to stderr so the run is reproducible.
assert_stderr_contains "shuf -v default prints 'RNG seed:'" \
    "RNG seed:" \
    "$YAMTK" shuf -v -i "$TESTDIR/dna.fa"
# Explicit -s 4 should report the user value
assert_stderr_contains "shuf -v -s 4 prints 'RNG seed: 4'" \
    "RNG seed: 4" \
    "$YAMTK" shuf -v -s 4 -i "$TESTDIR/dna.fa"
