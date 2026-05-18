#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# yamseed prints the chosen RNG seed under -v.
assert_stderr_contains "seed -v default prints 'RNG seed:'" \
    "RNG seed:" \
    "$YAMTK" seed -v -1 CACGTG -X seq1:10-16 \
                  -i "$TESTDIR/fixtures/me_implanted.fa"
# Explicit -s reports the user value.
assert_stderr_contains "seed -v -s 7 prints 'RNG seed: 7'" \
    "RNG seed: 7" \
    "$YAMTK" seed -v -s 7 -1 CACGTG -X seq1:10-16 \
                  -i "$TESTDIR/fixtures/me_implanted.fa"
