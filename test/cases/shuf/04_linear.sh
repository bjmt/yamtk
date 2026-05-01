#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf linear (-l)" \
    "$TESTDIR/expected/shuf_linear.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -l -k 3 -s 42
