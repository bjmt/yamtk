#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf k=2 (Eulerian) golden output" \
    "$TESTDIR/expected/shuf_k2.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -k 2 -s 42
