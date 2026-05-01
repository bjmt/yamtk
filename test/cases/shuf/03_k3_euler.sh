#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf k=3 (Eulerian) golden output" \
    "$TESTDIR/expected/shuf_k3.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -k 3 -s 42
