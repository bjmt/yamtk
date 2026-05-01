#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf k=1 (Fisher-Yates)" \
    "$TESTDIR/expected/shuf_k1.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -k 1 -s 42
