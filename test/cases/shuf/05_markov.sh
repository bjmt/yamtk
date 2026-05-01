#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "shuf Markov (-m)" \
    "$TESTDIR/expected/shuf_markov.txt" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -m -k 2 -s 42
