#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

# Range of length 4 with k=3 (2k=6): range should be skipped with a warning.
bed=$(mktemp /tmp/yamtk_test_XXXXXX)
printf "1\t0\t4\tshort\t.\t+\n" > "$bed"

assert_exit0 "shuf -x too-short range exits 0" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -x "$bed" -k 3 -s 42

assert_stderr_contains "shuf -x too-short range emits warning" "Warning" \
    "$YAMTK" shuf -v -i "$TESTDIR/dna.fa" -x "$bed" -k 3 -s 42

rm -f "$bed"
