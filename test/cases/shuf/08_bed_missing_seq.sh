#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

bed=$(mktemp /tmp/yamtk_test_XXXXXX)
printf "notasequence\t0\t10\tX\t.\t+\n" > "$bed"

assert_exit_nonzero "shuf -x missing seq name exits non-zero" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -x "$bed" -k 1 -s 42

assert_stderr_contains "shuf -x missing seq name error message" "not in input" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -x "$bed" -k 1 -s 42

rm -f "$bed"
