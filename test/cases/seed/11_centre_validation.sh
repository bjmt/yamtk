#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Rejections.
assert_exit_nonzero "seed -c 0 fails" \
    "$YAMTK" seed -1 CACGTG -i "$TESTDIR/fixtures/me_implanted.fa" -f 0.01 -c 0
assert_exit_nonzero "seed -c 100 fails (cap)" \
    "$YAMTK" seed -1 CACGTG -i "$TESTDIR/fixtures/me_implanted.fa" -f 0.01 -c 100
assert_exit_nonzero "seed -c bad value fails" \
    "$YAMTK" seed -1 CACGTG -i "$TESTDIR/fixtures/me_implanted.fa" -f 0.01 -c not_a_number
assert_exit_nonzero "seed -c > 1 without -f fails (single-range)" \
    "$YAMTK" seed -1 CACGTG -i "$TESTDIR/fixtures/me_implanted.fa" -X "1:0-50" -c 2
