#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

assert_golden "bkg pool basic matches golden" \
    "$TESTDIR/expected/bkg_pool_basic.fa" \
    "$YAMTK" bkg -i "$TESTDIR/fixtures/enr_pos.fa" -p "$TESTDIR/fixtures/enr_neg.fa" \
        -s 1 -G 0.1 -T 50
