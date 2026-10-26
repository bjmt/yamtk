#!/usr/bin/env bash
# Audit: reverse-complementing a zero-length window dereferenced a NULL buffer.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
for s in 1 2 3 4 5 6; do
    assert_exit0 "bkg empty target -r -N 1 -s $s" "$YAMTK" bkg -i "$A/header_only.fa" -g "$A/genome_short.fa" -r -N 1 -s $s
done
