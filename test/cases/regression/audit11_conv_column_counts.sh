#!/usr/bin/env bash
# Audit: PCM round trips rescaled every column to column 0's depth.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
assert_stdout_contains "jaspar keeps column 2 depth" "150" "$YAMTK" conv -m "$A/unequal_cols.jaspar" -t jaspar
assert_stdout_contains "hocomoco keeps column 2 depth" "0 150 0 50" "$YAMTK" conv -m "$A/unequal_cols.jaspar" -t hocomoco
