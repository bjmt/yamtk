#!/usr/bin/env bash
# Audit: the dedup hash key concatenated strand+motif+seq with no separator, so
# ("M1","1chr") and ("M11","chr") collided and were deduplicated together.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
assert_line_count_eq "dedup keeps both" 2 "$YAMTK" dedup -i "$A/dedup_collide.bed" -b
