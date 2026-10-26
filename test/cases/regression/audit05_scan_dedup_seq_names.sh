#!/usr/bin/env bash
# Audit: scan -d appended the "__N<i>" suffix past the exact-length name buffer.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
assert_exit0 "scan -d duplicate names" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$A/dup_names.fa" -t 0.9 -d
assert_stdout_contains "renamed duplicate" "d__N2" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$A/dup_names.fa" -t 0.9 -d
assert_exit0 "scan -d -l duplicate names" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$A/dup_names.fa" -t 0.9 -d -l
