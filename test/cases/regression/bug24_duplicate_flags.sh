#!/usr/bin/env bash
# Bug #24: yamshuf and yamdedup would silently reopen -i/-o on duplicate flags,
# leaking the first gzFile/FILE*. yamscan had the same issue for -s/-o.
# Fix: error out with "specified more than once" if the flag is already set.
source "$TESTDIR/lib.sh"

assert_exit_nonzero "shuf -i twice exits non-zero" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -i "$TESTDIR/dna.fa" -s 42

assert_stderr_contains "shuf -i twice prints 'more than once'" \
    "more than once" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -i "$TESTDIR/dna.fa" -s 42

tmp=$(mktemp /tmp/yamtk_XXXXXX)
trap 'rm -f "$tmp"' EXIT

assert_exit_nonzero "shuf -o twice exits non-zero" \
    "$YAMTK" shuf -i "$TESTDIR/dna.fa" -s 42 -o "$tmp" -o "$tmp"

assert_exit_nonzero "dedup -i twice exits non-zero" \
    "$YAMTK" dedup -i "$TESTDIR/fixtures/scan_output_for_dedup.txt" \
                   -i "$TESTDIR/fixtures/scan_output_for_dedup.txt"

assert_exit_nonzero "scan -s twice exits non-zero" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" \
                  -s "$TESTDIR/dna.fa" -s "$TESTDIR/dna.fa" -t 0.05 -j 1

assert_exit_nonzero "scan -o twice exits non-zero" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" \
                  -t 0.05 -j 1 -o "$tmp" -o "$tmp"
