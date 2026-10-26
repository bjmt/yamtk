#!/usr/bin/env bash
# Audit: an over-long FASTA header made enr free an unset pointer on its error
# path, and streaming -l -x never checked BED ends against the sequence.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
assert_stderr_contains "enr long header (-i)" "too large" "$YAMTK" enr -i "$A/long_header.fa" -m "$TESTDIR/motif.meme" -s 1
assert_stderr_contains "enr long header (-n)" "too large" "$YAMTK" enr -i "$TESTDIR/dna.fa" -n "$A/long_header.fa" -m "$TESTDIR/motif.meme"
assert_stderr_contains "enr -x past end" "exceeds seq length" "$YAMTK" enr -i "$TESTDIR/dna.fa" -n "$TESTDIR/dna.fa" -m "$TESTDIR/motif.meme" -x "$A/bed_past_end.bed"
assert_stderr_contains "enr -l -x past end" "exceeds seq length" "$YAMTK" enr -i "$TESTDIR/dna.fa" -n "$TESTDIR/dna.fa" -m "$TESTDIR/motif.meme" -x "$A/bed_past_end.bed" -l
# -r keeps the full header; BED names are matched on the first word.
assert_exit0 "enr -r -x with description header" "$YAMTK" enr -r -i "$A/desc_header.fa" -n "$A/desc_header.fa" -m "$TESTDIR/motif.meme" -q 1 -x "$A/desc_header.bed"
