#!/usr/bin/env bash
# Audit: nan/inf tokens and negative -p/-N or -b values reached log2() and a
# NaN-to-int cast; all of them must now be rejected up front.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
assert_exit_nonzero "scan HOMER nan" "$YAMTK" scan -m "$A/nan_values.homer" -s "$TESTDIR/dna.fa"
assert_stderr_contains "scan -p -1" "positive" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -p -1
assert_stderr_contains "scan -N -1" "positive" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -N -1
assert_stderr_contains "enr -p -1000" "positive" "$YAMTK" enr -i "$TESTDIR/dna.fa" -m "$TESTDIR/motif.meme" -s 1 -p -1000 -N 1000
assert_stderr_contains "scan -b nan" "background" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -b nan,0.25,0.25,0.25
assert_stderr_contains "ref -b negative" "non-negative" "$YAMTK" ref -1 CTCGA -i "$TESTDIR/dna.fa" -b -0.5,0.5,0.5,0.5
assert_stderr_contains "me -b nan" "background" "$YAMTK" me -i "$TESTDIR/dna.fa" -o - -O - -s 1 -N 1 -k 3 -K 3 -b nan,0.25,0.25,0.25
assert_exit_nonzero "scan -t nan" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t nan
