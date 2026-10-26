#!/usr/bin/env bash
# Audit: a single numeric token longer than 255 chars used to overflow a
# fixed 256-byte stack buffer in every motif parser and in the -b parser.
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
for f in long_token.meme long_bkg.meme long_token.jaspar long_token.homer; do
    assert_stderr_contains "scan $f" "too long" "$YAMTK" scan -m "$A/$f" -s "$TESTDIR/dna.fa"
    assert_stderr_contains "enr $f"  "too long" "$YAMTK" enr -i "$TESTDIR/dna.fa" -m "$A/$f" -s 1
    assert_stderr_contains "ref $f"  "too long" "$YAMTK" ref -m "$A/$f" -i "$TESTDIR/dna.fa"
    assert_stderr_contains "cmp $f"  "too long" "$YAMTK" cmp -m "$A/$f" -t "$TESTDIR/motif.meme"
done
for f in long_token.meme long_token.jaspar long_token.homer; do
    assert_stderr_contains "seed $f" "too long" "$YAMTK" seed -m "$A/$f" -i "$TESTDIR/dna.fa" -f 0.5 -s 1
    assert_stderr_contains "conv $f" "too long" "$YAMTK" conv -m "$A/$f" -t jaspar
done
long=$(head -c 400 /dev/zero | tr '\0' 0).25
assert_stderr_contains "scan -b long" "too long" "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -b "$long,0.25,0.25,0.25"
assert_stderr_contains "enr -b long"  "too long" "$YAMTK" enr -i "$TESTDIR/dna.fa" -m "$TESTDIR/motif.meme" -s 1 -b "$long,0.25,0.25,0.25"
assert_stderr_contains "me -b long"   "too long" "$YAMTK" me -i "$TESTDIR/dna.fa" -o - -O - -s 1 -N 1 -k 3 -K 3 -b "$long,0.25,0.25,0.25"
assert_stderr_contains "ref -b long"  "too long" "$YAMTK" ref -m "$TESTDIR/motif.meme" -i "$TESTDIR/dna.fa" -b "$long,0.25,0.25,0.25"
assert_stderr_contains "cmp -b long"  "too long" "$YAMTK" cmp -m "$TESTDIR/motif.meme" -t "$TESTDIR/motif.meme" -b "$long,0.25,0.25,0.25"
assert_stderr_contains "scan -1 300-mer" "too large" "$YAMTK" scan -1 "$(head -c 300 /dev/zero | tr '\0' A)" -s "$TESTDIR/dna.fa"
