#!/usr/bin/env bash
# Audit: a letter-probability matrix before any MOTIF line indexed motifs[-1].
source "$TESTDIR/lib.sh"
A="$TESTDIR/fixtures/audit"
assert_stderr_contains "conv" "before any MOTIF" "$YAMTK" conv -m "$A/no_motif_line.meme" -t jaspar
assert_stderr_contains "ref"  "before any MOTIF" "$YAMTK" ref -m "$A/no_motif_line.meme" -i "$TESTDIR/dna.fa"
assert_stderr_contains "cmp"  "before any MOTIF" "$YAMTK" cmp -m "$A/no_motif_line.meme" -t "$TESTDIR/motif.meme"
assert_stderr_contains "seed" "before any MOTIF" "$YAMTK" seed -m "$A/no_motif_line.meme" -i "$TESTDIR/dna.fa" -f 0.5 -s 1
