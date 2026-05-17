#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seed -x basic FASTA matches golden" \
    "$TESTDIR/expected/seed_bed_basic.fa" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
                  -i "$TESTDIR/fixtures/me_implanted.fa" \
                  -x "$TESTDIR/fixtures/seed_basic.bed" -s 1
