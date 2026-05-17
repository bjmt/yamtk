#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_golden "seed random (-f 0.01 -s 1) FASTA matches golden" \
    "$TESTDIR/expected/seed_random_s1.fa" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
                  -i "$TESTDIR/fixtures/me_implanted.fa" \
                  -f 0.01 -s 1
