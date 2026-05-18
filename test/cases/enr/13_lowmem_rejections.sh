#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -l + -T ranksum -> error
assert_exit_nonzero "enr -l + -T ranksum is rejected" \
    "$YAMTK" enr -l -T ranksum -i "$TESTDIR/fixtures/enr_pos.fa" \
                                -n "$TESTDIR/fixtures/enr_neg.fa" \
                                -m "$TESTDIR/fixtures/enr_motifs.meme"
assert_stderr_contains "enr -l + ranksum error mentions ranksum" \
    "ranksum" \
    "$YAMTK" enr -l -T ranksum -i "$TESTDIR/fixtures/enr_pos.fa" \
                                -n "$TESTDIR/fixtures/enr_neg.fa" \
                                -m "$TESTDIR/fixtures/enr_motifs.meme"
# -l + stdin -i -  -> error
assert_exit_nonzero "enr -l + stdin (-i -) is rejected" \
    "$YAMTK" enr -l -i - -n "$TESTDIR/fixtures/enr_neg.fa" \
                         -m "$TESTDIR/fixtures/enr_motifs.meme"
# -l + -j 4 should NOT error (auto-downgrade); just succeed
assert_exit0 "enr -l + -j 4 succeeds (auto-downgrade)" \
    "$YAMTK" enr -l -j 4 -i "$TESTDIR/fixtures/enr_pos.fa" \
                          -n "$TESTDIR/fixtures/enr_neg.fa" \
                          -m "$TESTDIR/fixtures/enr_motifs.meme"
