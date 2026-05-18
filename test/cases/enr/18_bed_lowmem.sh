#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -l + -x with -n should yield the same counts as the fully-loaded path
# (modulo stat noise that's a function of identical counts: identical).
nonl=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -x "$TESTDIR/fixtures/enr_pos.bed" 2>/dev/null | grep -v '^#')
lm=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -x "$TESTDIR/fixtures/enr_pos.bed" -l 2>/dev/null | grep -v '^#')
if [ "$nonl" != "$lm" ]; then
    echo "Non-l output:" >&2; echo "$nonl" >&2
    echo "-l output:"   >&2; echo "$lm"   >&2
    FAIL "low-mem + -x matches fully-loaded + -x"
fi
PASS "low-mem + -x matches fully-loaded + -x"
