#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmpdir=$(mktemp -d /tmp/yamseed_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

# -X with a non-ambiguous consensus → deterministic seeding at the given range.
"$YAMTK" seed -1 CACGTG -X seq1:10-16 \
  -i "$TESTDIR/fixtures/me_implanted.fa" > "$tmpdir/s.fa" 2>/dev/null
bases=$(awk '/^>seq1$/{p=1;next} p && /^>/{p=0} p' "$tmpdir/s.fa" | tr -d '\n' | awk '{print substr($0,11,6)}')
if [ "$bases" != "CACGTG" ]; then
    FAIL "seed -X fwd: expected 'CACGTG' at seq1[10..16), got '$bases'"
fi
PASS "seed -X forward strand places consensus exactly"

# -X with explicit '-' strand → RC sampled bases.
"$YAMTK" seed -1 GGGCGG -X seq1:50-56:- \
  -i "$TESTDIR/fixtures/me_implanted.fa" > "$tmpdir/s2.fa" 2>/dev/null
bases=$(awk '/^>seq1$/{p=1;next} p && /^>/{p=0} p' "$tmpdir/s2.fa" | tr -d '\n' | awk '{print substr($0,51,6)}')
if [ "$bases" != "CCGCCC" ]; then
    FAIL "seed -X -strand: expected 'CCGCCC' (RC of GGGCGG), got '$bases'"
fi
PASS "seed -X '-' strand places RC of consensus"

# -X is mutually exclusive with -f and with -x
assert_exit_nonzero "seed -X + -f fails" \
    "$YAMTK" seed -1 CACGTG -X seq1:10-16 -f 0.01 \
                  -i "$TESTDIR/fixtures/me_implanted.fa"

# -X with multi-motif -m must error
assert_exit_nonzero "seed -X with multi-motif -m fails" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" -X seq1:10-16 \
                  -i "$TESTDIR/fixtures/me_implanted.fa"

# Malformed -X
assert_exit_nonzero "seed -X malformed (missing dash)" \
    "$YAMTK" seed -1 CACGTG -X "seq1:1016" \
                  -i "$TESTDIR/fixtures/me_implanted.fa"
assert_exit_nonzero "seed -X start >= end" \
    "$YAMTK" seed -1 CACGTG -X "seq1:16-10" \
                  -i "$TESTDIR/fixtures/me_implanted.fa"
