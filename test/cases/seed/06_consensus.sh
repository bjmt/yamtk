#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmpdir=$(mktemp -d /tmp/yamseed_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

# Deterministic placement: with a non-ambiguous IUPAC consensus, the seeded
# bases must equal the consensus exactly.
cat > "$tmpdir/site.bed" <<'EOF'
seq1	10	16	CACGTG	0	+
EOF
"$YAMTK" seed -1 CACGTG \
  -i "$TESTDIR/fixtures/me_implanted.fa" \
  -x "$tmpdir/site.bed" > "$tmpdir/seeded.fa" 2>/dev/null
bases=$(awk '/^>seq1$/{p=1;next} p && /^>/{p=0} p' "$tmpdir/seeded.fa" | tr -d '\n' | awk '{print substr($0,11,6)}')
if [ "$bases" != "CACGTG" ]; then
    FAIL "seed -1 CACGTG did not deterministically seed CACGTG; got '$bases'"
fi
PASS "seed -1 with non-ambiguous consensus deterministically seeds the consensus"

# Mutex: -m and -1 cannot be combined
assert_exit_nonzero "seed -m + -1 fails" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" -1 CACGTG \
                  -i "$TESTDIR/fixtures/me_implanted.fa" -f 0.01

# Invalid IUPAC character rejected
assert_exit_nonzero "seed -1 with invalid IUPAC char fails" \
    "$YAMTK" seed -1 CAxxTG -i "$TESTDIR/fixtures/me_implanted.fa" -f 0.01
