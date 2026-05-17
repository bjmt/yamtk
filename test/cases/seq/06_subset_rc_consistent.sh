#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# A '-' strand subset should equal the RC of the forward slice.
tmpdir=$(mktemp -d /tmp/yamseq_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
# Forward extraction (BED with same range but '+')
cat > "$tmpdir/fwd.bed" <<'EOF'
2	10	48	X	.	+
EOF
cat > "$tmpdir/rev.bed" <<'EOF'
2	10	48	X	.	-
EOF
fwd=$("$YAMTK" seq -a subset -x "$tmpdir/fwd.bed" -i "$TESTDIR/dna.fa" 2>/dev/null | grep -v '^>' | tr -d '\n')
rev=$("$YAMTK" seq -a subset -x "$tmpdir/rev.bed" -i "$TESTDIR/dna.fa" 2>/dev/null | grep -v '^>' | tr -d '\n')
fwd_rc=$(printf ">x\n%s\n" "$fwd" | "$YAMTK" seq -a rc -i - 2>/dev/null | grep -v '^>' | tr -d '\n')
if [ "$rev" != "$fwd_rc" ]; then
    FAIL "subset '-' strand differs from RC of forward subset: rev='$rev' fwd_rc='$fwd_rc'"
fi
PASS "subset '-' strand equals RC of forward slice"
