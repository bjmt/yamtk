#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# T -> U -> T round-trip preserves sequence content.
tmpdir=$(mktemp -d /tmp/yamseq_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seq -a rna -i "$TESTDIR/dna.fa" 2>/dev/null \
  | "$YAMTK" seq -a dna -i - 2>/dev/null \
  | awk '/^>/{if(s)print s; print $1; s="";next} {s=s$0} END{if(s)print s}' \
  > "$tmpdir/rt.fa"
awk '/^>/{if(s)print s; print $1; s="";next} {s=s$0} END{if(s)print s}' \
  "$TESTDIR/dna.fa" > "$tmpdir/in.fa"
if ! diff -u "$tmpdir/in.fa" "$tmpdir/rt.fa" >/dev/null 2>&1; then
    diff -u "$tmpdir/in.fa" "$tmpdir/rt.fa" >&2 || true
    FAIL "rna dna round-trip differs from input"
fi
PASS "rna dna round-trip is identity"

# Also: rna output must contain no T (it should be U).
if "$YAMTK" seq -a rna -i "$TESTDIR/dna.fa" 2>/dev/null | grep -v '^>' | grep -q 'T'; then
    FAIL "rna output still contains uppercase T"
fi
PASS "rna output contains no uppercase T"
