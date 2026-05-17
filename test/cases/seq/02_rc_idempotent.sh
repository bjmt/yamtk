#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Double RC should give back the input at the sequence-content level
# (line wrap and headers may differ; compare linearized seqs only).
tmpdir=$(mktemp -d /tmp/yamseq_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seq -a rc -i "$TESTDIR/dna.fa" 2>/dev/null \
  | "$YAMTK" seq -a rc -i - 2>/dev/null \
  | awk '/^>/{if(s)print s; print $1; s="";next} {s=s$0} END{if(s)print s}' \
  > "$tmpdir/rt.fa"
awk '/^>/{if(s)print s; print $1; s="";next} {s=s$0} END{if(s)print s}' \
  "$TESTDIR/dna.fa" > "$tmpdir/in.fa"
if ! diff -u "$tmpdir/in.fa" "$tmpdir/rt.fa" >/dev/null 2>&1; then
    diff -u "$tmpdir/in.fa" "$tmpdir/rt.fa" >&2 || true
    FAIL "rc rc round-trip differs from input"
fi
PASS "rc rc round-trip is identity"
