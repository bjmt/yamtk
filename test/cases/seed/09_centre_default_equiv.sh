#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# With -c 1 (or no -c), output must be byte-identical to the pre-flag
# behaviour for the same seed — guarantees no silent break of existing
# fixtures / user pipelines.
tmpdir=$(mktemp -d /tmp/yamseed_c_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
              -i "$TESTDIR/fixtures/me_implanted.fa" \
              -f 0.01 -s 1 -O "$tmpdir/truth_noc.bed" >/dev/null 2>&1
"$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
              -i "$TESTDIR/fixtures/me_implanted.fa" \
              -f 0.01 -s 1 -c 1 -O "$tmpdir/truth_c1.bed" >/dev/null 2>&1
if ! diff -q "$tmpdir/truth_noc.bed" "$tmpdir/truth_c1.bed" >/dev/null 2>&1; then
    diff -u "$tmpdir/truth_noc.bed" "$tmpdir/truth_c1.bed" >&2 || true
    FAIL "seed -c 1 must produce byte-identical output to no -c"
fi
PASS "seed -c 1 == no -c (determinism preserved)"
