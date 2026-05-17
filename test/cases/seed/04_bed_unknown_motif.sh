#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmpdir=$(mktemp -d /tmp/yamseed_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
cat > "$tmpdir/bad.bed" <<EOF
seq1	10	16	xyz_unknown	0	+
EOF
assert_exit_nonzero "seed -x with unknown motif name aborts" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
                  -i "$TESTDIR/fixtures/me_implanted.fa" \
                  -x "$tmpdir/bad.bed"
assert_stderr_contains "seed -x unknown motif error message" \
    "no matching motif" \
    "$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
                  -i "$TESTDIR/fixtures/me_implanted.fa" \
                  -x "$tmpdir/bad.bed"
