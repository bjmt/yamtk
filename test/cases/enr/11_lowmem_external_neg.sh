#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -l with -n must produce byte-identical output to the default mode (modulo
# the ##yamenr command-line header line and the ##MotifCount header line,
# which assert_enr_golden strips).
tmpdir=$(mktemp -d /tmp/yamenr_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" enr    -i "$TESTDIR/fixtures/enr_pos.fa" -n "$TESTDIR/fixtures/enr_neg.fa" \
                -m "$TESTDIR/fixtures/enr_motifs.meme" -T seqs 2>/dev/null \
  | grep -v '^##yamenr \|^##MotifCount' > "$tmpdir/default.tsv"
"$YAMTK" enr -l -i "$TESTDIR/fixtures/enr_pos.fa" -n "$TESTDIR/fixtures/enr_neg.fa" \
                -m "$TESTDIR/fixtures/enr_motifs.meme" -T seqs 2>/dev/null \
  | grep -v '^##yamenr \|^##MotifCount' > "$tmpdir/lowmem.tsv"
if ! diff -u "$tmpdir/default.tsv" "$tmpdir/lowmem.tsv" >/dev/null; then
    diff -u "$tmpdir/default.tsv" "$tmpdir/lowmem.tsv" >&2 || true
    FAIL "-l with -n diverges from default mode (-T seqs)"
fi
PASS "-l with -n matches default mode (-T seqs)"

# Same check with -T sites
"$YAMTK" enr    -i "$TESTDIR/fixtures/enr_pos.fa" -n "$TESTDIR/fixtures/enr_neg.fa" \
                -m "$TESTDIR/fixtures/enr_motifs.meme" -T sites 2>/dev/null \
  | grep -v '^##yamenr \|^##MotifCount' > "$tmpdir/default_sites.tsv"
"$YAMTK" enr -l -i "$TESTDIR/fixtures/enr_pos.fa" -n "$TESTDIR/fixtures/enr_neg.fa" \
                -m "$TESTDIR/fixtures/enr_motifs.meme" -T sites 2>/dev/null \
  | grep -v '^##yamenr \|^##MotifCount' > "$tmpdir/lowmem_sites.tsv"
if ! diff -u "$tmpdir/default_sites.tsv" "$tmpdir/lowmem_sites.tsv" >/dev/null; then
    diff -u "$tmpdir/default_sites.tsv" "$tmpdir/lowmem_sites.tsv" >&2 || true
    FAIL "-l with -n diverges from default mode (-T sites)"
fi
PASS "-l with -n matches default mode (-T sites)"
