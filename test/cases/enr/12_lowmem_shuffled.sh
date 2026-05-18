#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -l with shuffled negatives must produce byte-identical output to the
# default mode for the same -k and -s.
tmpdir=$(mktemp -d /tmp/yamenr_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

for k in 1 2; do
  "$YAMTK" enr    -i "$TESTDIR/fixtures/enr_pos.fa" \
                  -m "$TESTDIR/fixtures/enr_motifs.meme" -k $k -s 42 2>/dev/null \
    | grep -v '^##yamenr \|^##MotifCount' > "$tmpdir/d_k$k.tsv"
  "$YAMTK" enr -l -i "$TESTDIR/fixtures/enr_pos.fa" \
                  -m "$TESTDIR/fixtures/enr_motifs.meme" -k $k -s 42 2>/dev/null \
    | grep -v '^##yamenr \|^##MotifCount' > "$tmpdir/l_k$k.tsv"
  if ! diff -u "$tmpdir/d_k$k.tsv" "$tmpdir/l_k$k.tsv" >/dev/null; then
      diff -u "$tmpdir/d_k$k.tsv" "$tmpdir/l_k$k.tsv" >&2 || true
      FAIL "-l shuffled-neg diverges from default (k=$k)"
  fi
  PASS "-l shuffled-neg matches default (k=$k -s 42)"
done
