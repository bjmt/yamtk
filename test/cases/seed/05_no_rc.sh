#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmpdir=$(mktemp -d /tmp/yamseed_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT
"$YAMTK" seed -m "$TESTDIR/fixtures/enr_motifs.meme" \
              -i "$TESTDIR/fixtures/me_implanted.fa" \
              -f 0.05 -s 1 -R -O "$tmpdir/truth.bed" >/dev/null 2>&1
neg=$(awk '$6 == "-"' "$tmpdir/truth.bed" | wc -l | tr -d ' ')
if [ "$neg" != "0" ]; then
    FAIL "-R should produce zero '-' strand insertions, got $neg"
fi
PASS "-R produces zero '-' strand insertions"
