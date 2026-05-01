#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

gz=$(mktemp /tmp/yamtk_test_XXXXXX).fa.gz
trap 'rm -f "$gz"' EXIT
gzip -c "$TESTDIR/dna.fa" > "$gz"

# gzipped and plain should produce identical output
plain=$("$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 1 2>/dev/null | grep -v '^##yamscan ')
gzout=$("$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$gz"               -t 0.05 -j 1 2>/dev/null | grep -v '^##yamscan ')

if [ "$plain" != "$gzout" ]; then
    diff <(echo "$plain") <(echo "$gzout") >&2
    echo "not ok - gzipped vs plain output differs" >&2; exit 1
fi
echo "ok - gzipped sequences produce identical output"
