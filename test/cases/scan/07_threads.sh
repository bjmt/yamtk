#!/usr/bin/env bash
# Multi-threaded output (sorted) must match single-threaded output
source "$TESTDIR/lib.sh"

single=$("$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 1 2>/dev/null | grep -v '^##yamscan ' | sort)
multi=$( "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 4 2>/dev/null | grep -v '^##yamscan ' | sort)

if [ "$single" != "$multi" ]; then
    diff <(echo "$single") <(echo "$multi") >&2
    echo "not ok - -j 1 vs -j 4 output differs after sort" >&2; exit 1
fi
echo "ok - -j 1 and -j 4 produce same hits (modulo order)"
