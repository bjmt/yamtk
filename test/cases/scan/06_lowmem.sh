#!/usr/bin/env bash
# -l disables low-memory mode (loads all seqs into memory); output must match default
source "$TESTDIR/lib.sh"

default_out=$("$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 1 2>/dev/null | grep -v '^##yamscan ')
lowmem_off=$( "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 1 -l 2>/dev/null | grep -v '^##yamscan ')

if [ "$default_out" != "$lowmem_off" ]; then
    diff <(echo "$default_out") <(echo "$lowmem_off") >&2
    echo "not ok - low-mem mode vs in-memory output differs" >&2; exit 1
fi
echo "ok - default (low-mem) and -l (in-memory) produce identical output"
