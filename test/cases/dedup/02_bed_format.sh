#!/usr/bin/env bash
# yamdedup also accepts BED input: 6 columns (name, start, end, motif, score, strand)
source "$TESTDIR/lib.sh"

bed_in=$(mktemp /tmp/yamtk_dedup_XXXXXX.bed)
trap 'rm -f "$bed_in"' EXIT

# Three overlapping hits on seq1/motif1/+ - the highest score (3.0) should survive
printf 'seq1\t1\t6\tmotif1\t3.0\t+\n' >  "$bed_in"
printf 'seq1\t2\t7\tmotif1\t2.0\t+\n' >> "$bed_in"
printf 'seq1\t3\t8\tmotif1\t1.0\t+\n' >> "$bed_in"

out=$("$YAMTK" dedup -i "$bed_in" 2>/dev/null)
count=$(echo "$out" | grep -c '^seq1' || true)

if [ "$count" -ne 1 ]; then
    echo "Output:" >&2; echo "$out" >&2
    echo "not ok - expected 1 survivor, got $count" >&2; exit 1
fi
winner=$(echo "$out" | grep '^seq1' | awk '{print $5}')
if [ "$winner" != "3.0" ]; then
    echo "not ok - expected winner score 3.0, got $winner" >&2; exit 1
fi
echo "ok - dedup BED format keeps highest-scoring hit"
