#!/usr/bin/env bash
# Bug #14: verbose BED log used motif loop var `i` instead of bed region var `k`
# for bed.starts[]/bed.ends[], which read garbage when i >= bed.n_regions.
# Fix: bed.starts[k] + 1, bed.ends[k].
# Test: every "Scanning range: A-B" line in stderr must match a real interval
# from the BED file.
source "$TESTDIR/lib.sh"

ranges=$("$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" \
    -x "$TESTDIR/dna.bed" -t 0.05 -j 1 -w 2>&1 >/dev/null \
    | grep 'Scanning range:' | sed 's/.*Scanning range: //')

bed_starts=$(awk '{print $2+1}' "$TESTDIR/dna.bed")
bed_ends=$(awk '{print $3}' "$TESTDIR/dna.bed")

# Each printed A-B must match exactly one BED interval's (start+1)-(end)
bad=0
while IFS=- read -r a b; do
    matched=0
    while IFS= read -r s && IFS= read -r e <&3; do
        if [ "$a" = "$s" ] && [ "$b" = "$e" ]; then matched=1; break; fi
    done < <(echo "$bed_starts") 3< <(echo "$bed_ends")
    if [ "$matched" -eq 0 ]; then bad=1; fi
done <<< "$ranges"

if [ "$bad" -eq 1 ]; then
    FAIL "verbose BED range log printed wrong coordinates"
fi
PASS "verbose BED scan ranges match dna.bed intervals"
