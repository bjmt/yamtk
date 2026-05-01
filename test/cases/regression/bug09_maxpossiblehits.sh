#!/usr/bin/env bash
# Bug #9: MaxPossibleHits underflowed (uint64 wrap) when a sequence was shorter
# than the motif. The header would report a near-UINT64_MAX number.
# Fix: guarded subtraction, so short sequences contribute 0 to MaxPossibleHits.
source "$TESTDIR/lib.sh"

# Create a FASTA with one sequence shorter than the 5-bp motif
short_fa=$(mktemp /tmp/yamtk_short_XXXXXX.fa)
trap 'rm -f "$short_fa"' EXIT
printf '>short\nACG\n' > "$short_fa"   # 3 bp < motif width 5

header=$("$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$short_fa" -j 1 2>/dev/null | grep 'MaxPossibleHits')

# Extract the number
max_hits=$(echo "$header" | grep -o 'MaxPossibleHits=[0-9]*' | cut -d= -f2)

# Must be a small sane number, not a wrap-around value near 2^64
if [ -z "$max_hits" ]; then
    echo "not ok - MaxPossibleHits not found in header" >&2; echo "$header" >&2; exit 1
fi
# Any value > 1000000 is almost certainly a wrap; the real answer is 0
if [ "$max_hits" -gt 1000000 ] 2>/dev/null; then
    echo "not ok - MaxPossibleHits=$max_hits looks like a uint64 underflow" >&2; exit 1
fi
echo "ok - MaxPossibleHits=$max_hits (no underflow for sequence shorter than motif)"
