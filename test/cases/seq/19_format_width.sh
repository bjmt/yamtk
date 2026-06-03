#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -l <int> sets the wrap width: no sequence line exceeds it.
max80=$("$YAMTK" seq -a format -l 80 -i "$TESTDIR/dna.fa" 2>/dev/null \
    | grep -v '^>' | awk '{ if (length($0) > m) m = length($0) } END { print m }')
if [ "$max80" -gt 80 ]; then
    FAIL "seq format -l 80: found a line longer than 80 ($max80)"
fi
PASS "seq format -l 80 wraps at <= 80"

# -l 0 emits each sequence on a single line: one sequence line per record.
nseq=$(grep -c '^>' "$TESTDIR/dna.fa")
nlines=$("$YAMTK" seq -a format -l 0 -i "$TESTDIR/dna.fa" 2>/dev/null \
    | grep -vc '^>')
if [ "$nlines" -ne "$nseq" ]; then
    FAIL "seq format -l 0: expected $nseq sequence lines, got $nlines"
fi
PASS "seq format -l 0 emits one line per sequence"
