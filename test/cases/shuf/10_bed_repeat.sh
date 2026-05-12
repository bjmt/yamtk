#!/usr/bin/env bash
source "$TESTDIR/lib.sh"

# -r 2 with -x: 3 total reps (0, -1, -2) per sequence; 3 seqs = 9 FASTA headers.
out=$("$YAMTK" shuf -i "$TESTDIR/dna.fa" -x "$TESTDIR/dna.bed" -k 1 -s 42 -r 2 2>/dev/null)

n_headers=$(echo "$out" | grep -c '^>' || true)
if [ "$n_headers" -ne 9 ]; then
    FAIL "shuf -x -r 2: expected 9 FASTA headers, got $n_headers"
fi
PASS "shuf -x -r 2 produces 9 FASTA headers (3 reps * 3 seqs)"

if ! echo "$out" | grep -q '^>2-1$'; then
    FAIL "shuf -x -r 2: missing repeat-1 suffix"
fi
PASS "shuf -x -r 2 has -1 suffix"

if ! echo "$out" | grep -q '^>2-2$'; then
    FAIL "shuf -x -r 2: missing repeat-2 suffix"
fi
PASS "shuf -x -r 2 has -2 suffix"
