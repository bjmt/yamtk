#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
tmpin=$(mktemp /tmp/yamconv_test.XXXXXX)
trap 'rm -f "$tmpin"' EXIT
# Width-8 motif with two flanking low-IC columns on each side of a CACGTG core.
cat > "$tmpin" <<'EOF'
MEME version 4

ALPHABET= ACGT

MOTIF padded
letter-probability matrix: alength= 4 w= 8 nsites= 100 E= 0
0.25 0.25 0.25 0.25
0.02 0.96 0.01 0.01
0.96 0.01 0.01 0.02
0.01 0.96 0.01 0.02
0.01 0.01 0.96 0.02
0.01 0.96 0.02 0.01
0.25 0.25 0.25 0.25
0.25 0.25 0.25 0.25
EOF

w_before=$("$YAMTK" conv -m "$tmpin" -t meme 2>/dev/null | awk '/letter-probability/ { for (i=1; i<=NF; i++) if ($i=="w=") { print $(i+1); exit } }')
w_after=$("$YAMTK" conv -m "$tmpin" -t meme -T 0.5 2>/dev/null | awk '/letter-probability/ { for (i=1; i<=NF; i++) if ($i=="w=") { print $(i+1); exit } }')

if [ "$w_before" != "8" ]; then
    FAIL "expected pre-trim width 8, got $w_before"
fi
if [ "$w_after" -ge "$w_before" ]; then
    FAIL "IC trim did not shrink motif (pre=$w_before, post=$w_after)"
fi
PASS "IC trim shrinks motif from $w_before to $w_after columns"
