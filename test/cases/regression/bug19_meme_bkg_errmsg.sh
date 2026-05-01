#!/usr/bin/env bash
# Bug #19: get_meme_bkg error messages referenced 'C' while parsing G and T/U.
# Fix: change error strings to name the correct character ('G', 'T/U').
# Test: MEME files with missing whitespace before G/T in the background block
# must produce a diagnostic naming the right character.
source "$TESTDIR/lib.sh"

meme_g=$(mktemp /tmp/yamtk_XXXXXX)
meme_t=$(mktemp /tmp/yamtk_XXXXXX)
trap 'rm -f "$meme_g" "$meme_t"' EXIT

# No space between the C value and G token triggers the whitespace error for G
cat > "$meme_g" << 'EOF'
MEME version 4

ALPHABET= ACGT

Background letter frequencies
A 0.30 C 0.20G 0.25 T 0.25
EOF

# No space between the G value and T token triggers the whitespace error for T/U
cat > "$meme_t" << 'EOF'
MEME version 4

ALPHABET= ACGT

Background letter frequencies
A 0.30 C 0.20 G 0.25T 0.25
EOF

assert_stderr_contains "bad G background line mentions 'G' not 'C'" \
    "whitespace before 'G'" \
    "$YAMTK" scan -m "$meme_g" -s "$TESTDIR/dna.fa" -t 0.05 -j 1

assert_stderr_contains "bad T background line mentions 'T/U' not 'C'" \
    "whitespace before 'T/U'" \
    "$YAMTK" scan -m "$meme_t" -s "$TESTDIR/dna.fa" -t 0.05 -j 1
