#!/usr/bin/env bash
# Audit: with -n 0 the written nsites= was the PPM rebuild scale (1000), not a
# site count.
source "$TESTDIR/lib.sh"
# -t 0.05 so the seed has hits on the small fixture (16 at this threshold).
out=$("$YAMTK" ref -m "$TESTDIR/motif.meme" -i "$TESTDIR/dna.fa" -n 0 -t 0.05 2>/dev/null | grep -oE 'nsites= *[0-9]+')
if [ -z "$out" ]; then FAIL "ref -n 0: no nsites= in output"; fi
if echo "$out" | grep -qE 'nsites= *1000$'; then FAIL "ref -n 0: nsites= is still the rebuild scale (1000)"; fi
PASS "ref -n 0 writes a real site count ($out)"
