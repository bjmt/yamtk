#!/usr/bin/env bash
# Bug #8: negative integers were accepted for -j (scan) and -k/-r/-s (shuf),
# causing signed-to-unsigned promotion yielding huge values.
source "$TESTDIR/lib.sh"

# scan -j -1
if "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -j -1 >/dev/null 2>&1; then
    echo "not ok - scan -j -1 should exit non-zero" >&2; exit 1
fi
echo "ok - scan -j -1 is rejected"

# shuf -k -1
if "$YAMTK" shuf -i "$TESTDIR/dna.fa" -k -1 >/dev/null 2>&1; then
    echo "not ok - shuf -k -1 should exit non-zero" >&2; exit 1
fi
echo "ok - shuf -k -1 is rejected"

# shuf -r -1
if "$YAMTK" shuf -i "$TESTDIR/dna.fa" -r -1 >/dev/null 2>&1; then
    echo "not ok - shuf -r -1 should exit non-zero" >&2; exit 1
fi
echo "ok - shuf -r -1 is rejected"

# shuf -s -1
if "$YAMTK" shuf -i "$TESTDIR/dna.fa" -s -1 >/dev/null 2>&1; then
    echo "not ok - shuf -s -1 should exit non-zero" >&2; exit 1
fi
echo "ok - shuf -s -1 is rejected"

# shuf -r INT_MAX (overflow in shuf_repeats + 1)
if "$YAMTK" shuf -i "$TESTDIR/dna.fa" -r 2147483647 >/dev/null 2>&1; then
    echo "not ok - shuf -r 2147483647 should exit non-zero" >&2; exit 1
fi
echo "ok - shuf -r INT_MAX is rejected"
