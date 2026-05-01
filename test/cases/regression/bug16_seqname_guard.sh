#!/usr/bin/env bash
# Bug #16: add_seq_name wrote into `name[]` before checking the length guard.
# The buffer is correctly sized so no overflow happened in practice, but the
# guard was dead code. Fix: move the check above the copy loop.
# Test: a FASTA with a >512-char sequence name must exit 1 cleanly with a
# diagnostic — not silently truncate and not crash under ASan.
source "$TESTDIR/lib.sh"

fa=$(mktemp /tmp/yamtk_XXXXXX)
trap 'rm -f "$fa"' EXIT

# Build a 600-char sequence name
long_name=$(python3 -c "print('A' * 600)")
printf '>%s\nACGTACGT\n' "$long_name" > "$fa"

assert_exit_nonzero "scan with 600-char seq name exits non-zero" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$fa" -t 0.05 -j 1

assert_stderr_contains "scan with 600-char seq name prints size error" \
    "too large" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$fa" -t 0.05 -j 1
