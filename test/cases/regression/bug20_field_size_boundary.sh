#!/usr/bin/env bash
# Bug #20: check_field_size used `size > FIELD_MAX_CHAR` (512) so a field of
# exactly 512 chars was silently accepted but extract_field had already dropped
# it to empty (its guard is `size < FIELD_MAX_CHAR`).
# Fix: change to `size >= FIELD_MAX_CHAR` so a 512-char field is rejected.
# Test: dedup input with a 512-char seq_name field must exit non-zero.
source "$TESTDIR/lib.sh"

input=$(mktemp /tmp/yamtk_XXXXXX)
trap 'rm -f "$input"' EXIT

# Build a 512-char sequence name (exactly FIELD_MAX_CHAR)
seq_name=$(python3 -c "print('A' * 512, end='')")

# Minimal yamscan-style line: seq_name start end strand motif pvalue score score_pct match
printf '%s\t1\t10\t+\tmotif1\t0.001\t5.0\t80.0\tACGTACGTA\n' "$seq_name" > "$input"

assert_exit_nonzero "dedup with 512-char field exits non-zero" \
    "$YAMTK" dedup -i "$input"

assert_stderr_contains "dedup with 512-char field prints size error" \
    "exceeds max" \
    "$YAMTK" dedup -i "$input"
