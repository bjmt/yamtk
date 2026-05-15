#!/usr/bin/env bash
# Test 12: -I tunes the IC threshold for auto-extend stopping and IC trim
source "${TESTDIR}/lib.sh"

# A low -I should let auto-extend grow wider than the default (0.5) before stopping.
# Default stops at width 8 on this fixture; -I 0.05 walks several more steps.
default_w=$("$YAMTK" ref -1 CGTACG -E -i "${TESTDIR}/fixtures/me_implanted.fa" 2>/dev/null \
            | awk '/^letter-probability matrix/ {print $7}')
lowered_w=$("$YAMTK" ref -1 CGTACG -E -I 0.05 -i "${TESTDIR}/fixtures/me_implanted.fa" 2>/dev/null \
            | awk '/^letter-probability matrix/ {print $7}')

if [ -z "$default_w" ] || [ -z "$lowered_w" ]; then
    FAIL "could not parse refined widths (default=$default_w lowered=$lowered_w)"
fi
if [ "$lowered_w" -le "$default_w" ]; then
    FAIL "expected -I 0.05 to extend further than default; got default=$default_w lowered=$lowered_w"
fi
PASS "-I 0.05 extends further than default (default=$default_w lowered=$lowered_w)"

# Raising -I should not extend further than the default.
raised_w=$("$YAMTK" ref -1 CGTACG -E -I 1.5 -i "${TESTDIR}/fixtures/me_implanted.fa" 2>/dev/null \
           | awk '/^letter-probability matrix/ {print $7}')
if [ "$raised_w" -gt "$default_w" ]; then
    FAIL "expected -I 1.5 to extend no further than default; got default=$default_w raised=$raised_w"
fi
PASS "-I 1.5 does not extend past default (default=$default_w raised=$raised_w)"

assert_exit_nonzero "-I out of range (negative)" \
  "$YAMTK" ref -1 CGTACG -E -I -1 -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_exit_nonzero "-I out of range (>2)" \
  "$YAMTK" ref -1 CGTACG -E -I 3 -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_stderr_contains "-I error message names the range" "must be in \\[0, 2\\]" \
  "$YAMTK" ref -1 CGTACG -E -I 3 -i "${TESTDIR}/fixtures/me_implanted.fa"
