#!/usr/bin/env bash
# Test 07: invalid arguments exit non-zero with informative stderr
source "${TESTDIR}/lib.sh"

assert_exit_nonzero "no -m flag" \
  "$YAMTK" ref -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_stderr_contains "no -m stderr" "Missing required argument -m" \
  "$YAMTK" ref -i "${TESTDIR}/fixtures/me_implanted.fa"

assert_exit_nonzero "no -i flag" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme"

assert_stderr_contains "no -i stderr" "Missing required argument -i" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme"

assert_exit_nonzero "extension exceeds max width" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 30

assert_stderr_contains "oversized extension stderr" "max" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa" -e 30

assert_exit_nonzero "bad -i path" \
  "$YAMTK" ref -m "${TESTDIR}/fixtures/ref_seed.meme" -i /nonexistent.fa

assert_exit_nonzero "negative -n" \
  "$YAMTK" ref -n -1 -m "${TESTDIR}/fixtures/ref_seed.meme" \
  -i "${TESTDIR}/fixtures/me_implanted.fa"
