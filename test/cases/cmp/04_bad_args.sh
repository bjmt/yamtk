#!/usr/bin/env bash
# Test 04: missing required args exits non-zero with informative stderr.
source "${TESTDIR}/lib.sh"

assert_exit_nonzero "no -m flag" \
  "$YAMTK" cmp -t "${TESTDIR}/fixtures/cmp_targets.meme"

assert_stderr_contains "no -m stderr" "Missing -m" \
  "$YAMTK" cmp -t "${TESTDIR}/fixtures/cmp_targets.meme"

assert_exit_nonzero "no -t flag" \
  "$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme"

assert_stderr_contains "no -t stderr" "Missing -t" \
  "$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme"

assert_exit_nonzero "unknown metric" \
  "$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
               -t "${TESTDIR}/fixtures/cmp_targets.meme" -d bogus

assert_stderr_contains "unknown metric stderr" "Unknown metric" \
  "$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
               -t "${TESTDIR}/fixtures/cmp_targets.meme" -d bogus
