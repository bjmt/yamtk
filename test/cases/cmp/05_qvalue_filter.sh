#!/usr/bin/env bash
# Test 05: -Q filters rows by q-value threshold.
source "${TESTDIR}/lib.sh"

# With -Q 1 (effectively no filter): should see all 8 pairs (2 queries x 4 targets).
ALL=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/cmp_targets.meme" -Q 1 2>/dev/null \
      | grep -v '^##\|^Query_ID' | wc -l | tr -d ' ')
if [ "$ALL" -ne 8 ]; then
    FAIL "filter: expected 8 rows with -Q 1, got $ALL"
fi
PASS "filter: 8 rows with -Q 1"

# With default filter (-Q 0.01): only strong matches survive.
DEF=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/cmp_targets.meme" 2>/dev/null \
      | grep -v '^##\|^Query_ID' | wc -l | tr -d ' ')
if [ "$DEF" -lt 2 ] || [ "$DEF" -gt 5 ]; then
    FAIL "filter: expected ~3 rows with default -Q 0.01, got $DEF"
fi
PASS "filter: default -Q 0.01 keeps $DEF (~3) strong hits"

# Tighter than default (-Q 0): exactly zero rows survive.
ZERO=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                    -t "${TESTDIR}/fixtures/cmp_targets.meme" -Q 0 2>/dev/null \
       | awk '!/^##/' | wc -l | tr -d ' ')
if [ "$ZERO" -ne 0 ]; then
    FAIL "filter: expected 0 rows with -Q 0, got $ZERO"
fi
PASS "filter: -Q 0 keeps 0 rows"
