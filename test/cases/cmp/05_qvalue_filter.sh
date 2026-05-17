#!/usr/bin/env bash
# Test 05: -q filters rows by q-value threshold.
source "${TESTDIR}/lib.sh"

# With -q 1 (effectively no filter): should see all 8 pairs (2 queries x 4 targets).
ALL=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/cmp_targets.meme" -q 1 2>/dev/null \
      | grep -v '^##\|^Query_ID' | wc -l | tr -d ' ')
if [ "$ALL" -ne 8 ]; then
    FAIL "filter: expected 8 rows with -q 1, got $ALL"
fi
PASS "filter: 8 rows with -q 1"

# With default filter (-Q 0.01): only strong matches survive.
DEF=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/cmp_targets.meme" 2>/dev/null \
      | grep -v '^##\|^Query_ID' | wc -l | tr -d ' ')
if [ "$DEF" -lt 2 ] || [ "$DEF" -gt 5 ]; then
    FAIL "filter: expected ~3 rows with default -q 0.01, got $DEF"
fi
PASS "filter: default -q 0.01 keeps $DEF (~3) strong hits"

# Tighter than default (-Q 0): exactly zero rows survive.
ZERO=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                    -t "${TESTDIR}/fixtures/cmp_targets.meme" -q 0 2>/dev/null \
       | awk '!/^##/' | wc -l | tr -d ' ')
if [ "$ZERO" -ne 0 ]; then
    FAIL "filter: expected 0 rows with -q 0, got $ZERO"
fi
PASS "filter: -q 0 keeps 0 rows"
