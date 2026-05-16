#!/usr/bin/env bash
# Test 03: yamcmp is deterministic — repeat runs of the same command produce
# byte-identical output (no RNG anywhere).
source "${TESTDIR}/lib.sh"

OUT1=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                    -t "${TESTDIR}/fixtures/cmp_targets.meme" 2>/dev/null \
       | grep -v '^##yamcmp ')
OUT2=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                    -t "${TESTDIR}/fixtures/cmp_targets.meme" 2>/dev/null \
       | grep -v '^##yamcmp ')

if [ "$OUT1" != "$OUT2" ]; then
    FAIL "determinism: repeat runs produced different output"
fi
PASS "determinism: repeat runs (empirical mode) produce identical output"

# Same check under -e (enum mode).
OUT3=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                    -t "${TESTDIR}/fixtures/cmp_targets.meme" -e 2>/dev/null \
       | grep -v '^##yamcmp ')
OUT4=$("$YAMTK" cmp -q "${TESTDIR}/fixtures/cmp_query.meme" \
                    -t "${TESTDIR}/fixtures/cmp_targets.meme" -e 2>/dev/null \
       | grep -v '^##yamcmp ')
if [ "$OUT3" != "$OUT4" ]; then
    FAIL "determinism: repeat -e runs produced different output"
fi
PASS "determinism: repeat -e runs produce identical output"
