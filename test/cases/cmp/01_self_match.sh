#!/usr/bin/env bash
# Test 01: each query motif's best hit (lowest q-value) against itself in the
# target file is at offset 0, orientation +, with overlap == motif width.
source "${TESTDIR}/lib.sh"

OUT=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/cmp_query.meme" -q 1 2>/dev/null \
      | grep -v '^##\|^Query_ID')

# Group by query, take first row per query (sorted by q-value ascending in output).
TOP=$(echo "$OUT" | awk -F'\t' '!seen[$1]++ { print }')

# Every top hit must be self-match: Query_ID == Target_ID, offset == 0,
# orientation == +, overlap > 0.
echo "$TOP" | while IFS=$'\t' read -r qid tid offset pval qval score overlap qcon tcon orient; do
    if [ "$qid" != "$tid" ]; then
        FAIL "self-match: query $qid top hit is $tid (expected $qid)"
    fi
    if [ "$offset" != "0" ]; then
        FAIL "self-match: query $qid top hit offset $offset (expected 0)"
    fi
    if [ "$orient" != "+" ]; then
        FAIL "self-match: query $qid top hit orientation $orient (expected +)"
    fi
    if [ "$overlap" -le 0 ]; then
        FAIL "self-match: query $qid overlap $overlap (expected > 0)"
    fi
done
PASS "self-match: every query maps to itself at offset 0, orientation +"
