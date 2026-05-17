#!/usr/bin/env bash
# Test 02: the non-palindromic `gata` query (AGATAA, RC=TTATCT) should hit
# `gata_rc_db` with orientation '-'. With default -R disabled (i.e. RC scanning
# on), the gata_db (forward) and gata_rc_db (reverse) targets should tie for
# best hit.
source "${TESTDIR}/lib.sh"

OUT=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/cmp_targets.meme" -q 1 2>/dev/null \
      | grep -v '^##\|^Query_ID')

# gata vs gata_rc_db row: must have orientation '-' and best-tier q-value.
GATA_RC=$(echo "$OUT" | awk -F'\t' '$1=="gata" && $2=="gata_rc_db"')
if [ -z "$GATA_RC" ]; then
    FAIL "rc-match: no gata vs gata_rc_db row"
fi
ORIENT=$(echo "$GATA_RC" | awk -F'\t' '{print $10}')
if [ "$ORIENT" != "-" ]; then
    FAIL "rc-match: gata vs gata_rc_db orientation $ORIENT (expected -)"
fi
PASS "rc-match: gata vs gata_rc_db is orientation '-'"

# With -R (forward only), gata_rc_db should NOT be a strong hit.
OUT_FWD=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                       -t "${TESTDIR}/fixtures/cmp_targets.meme" -R -q 1 2>/dev/null \
          | grep -v '^##\|^Query_ID')

GATA_RC_FWD_QVAL=$(echo "$OUT_FWD" | awk -F'\t' '$1=="gata" && $2=="gata_rc_db" {print $5}')
GATA_SELF_FWD_QVAL=$(echo "$OUT_FWD" | awk -F'\t' '$1=="gata" && $2=="gata_db" {print $5}')

# Compare numerically: forward-only gata_db should be MUCH better than gata_rc_db.
awk -v rc="$GATA_RC_FWD_QVAL" -v self="$GATA_SELF_FWD_QVAL" \
    'BEGIN { if (rc <= self) exit 1; else exit 0 }' \
    || FAIL "rc-match: with -R, gata vs gata_rc_db qval=$GATA_RC_FWD_QVAL not worse than self qval=$GATA_SELF_FWD_QVAL"
PASS "rc-match: -R correctly excludes reverse-complement targets"
