#!/usr/bin/env bash
# Test 11: multi-motif input file is processed motif-by-motif, with
# under-represented motifs dropped and the survivors emitted.
# Explicit -t 1e-3 keeps gcbox under MIN_REFINE_HITS=20 (the default
# 1e-4 is unreachable for 6-bp motifs and falls back to p=0.01, which
# pulls in too many positions to trigger the under-represented drop).
source "${TESTDIR}/lib.sh"

# enr_motifs.meme contains ebox (CACGTG, ~104 hits in enr_pos.fa at -t 1e-3)
# and gcbox (GGGCGG, ~8 hits — below MIN_REFINE_HITS=20).
tmp="$(mktemp -d -t yamref-multi-XXXXXX)"
trap 'rm -rf "$tmp"' EXIT

assert_exit0 "multi-motif input completes" \
  "$YAMTK" ref -t 1e-3 -m "${TESTDIR}/fixtures/enr_motifs.meme" \
  -i "${TESTDIR}/fixtures/enr_pos.fa" -o "$tmp/out.meme"

assert_stderr_contains "ebox summary printed" "\\[ebox CACGTG\\]" \
  "$YAMTK" ref -t 1e-3 -m "${TESTDIR}/fixtures/enr_motifs.meme" \
  -i "${TESTDIR}/fixtures/enr_pos.fa" -o /dev/null

assert_stderr_contains "gcbox drop warning printed" "\\[gcbox GGGCGG\\] dropped" \
  "$YAMTK" ref -t 1e-3 -m "${TESTDIR}/fixtures/enr_motifs.meme" \
  -i "${TESTDIR}/fixtures/enr_pos.fa" -o /dev/null

# Output should contain only the surviving motif (ebox), not gcbox.
n_motifs=$(grep -c '^MOTIF ' "$tmp/out.meme")
if [ "$n_motifs" -ne 1 ]; then
    FAIL "expected 1 MOTIF in output, got $n_motifs"
fi
PASS "only surviving motif written"

if ! grep -q '^MOTIF ebox ' "$tmp/out.meme"; then
    FAIL "ebox missing from output"
fi
PASS "ebox is in output"

if grep -q '^MOTIF gcbox ' "$tmp/out.meme"; then
    FAIL "gcbox should have been dropped"
fi
PASS "gcbox not in output"
