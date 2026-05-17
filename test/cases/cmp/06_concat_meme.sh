#!/usr/bin/env bash
# Test 06: concatenated MEME files (multiple chunks with repeated headers and
# 'URL' trailer lines, like JASPAR2026_CORE_*_meme.txt) parse successfully,
# emit one warning per repeated header type, and produce both motifs as targets.
source "${TESTDIR}/lib.sh"

# Run capturing stderr separately. Use -q 1 so the test only checks parsing,
# not the (db-dependent) p-values from such a tiny target file.
OUT=$("$YAMTK" cmp -m "${TESTDIR}/fixtures/cmp_query.meme" \
                   -t "${TESTDIR}/fixtures/concat_meme.meme" -q 1 2>/tmp/yamtk_cmp_06.err)
STDERR=$(cat /tmp/yamtk_cmp_06.err); rm -f /tmp/yamtk_cmp_06.err

# Both motifs from the concatenated file should appear as targets.
for tid in chunk1_motif chunk2_motif; do
    if ! echo "$OUT" | awk -F'\t' -v t="$tid" '$2==t {found=1} END{exit !found}'; then
        FAIL "concat-meme: target $tid missing from output"
    fi
done
PASS "concat-meme: both motifs parsed from concatenated MEME file"

# URL line is not interpreted as a PPM row (would produce a parse error otherwise).
if echo "$STDERR" | grep -q -i "Failed to parse"; then
    FAIL "concat-meme: parser tried to read URL line as PPM data"
fi
PASS "concat-meme: URL trailer line skipped"

# Each warning type fires exactly once.
N_BKG_WARN=$(echo "$STDERR" | grep -c "Multiple 'Background letter frequencies' lines")
N_ALPH_WARN=$(echo "$STDERR" | grep -c "Multiple ALPHABET lines")
if [ "$N_BKG_WARN" -ne 1 ]; then
    FAIL "concat-meme: expected 1 background warning, got $N_BKG_WARN"
fi
if [ "$N_ALPH_WARN" -ne 1 ]; then
    FAIL "concat-meme: expected 1 alphabet warning, got $N_ALPH_WARN"
fi
PASS "concat-meme: each duplicate-header warning emitted exactly once"
