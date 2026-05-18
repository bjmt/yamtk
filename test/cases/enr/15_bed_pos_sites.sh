#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# -x with -T sites: position denominator is the sum of region widths (* strands).
assert_enr_golden "enr -x positives (sites)" \
    "$TESTDIR/expected/enr_bed_pos_sites.txt" \
    "$YAMTK" enr -q 1.0 \
        -i "$TESTDIR/fixtures/enr_pos.fa" \
        -n "$TESTDIR/fixtures/enr_neg.fa" \
        -m "$TESTDIR/fixtures/enr_motifs.meme" \
        -t 5e-4 -T sites \
        -x "$TESTDIR/fixtures/enr_pos.bed"

# The pos_site_hits (column 6) for ebox under -x should be <= the full-seq value
# (we're only scanning 30bp of each 100bp seq).
hits_full=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -T sites 2>/dev/null | grep -v '^#' | awk '$1=="ebox"{print $6; exit}')
hits_bed=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 -T sites -x "$TESTDIR/fixtures/enr_pos.bed" 2>/dev/null | grep -v '^#' | awk '$1=="ebox"{print $6; exit}')
if [ "$hits_bed" -gt "$hits_full" ]; then
    FAIL "BED-restricted hits not <= full-seq hits (bed=$hits_bed, full=$hits_full)"
fi
PASS "BED-restricted sites <= full-seq sites"
