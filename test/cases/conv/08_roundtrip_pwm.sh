#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# MEME -> JASPAR -> MEME, then check PPM values match the original to within
# pseudocount/rounding tolerance. Compares all numeric tokens line by line via
# awk so trivial formatting noise doesn't matter.
tmp1=$(mktemp /tmp/yamconv_test.XXXXXX)
tmp2=$(mktemp /tmp/yamconv_test.XXXXXX)
trap 'rm -f "$tmp1" "$tmp2"' EXIT
"$YAMTK" conv -m "$TESTDIR/fixtures/enr_motifs.meme" -t jaspar > "$tmp1" 2>/dev/null
"$YAMTK" conv -m "$tmp1" -t meme > "$tmp2" 2>/dev/null
# Extract probability matrices from both and compare with epsilon = 0.005
extract_pwm() {
    awk '/^ / && NF==4 {print $1, $2, $3, $4}' "$1"
}
orig=$(extract_pwm "$TESTDIR/fixtures/enr_motifs.meme")
rt=$(extract_pwm "$tmp2")
n=$(echo "$orig" | wc -l | tr -d ' ')
n_rt=$(echo "$rt" | wc -l | tr -d ' ')
if [ "$n" != "$n_rt" ]; then
    FAIL "roundtrip: row count differs ($n vs $n_rt)"
fi
# Tolerance 0.02: rounding at nsites=100 plus the source rows summing slightly
# over 1.00 (e.g. 0.01+0.01+0.96+0.03=1.01) can move a single column by ~1/N.
set +e
err=$(paste <(echo "$orig") <(echo "$rt") | awk '
{
    for (i=1; i<=4; i++) {
        d = $i - $(i+4); if (d < 0) d = -d;
        if (d > 0.02) { print "diff too large at col "i" (row "NR"): "$i" vs "$(i+4); exit 1 }
    }
}')
rc=$?
set -e
if [ "$rc" -ne 0 ]; then
    echo "$err" >&2
    FAIL "roundtrip MEME -> JASPAR -> MEME drifted > 0.02"
fi
PASS "roundtrip MEME -> JASPAR -> MEME within tolerance"
