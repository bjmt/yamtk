#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Statistical check: -c 20 must make insertion midpoints cluster
# meaningfully closer to the sequence centre than -c 1.
tmpdir=$(mktemp -d /tmp/yamseed_c_test.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

# Generate a 100-seq x 200-bp synthetic FASTA so we have enough samples
# to compare distributions stably (motif width 6).
python3 - > "$tmpdir/in.fa" <<'PY'
import random
random.seed(11)
b = 'ACGT'
for i in range(200):
    print(f'>s{i}')
    print(''.join(random.choice(b) for _ in range(200)))
PY

"$YAMTK" seed -1 CACGTG -i "$tmpdir/in.fa" \
              -f 0.05 -s 7 -c 1 -O "$tmpdir/truth_c1.bed" >/dev/null 2>&1
"$YAMTK" seed -1 CACGTG -i "$tmpdir/in.fa" \
              -f 0.05 -s 7 -c 20 -O "$tmpdir/truth_c20.bed" >/dev/null 2>&1

# Mean |midpoint - L/2| over the truth BED. With L=200, baseline uniform
# expectation is around 50; -c 20 should be well below that.
mean_d() { awk '{m=($2+$3)/2; d=m-100; if(d<0)d=-d; s+=d; n++} END{if(n>0)printf "%.3f", s/n; else print "0"}' "$1"; }
d1=$(mean_d "$tmpdir/truth_c1.bed")
d20=$(mean_d "$tmpdir/truth_c20.bed")
# Compare with awk to avoid bc dependency.
ok=$(awk -v a="$d20" -v b="$d1" 'BEGIN{print (a < b*0.5) ? 1 : 0}')
if [ "$ok" -ne 1 ]; then
    echo "Expected -c 20 mean distance ($d20) to be < 50% of -c 1 mean distance ($d1)." >&2
    FAIL "seed -c 20 concentrates positions"
fi
PASS "seed -c 20 concentrates positions (mean dist: $d20 vs $d1)"
