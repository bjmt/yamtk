#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Build a 1000-sequence fixture so -f sampling is statistically meaningful.
fa=$(mktemp -t yamseq_subXXXXXX)
trap 'rm -f "$fa"' EXIT
for i in $(seq 1 1000); do
    printf '>seq%d\nACGT\n' "$i"
done > "$fa"

# -f 0.5 -s 42 should keep close to 500. Allow generous slack (450..550).
out=$("$YAMTK" seq -a subsample -f 0.5 -s 42 -i "$fa" 2>/dev/null)
n=$(echo "$out" | grep -c '^>')
if [ "$n" -lt 450 ] || [ "$n" -gt 550 ]; then
    FAIL "subsample -f 0.5 keeps ~500/1000 (got $n)"
fi
PASS "subsample -f 0.5 keeps ~half"

# Bernoulli output must already be in input order (we never buffer).
ordered=$(echo "$out" | awk '/^>/ {print substr($1,5)+0}' | awk 'NR==1 || $1>prev {prev=$1; ok=1; next} {ok=0; exit} END {print ok}')
if [ "$ordered" != "1" ]; then
    FAIL "subsample -f preserves input order"
fi
PASS "subsample -f preserves input order"
