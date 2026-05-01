#!/usr/bin/env bash
# Bug #12: shuffle_linear was biased and skipped tail characters.
# Test 1: letter histogram is preserved (no characters lost or gained).
# Test 2: output is deterministic with a fixed seed.
# Fixtures: seq_22 (22 bp, 22%3=1 residual) and seq_21 (21 bp, multiple of 3).
source "$TESTDIR/lib.sh"

out1=$(mktemp /tmp/yamtk_lin_XXXXXX)
out2=$(mktemp /tmp/yamtk_lin_XXXXXX)
trap 'rm -f "$out1" "$out2"' EXIT

"$YAMTK" shuf -i "$TESTDIR/fixtures/linear_shuf.fa" -l -k 3 -s 42 2>/dev/null > "$out1"
"$YAMTK" shuf -i "$TESTDIR/fixtures/linear_shuf.fa" -l -k 3 -s 42 2>/dev/null > "$out2"

# Determinism
if ! diff -q "$out1" "$out2" >/dev/null; then
    echo "not ok - linear shuffle is non-deterministic with fixed seed" >&2; exit 1
fi
echo "ok - linear shuffle is deterministic"

# Letter histogram preserved for each sequence
for seq_name in seq_22 seq_21; do
    in_hist=$(awk -v s=">$seq_name" '$0==s{f=1;next} f&&/^>/{exit} f' \
        "$TESTDIR/fixtures/linear_shuf.fa" | tr -d '\n' | fold -w1 | sort | uniq -c)
    out_hist=$(awk -v s=">$seq_name" '$0==s{f=1;next} f&&/^>/{exit} f' \
        "$out1" | tr -d '\n' | fold -w1 | sort | uniq -c)
    if [ "$in_hist" != "$out_hist" ]; then
        echo "not ok - linear shuffle changed letter histogram for $seq_name" >&2
        echo "  input:  $in_hist" >&2
        echo "  output: $out_hist" >&2; exit 1
    fi
    echo "ok - linear shuffle preserves letter histogram for $seq_name"
done
