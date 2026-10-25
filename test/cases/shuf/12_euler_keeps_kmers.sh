#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Regression test: Euler shuffling must keep every k-mer count exactly.
# (The exit-edge tree used to be rooted at the last (k-1)-mer instead of the
# start of the reserved last k-mer, which changed counts and could hang.)

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

run() { if command -v timeout >/dev/null 2>&1; then timeout 30 "$@"; else "$@"; fi; }

cp "$TESTDIR/dna.fa" "$tmp/in.fa"
awk 'BEGIN { srand(7); split("A C G T", b, " ");
  for (i = 1; i <= 200; i++) { L = 20 + int(rand() * 180); s = "";
    for (j = 0; j < L; j++) s = s b[1 + int(rand() * 4)];
    print ">r" i; print s } }' >> "$tmp/in.fa"

for k in 2 3 4 5; do
  "$YAMTK" shuf -p -k "$k" -i "$tmp/in.fa" > "$tmp/before.tsv"
  for seed in 1 2 42; do
    run "$YAMTK" shuf -k "$k" -s "$seed" -i "$tmp/in.fa" > "$tmp/shuf.fa" \
      || FAIL "shuf Euler k=$k seed=$seed: exited non-zero or timed out"
    "$YAMTK" shuf -p -k "$k" -i "$tmp/shuf.fa" > "$tmp/after.tsv"
    if ! diff -q "$tmp/before.tsv" "$tmp/after.tsv" >/dev/null; then
      FAIL "shuf Euler k=$k seed=$seed: k-mer counts changed"
    fi
  done
done
PASS "shuf Euler keeps k-mer counts exactly (k=2..5, 3 seeds, 203 sequences)"
