#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Regression test: every shuffling method should be able to reach every valid
# output. (k = 1 used to never move the last letter, -l could never leave a
# block in place, and the Euler method only reached sequences ending in the
# same last k-mer.) With 3,000 draws, missing even one output by chance is
# vanishingly unlikely.

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

check() {  # desc, sequence, expected number of distinct outputs, shuf args...
  local desc="$1" seq="$2" want="$3"; shift 3
  awk -v s="$seq" 'BEGIN { for (i = 0; i < 3000; i++) print ">s" i "\n" s }' > "$tmp/in.fa"
  local got
  got=$("$YAMTK" shuf -s 1 -i "$tmp/in.fa" "$@" | awk '!/^>/' | sort -u | awk 'END { print NR }')
  [ "$got" -eq "$want" ] || FAIL "$desc: reached $got of $want possible outputs"
}

check "Fisher-Yates (k=1)" ACGTN 120 -k 1
check "linear (-l, k=2)" AACCGGTT 24 -k 2 -l
check "Euler (k=2)" ACGTTGCAACGT 16 -k 2
check "Euler (k=3)" AACGTTACGGTA 4 -k 3
PASS "shuffles reach every valid output (k=1, -l, Euler k=2/3)"
