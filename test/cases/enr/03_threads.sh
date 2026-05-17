#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
# Threading can reorder output; sort data rows before diffing.
actual=$("$YAMTK" enr -q 1.0 \
    -i "$TESTDIR/fixtures/enr_pos.fa" \
    -n "$TESTDIR/fixtures/enr_neg.fa" \
    -m "$TESTDIR/fixtures/enr_motifs.meme" \
    -t 5e-4 \
    -j 2 2>/dev/null \
    | grep -v '^##yamenr \|^##MotifCount')
expected=$(cat "$TESTDIR/expected/enr_basic.txt")
if ! diff <(echo "$expected" | sort) <(echo "$actual" | sort) >/dev/null 2>&1; then
    diff <(echo "$expected" | sort) <(echo "$actual" | sort) >&2 || true
    FAIL "enr threads -j2 matches single-thread output"
fi
PASS "enr threads -j2 matches single-thread output"
