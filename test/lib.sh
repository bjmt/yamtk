#!/usr/bin/env bash
# Shared helpers for yamtk test cases.
# Source this file from every case script.
# Expects $YAMTK (path to binary) and $TESTDIR (path to test/) to be set.

set -euo pipefail

PASS() { echo "ok - $*"; }
FAIL() { echo "not ok - $*" >&2; exit 1; }

assert_exit0() {
    local desc="$1"; shift
    if ! "$@" >/dev/null 2>&1; then
        FAIL "$desc: expected exit 0, got non-zero"
    fi
    PASS "$desc exits 0"
}

assert_exit_nonzero() {
    local desc="$1"; shift
    if "$@" >/dev/null 2>&1; then
        FAIL "$desc: expected non-zero exit, got 0"
    fi
    PASS "$desc exits non-zero"
}

assert_stderr_contains() {
    local desc="$1"; local pattern="$2"; shift 2
    local err
    err=$("$@" 2>&1 >/dev/null || true)
    if ! echo "$err" | grep -q "$pattern"; then
        FAIL "$desc: stderr did not contain '$pattern'"
    fi
    PASS "$desc stderr contains '$pattern'"
}

assert_golden() {
    local desc="$1"; local expected="$2"; shift 2
    local actual
    actual=$("$@" 2>/dev/null | grep -v '^##yamscan ')
    if ! diff -u "$expected" <(echo "$actual") >/dev/null 2>&1; then
        diff -u "$expected" <(echo "$actual") >&2 || true
        FAIL "$desc: output differs from golden $expected"
    fi
    PASS "$desc matches golden"
}

assert_stdout_contains() {
    local desc="$1"; local pattern="$2"; shift 2
    local out
    out=$("$@" 2>/dev/null)
    if ! echo "$out" | grep -q "$pattern"; then
        FAIL "$desc: stdout did not contain '$pattern'"
    fi
    PASS "$desc stdout contains '$pattern'"
}

assert_line_count_ge() {
    local desc="$1"; local min="$2"; shift 2
    local n
    n=$("$@" 2>/dev/null | grep -v '^#' | wc -l | tr -d ' ')
    if [ "$n" -lt "$min" ]; then
        FAIL "$desc: expected >= $min data lines, got $n"
    fi
    PASS "$desc has >= $min data lines"
}

assert_enr_golden() {
    local desc="$1"; local expected="$2"; shift 2
    local actual
    actual=$("$@" 2>/dev/null | grep -v '^##yamenr \|^##MotifCount')
    if ! diff -u "$expected" <(echo "$actual") >/dev/null 2>&1; then
        diff -u "$expected" <(echo "$actual") >&2 || true
        FAIL "$desc: output differs from golden $expected"
    fi
    PASS "$desc matches golden"
}

assert_me_golden() {
    local desc="$1"; local expected="$2"; shift 2
    local actual
    actual=$("$@" 2>/dev/null | grep -v '^##yamme ')
    if ! diff -u "$expected" <(echo "$actual") >/dev/null 2>&1; then
        diff -u "$expected" <(echo "$actual") >&2 || true
        FAIL "$desc: output differs from golden $expected"
    fi
    PASS "$desc matches golden"
}

assert_ref_golden() {
    local desc="$1"; local expected="$2"; shift 2
    local actual
    actual=$("$@" 2>/dev/null | grep -v '^##yamref ')
    if ! diff -u "$expected" <(echo "$actual") >/dev/null 2>&1; then
        diff -u "$expected" <(echo "$actual") >&2 || true
        FAIL "$desc: output differs from golden $expected"
    fi
    PASS "$desc matches golden"
}

assert_line_count_eq() {
    local desc="$1"; local expected_n="$2"; shift 2
    local n
    n=$("$@" 2>/dev/null | grep -v '^#' | wc -l | tr -d ' ')
    if [ "$n" -ne "$expected_n" ]; then
        FAIL "$desc: expected $expected_n data lines, got $n"
    fi
    PASS "$desc has $expected_n data lines"
}

# Check that a FASTA output has the same letter histogram as input (useful for shuf tests)
assert_same_letters() {
    local desc="$1"
    local actual_fasta="$2"
    local expected_fasta="$3"
    local actual_letters expected_letters
    actual_letters=$(grep -v '^>' "$actual_fasta" | tr -d '\n' | fold -w1 | sort | uniq -c | sort -k2)
    expected_letters=$(grep -v '^>' "$expected_fasta" | tr -d '\n' | fold -w1 | sort | uniq -c | sort -k2)
    if [ "$actual_letters" != "$expected_letters" ]; then
        FAIL "$desc: letter histogram differs from input"
    fi
    PASS "$desc preserves letter histogram"
}
