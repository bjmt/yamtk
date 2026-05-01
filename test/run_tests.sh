#!/usr/bin/env bash
# yamtk test driver.
# Usage: ./test/run_tests.sh [filter]
#   filter: optional substring to match against case names (e.g. "regression")
# Sets YAMTK and TESTDIR for all case scripts.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

export YAMTK="$REPO_DIR/yamtk"
export TESTDIR="$SCRIPT_DIR"

if [ ! -x "$YAMTK" ]; then
    echo "ERROR: binary not found at $YAMTK — run 'make yamtk' first" >&2
    exit 1
fi

FILTER="${1:-}"
PASS=0
FAIL=0
FAILED=()

# Collect all case scripts using find (compatible with bash 3.2 on macOS)
mapfile -t CASES < <(find "$SCRIPT_DIR/cases" -name '*.sh' | sort) 2>/dev/null \
    || CASES=($(find "$SCRIPT_DIR/cases" -name '*.sh' | sort))

for case_path in "${CASES[@]}"; do
    name="${case_path#$SCRIPT_DIR/cases/}"
    name="${name%.sh}"

    if [ -n "$FILTER" ] && [[ "$name" != *"$FILTER"* ]]; then
        continue
    fi

    tmpdir=$(mktemp -d)
    if bash "$case_path" >"$tmpdir/out" 2>&1; then
        echo "PASS  $name"
        PASS=$((PASS + 1))
    else
        echo "FAIL  $name"
        sed 's/^/      /' "$tmpdir/out" >&2
        FAIL=$((FAIL + 1))
        FAILED+=("$name")
    fi
    rm -rf "$tmpdir"
done

echo ""
echo "Results: $PASS passed, $FAIL failed"

if [ $FAIL -gt 0 ]; then
    echo ""
    echo "Failed cases:"
    for f in "${FAILED[@]}"; do
        echo "  $f"
    done
    exit 1
fi
