#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_exit0        "version" "$YAMTK" version
assert_stderr_contains "version string" "v2" "$YAMTK" version
