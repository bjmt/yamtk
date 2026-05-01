#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_exit_nonzero "unknown subcommand" "$YAMTK" bogus_subcommand
