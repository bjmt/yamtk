#!/usr/bin/env bash
source "$TESTDIR/lib.sh"
assert_exit0 "scan -h"  "$YAMTK" scan  -h
assert_exit0 "shuf -h"  "$YAMTK" shuf  -h
assert_exit0 "dedup -h" "$YAMTK" dedup -h
assert_exit0 "enr -h"  "$YAMTK" enr  -h
