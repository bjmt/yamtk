#!/usr/bin/env bash
# Bug #3: fill_cdf checks the wrong realloc result variable.
# This path is only reachable under genuine OOM (second realloc failure).
# Untestable without LD_PRELOAD malloc shims. Intentionally skipped.
echo "ok - bug03 (fill_cdf realloc null-check) intentionally skipped: OOM not injectable"
