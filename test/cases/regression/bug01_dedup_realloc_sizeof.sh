#!/usr/bin/env bash
# Bug #1: sizeof(pointer) instead of sizeof(struct) in dedup realloc.
# A bucket with > 256 hits triggered heap corruption before the fix.
# Under ASan (make debug) the unfixed binary aborts; fixed binary exits 0.
source "$TESTDIR/lib.sh"

out=$("$YAMTK" dedup -i "$TESTDIR/fixtures/dedup_large_bucket.txt" 2>/dev/null)
count=$(echo "$out" | grep -c '^s1' || true)

if [ $? -ne 0 ] || [ "$count" -lt 1 ]; then
    echo "not ok - dedup failed or produced no output on large bucket (heap corruption?)" >&2
    exit 1
fi
echo "ok - dedup handles > 256 hits in one bucket without corruption"
