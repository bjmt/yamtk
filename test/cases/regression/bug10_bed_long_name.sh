#!/usr/bin/env bash
# Bug #10: BED names longer than MOTIF_VALUE_MAX_CHAR caused a stack overflow
# because the diagnostic was printed but execution continued into the stack buffer write.
# Fix: badexit after the diagnostic.
source "$TESTDIR/lib.sh"

if "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" \
        -x "$TESTDIR/fixtures/longname.bed" -j 1 >/dev/null 2>&1; then
    echo "not ok - long BED name: scan should exit non-zero" >&2; exit 1
fi
echo "ok - long BED name causes clean exit with error (no stack overflow)"
