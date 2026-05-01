#!/usr/bin/env bash
# Bug #11: motif name parsers (MEME/JASPAR/HOMER) had no bound on the copy loop,
# allowing OOB writes into name[256]. Fix: clamp to MAX_NAME_SIZE-1.
# The scan should complete cleanly (truncated name is acceptable) or exit 1.
# The key requirement: no heap/stack corruption under ASan.
source "$TESTDIR/lib.sh"

# scan should either succeed (with truncated name) or exit 1 — but NOT crash
rc=0
"$YAMTK" scan -m "$TESTDIR/fixtures/longname_motif.meme" -s "$TESTDIR/dna.fa" \
    -t 0.9 -j 1 >/dev/null 2>&1 || rc=$?

# rc 0 (truncated name) or 1 (rejected) are both acceptable; a signal (rc>128) is not
if [ "$rc" -gt 1 ]; then
    echo "not ok - long motif name caused crash (exit $rc)" >&2; exit 1
fi
echo "ok - long motif name handled cleanly (exit $rc)"
