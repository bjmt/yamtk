#!/usr/bin/env bash
# Bug #17: kseq_t was never kseq_destroy'd in the low-mem-with-motifs path
# and the non-low-mem (threaded) path; it was leaked on every normal scan run.
# Fix: call kseq_destroy(kseq) at the end of the low-mem scan loop; the
# threaded path has load_seqs() destroy it internally.
# Test: run leak detection on Linux only (macOS ASan doesn't do leak detection).
source "$TESTDIR/lib.sh"

if [ "$(uname)" != "Linux" ]; then
    echo "ok - kseq leak test skipped (leak detection only available on Linux)"
    exit 0
fi

# Only meaningful under a LeakSanitizer build (make check-debug)
if ! echo "$ASAN_OPTIONS" | grep -q "detect_leaks"; then
    export ASAN_OPTIONS="${ASAN_OPTIONS:-}:detect_leaks=1"
fi

# low-mem path
assert_exit0 "scan low-mem path exits 0 (no leak abort)" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 1

# threaded path
assert_exit0 "scan threaded path exits 0 (no leak abort)" \
    "$YAMTK" scan -m "$TESTDIR/motif.meme" -s "$TESTDIR/dna.fa" -t 0.05 -j 2
