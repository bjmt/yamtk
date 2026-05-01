#!/usr/bin/env bash
# Bug #7: -log(p_value) with p<=0 produced NaN/+inf, breaking qsort transitivity.
# p=0 should be handled (or rejected); p<0 must be rejected with an error.
source "$TESTDIR/lib.sh"

# p_value < 0: must exit non-zero with a diagnostic
bad_input=$(mktemp /tmp/yamtk_dedup_XXXXXX.txt)
trap 'rm -f "$bad_input"' EXIT
{
  printf '##yamscan v2.0.0 [ synthetic ]\n'
  printf '##MotifCount=1 MotifSize=5 SeqCount=1 SeqSize=100 GC=50%% Ns=0 MaxPossibleHits=96\n'
  printf '##seq_name\tstart\tend\tstrand\tmotif\tpvalue\tscore\tscore_pct\tmatch\n'
  printf 's1\t1\t5\t+\tm1\t-0.5\t3.0\t50.0\tACGTA\n'
  printf 's1\t6\t10\t+\tm1\t0.05\t3.0\t50.0\tACGTA\n'
} > "$bad_input"

if "$YAMTK" dedup -i "$bad_input" >/dev/null 2>&1; then
    echo "not ok - negative p-value input: expected non-zero exit" >&2; exit 1
fi
echo "ok - negative p-value in dedup input is rejected with error"
