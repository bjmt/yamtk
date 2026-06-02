#!/usr/bin/env bash
#
# yamenr_bench.sh — micro-benchmark for yamtk enr.
#
# Generates a deterministic positives FASTA + multi-motif MEME (cached
# under /tmp), runs `yamtk enr` at several thread counts (shuffled
# negatives, deterministic seed), and reports wall clock + bp/s + peak
# MB. Throughput unit = bp/s of positives (negatives are an equal-size
# shuffle, so total work is ~2x bp/s).
#
# Usage:
#   bash bench/yamenr_bench.sh
#   bash bench/yamenr_bench.sh > baseline.tsv
#   bash bench/yamenr_bench.sh --baseline baseline.tsv
#
# Run from the repo root.

set -euo pipefail

YAMTK="${YAMTK:-./yamtk}"
BENCH_DIR="${BENCH_DIR:-/tmp/yamenr_bench}"
FA="$BENCH_DIR/pos.fa"
MM="$BENCH_DIR/motifs.meme"
TRIALS=3
JS=(1 2 4 8)
BASELINE=""

while [ $# -gt 0 ]; do
    case "$1" in
        --baseline) BASELINE="$2"; shift 2 ;;
        --trials)   TRIALS="$2";   shift 2 ;;
        --js)       IFS=',' read -ra JS <<<"$2"; shift 2 ;;
        -h|--help)
            sed -n '2,15p' "$0" >&2; exit 0 ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

if [ ! -x "$YAMTK" ]; then
    echo "error: yamtk binary not found at '$YAMTK' (set YAMTK=path or build first)." >&2
    exit 1
fi

mkdir -p "$BENCH_DIR"

if [ ! -s "$FA" ]; then
    echo "# generating $FA (one-time)..." >&2
    python3 - >"$FA" <<'PY'
import random
random.seed(0x5EAF00D)
b = 'ACGT'
# 10k seqs * 500 bp = 5 MB. Smaller than yamscan_bench because yamenr
# does two passes (pos + shuffled neg) plus per-motif stats, so an
# all-trials sweep needs to stay under ~3 minutes.
for i in range(10000):
    s = ''.join(random.choice(b) for _ in range(500))
    print(f'>s{i}'); print(s)
PY
fi
if [ ! -s "$MM" ]; then
    echo "# generating $MM (one-time)..." >&2
    cat >"$MM" <<'EOF'
MEME version 4
ALPHABET= ACGT
strands: + -
Background letter frequencies:
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF ebox CACGTG
letter-probability matrix: alength= 4 w= 6 nsites= 100 E= 0
 0.91 0.03 0.03 0.03
 0.03 0.91 0.03 0.03
 0.91 0.03 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.03 0.91
 0.03 0.03 0.91 0.03

MOTIF gcbox GGGCGG
letter-probability matrix: alength= 4 w= 6 nsites= 100 E= 0
 0.03 0.03 0.91 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.91 0.03
 0.03 0.91 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.91 0.03

MOTIF tata TATAAA
letter-probability matrix: alength= 4 w= 6 nsites= 100 E= 0
 0.03 0.03 0.03 0.91
 0.91 0.03 0.03 0.03
 0.03 0.03 0.03 0.91
 0.91 0.03 0.03 0.03
 0.91 0.03 0.03 0.03
 0.91 0.03 0.03 0.03

MOTIF wide ACGTACGTACGTAC
letter-probability matrix: alength= 4 w= 14 nsites= 100 E= 0
 0.91 0.03 0.03 0.03
 0.03 0.91 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.03 0.91
 0.91 0.03 0.03 0.03
 0.03 0.91 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.03 0.91
 0.91 0.03 0.03 0.03
 0.03 0.91 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.03 0.91
 0.91 0.03 0.03 0.03
 0.03 0.91 0.03 0.03

MOTIF medium CGTAACGT
letter-probability matrix: alength= 4 w= 8 nsites= 100 E= 0
 0.03 0.91 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.03 0.91
 0.91 0.03 0.03 0.03
 0.91 0.03 0.03 0.03
 0.03 0.91 0.03 0.03
 0.03 0.03 0.91 0.03
 0.03 0.03 0.03 0.91

EOF
fi

TOTAL_BP=$(awk '/^>/{next} {n+=length($0)} END{print n}' "$FA")
cat "$FA" >/dev/null  # warm cache

parse_peak_mb() {
    awk '
        /Approx. peak memory usage:/ {
            n=split($0, parts, ":"); val=parts[2]
            gsub(/,/, "", val); gsub(/\.$/, "", val)
            gsub(/^[ \t]+|[ \t]+$/, "", val)
            split(val, kv, " "); num=kv[1]+0; unit=kv[2]
            if (unit == "GB")      printf "%.2f", num * 1024
            else if (unit == "KB") printf "%.4f", num / 1024
            else                   printf "%.2f", num
            exit
        }' "$1"
}

STDERR_LOG="$BENCH_DIR/last.stderr"
echo -e "j\ttrial\twall_s\tbp_per_s\tpeak_mb"
for j in "${JS[@]}"; do
    for ((t=1; t<=TRIALS; t++)); do
        start=$(python3 -c 'import time; print(time.time())')
        # Use shuffled negatives with a fixed seed for determinism.
        "$YAMTK" enr -i "$FA" -m "$MM" -k 2 -s 42 -j "$j" -v \
                     -o /dev/null 2>"$STDERR_LOG"
        end=$(python3 -c 'import time; print(time.time())')
        wall=$(python3 -c "print(f'{$end - $start:.4f}')")
        bps=$(python3 -c "print(int($TOTAL_BP / ($end - $start)))")
        peak=$(parse_peak_mb "$STDERR_LOG")
        [ -z "$peak" ] && peak="-"
        echo -e "$j\t$t\t$wall\t$bps\t$peak"
    done
done | tee "$BENCH_DIR/last.tsv"

echo -e "\n# Median per -j (positives bp/s, peak MB):" >&2
python3 - "$BENCH_DIR/last.tsv" <<'PY' >&2
import sys, statistics, collections
d_bps = collections.defaultdict(list); d_mem = collections.defaultdict(list)
with open(sys.argv[1]) as f:
    for ln in f:
        if ln.startswith('#'): continue
        parts = ln.strip().split('\t')
        if len(parts) < 5 or parts[0] == 'j': continue
        j, _, _, bps, mem = parts
        d_bps[int(j)].append(float(bps))
        if mem != '-': d_mem[int(j)].append(float(mem))
for j in sorted(d_bps):
    mb = statistics.median(d_mem[j]) if d_mem.get(j) else float('nan')
    print(f"  j={j}  median {int(statistics.median(d_bps[j])):,} bp/s   {mb:6.1f} MB")
PY

if [ -n "$BASELINE" ]; then
    if [ ! -s "$BASELINE" ]; then
        echo "error: baseline '$BASELINE' missing or empty" >&2; exit 2
    fi
    echo -e "\n# Delta vs $BASELINE (median bp/s, median peak MB):" >&2
    python3 - "$BASELINE" "$BENCH_DIR/last.tsv" <<'PY' >&2
import sys, statistics, collections
def parse(path):
    bps = collections.defaultdict(list); mem = collections.defaultdict(list)
    with open(path) as f:
        for ln in f:
            if ln.startswith('#'): continue
            parts = ln.strip().split('\t')
            if len(parts) < 4 or parts[0] == 'j': continue
            j = int(parts[0]); bps[j].append(float(parts[3]))
            if len(parts) >= 5 and parts[4] != '-':
                mem[j].append(float(parts[4]))
    return ({j: statistics.median(v) for j, v in bps.items()},
            {j: statistics.median(v) for j, v in mem.items()})
a_bps, a_mem = parse(sys.argv[1])
b_bps, b_mem = parse(sys.argv[2])
for j in sorted(set(a_bps) | set(b_bps)):
    ba = a_bps.get(j); bn = b_bps.get(j)
    if not (ba and bn): continue
    dbps = (bn - ba) / ba * 100
    if j in a_mem and j in b_mem:
        ma, mn = a_mem[j], b_mem[j]
        dmem = (mn - ma) / ma * 100 if ma else 0.0
        print(f"  j={j}  {ba:>12,.0f} -> {bn:>12,.0f} bp/s  ({dbps:+.1f}%)   "
              f"{ma:6.1f} -> {mn:6.1f} MB  ({dmem:+.1f}%)")
    else:
        print(f"  j={j}  {ba:>12,.0f} -> {bn:>12,.0f} bp/s  ({dbps:+.1f}%)   (no mem data)")
PY
fi
