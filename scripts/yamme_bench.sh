#!/usr/bin/env bash
#
# yamme_bench.sh — micro-benchmark for yamtk me (motif discovery).
#
# Generates a deterministic positives FASTA with one implanted motif
# (cached under /tmp), runs `yamtk me` at several thread counts
# (deterministic shuffle seed), and reports wall clock + bp/s + peak
# MB. Throughput unit = positives bp/s.
#
# yamme is much heavier per bp than yamscan/yamenr (it enumerates
# seeds and refines them across widths), so the fixture is smaller and
# the parameter sweep is narrower.

set -euo pipefail

YAMTK="${YAMTK:-./yamtk}"
BENCH_DIR="${BENCH_DIR:-/tmp/yamme_bench}"
FA="$BENCH_DIR/pos.fa"
TRIALS=2
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
# 2000 seqs * 200 bp = 400 KB. Discovery is heavy: across widths 6-10,
# multiple refinement passes, full re-scans. Even this small fixture
# takes a few seconds per run.
implant = 'CACGTG'
for i in range(2000):
    s = [random.choice(b) for _ in range(200)]
    if i % 2 == 0:  # implant in half the sequences
        pos = random.randrange(0, 200 - len(implant))
        for k, c in enumerate(implant):
            s[pos+k] = c
    print(f'>s{i}'); print(''.join(s))
PY
fi

TOTAL_BP=$(awk '/^>/{next} {n+=length($0)} END{print n}' "$FA")
cat "$FA" >/dev/null

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
OUT_TSV="$BENCH_DIR/last.discovery.tsv"
OUT_MEME="$BENCH_DIR/last.discovery.meme"
echo -e "j\ttrial\twall_s\tbp_per_s\tpeak_mb"
for j in "${JS[@]}"; do
    for ((t=1; t<=TRIALS; t++)); do
        start=$(python3 -c 'import time; print(time.time())')
        # Fixed seed so shuffled negatives are deterministic; narrow
        # width range so the run finishes in reasonable wall time.
        "$YAMTK" me -i "$FA" -k 6 -K 8 -N 3 -s 42 -j "$j" -v \
                    -o "$OUT_TSV" -O "$OUT_MEME" 2>"$STDERR_LOG"
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
        print(f"  j={j}  {ba:>10,.0f} -> {bn:>10,.0f} bp/s  ({dbps:+.1f}%)   "
              f"{ma:6.1f} -> {mn:6.1f} MB  ({dmem:+.1f}%)")
    else:
        print(f"  j={j}  {ba:>10,.0f} -> {bn:>10,.0f} bp/s  ({dbps:+.1f}%)   (no mem data)")
PY
fi
