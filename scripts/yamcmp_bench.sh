#!/usr/bin/env bash
#
# yamcmp_bench.sh — micro-benchmark for yamtk cmp.
#
# Generates deterministic query and target motif databases (cached under
# /tmp), runs `yamtk cmp` at several thread counts, and reports wall
# clock + comparisons/s + peak MB. Throughput unit = (n_queries x
# n_targets) per wall-second.
#
# Usage:
#   bash scripts/yamcmp_bench.sh
#   bash scripts/yamcmp_bench.sh > baseline.tsv
#   bash scripts/yamcmp_bench.sh --baseline baseline.tsv

set -euo pipefail

YAMTK="${YAMTK:-./yamtk}"
BENCH_DIR="${BENCH_DIR:-/tmp/yamcmp_bench}"
QM="$BENCH_DIR/queries.meme"
TM="$BENCH_DIR/targets.meme"
TRIALS=3
JS=(1 2 4 8)
BASELINE=""
N_QUERIES=50
N_TARGETS=500

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

# Generate a MEME file of N random PWMs of varied widths. Concentration
# of each column is randomized to mimic real motif IC distributions.
gen_meme() {
    local out="$1" n="$2" seed="$3"
    python3 - "$out" "$n" "$seed" <<'PY'
import sys, random
out_path, n, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
random.seed(seed)
with open(out_path, 'w') as f:
    f.write('MEME version 4\nALPHABET= ACGT\nstrands: + -\n')
    f.write('Background letter frequencies:\nA 0.25 C 0.25 G 0.25 T 0.25\n\n')
    for i in range(n):
        w = random.choice([6, 7, 8, 9, 10, 12, 14])
        f.write(f'MOTIF m{i}_w{w}\n')
        f.write(f'letter-probability matrix: alength= 4 w= {w} nsites= 100 E= 0\n')
        for _ in range(w):
            # Pick a dominant base with high probability ~0.7-0.95, rest split.
            dom = random.randrange(4)
            dom_p = random.uniform(0.7, 0.95)
            rest = (1.0 - dom_p) / 3
            probs = [rest]*4; probs[dom] = dom_p
            f.write(' ' + ' '.join(f'{p:.4f}' for p in probs) + '\n')
        f.write('\n')
PY
}

if [ ! -s "$QM" ]; then
    echo "# generating $QM (${N_QUERIES} motifs, one-time)..." >&2
    gen_meme "$QM" "$N_QUERIES" 17
fi
if [ ! -s "$TM" ]; then
    echo "# generating $TM (${N_TARGETS} motifs, one-time)..." >&2
    gen_meme "$TM" "$N_TARGETS" 42
fi

TOTAL_CMP=$((N_QUERIES * N_TARGETS))
cat "$QM" "$TM" >/dev/null  # warm cache

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
echo -e "j\ttrial\twall_s\tcmp_per_s\tpeak_mb"
for j in "${JS[@]}"; do
    for ((t=1; t<=TRIALS; t++)); do
        start=$(python3 -c 'import time; print(time.time())')
        "$YAMTK" cmp -m "$QM" -t "$TM" -q 1.0 -j "$j" -v -o /dev/null 2>"$STDERR_LOG"
        end=$(python3 -c 'import time; print(time.time())')
        wall=$(python3 -c "print(f'{$end - $start:.4f}')")
        cps=$(python3 -c "print(int($TOTAL_CMP / ($end - $start)))")
        peak=$(parse_peak_mb "$STDERR_LOG")
        [ -z "$peak" ] && peak="-"
        echo -e "$j\t$t\t$wall\t$cps\t$peak"
    done
done | tee "$BENCH_DIR/last.tsv"

echo -e "\n# Median per -j (cmp/s, peak MB; cmp = ${N_QUERIES}q x ${N_TARGETS}t = ${TOTAL_CMP}):" >&2
python3 - "$BENCH_DIR/last.tsv" <<'PY' >&2
import sys, statistics, collections
d_bps = collections.defaultdict(list); d_mem = collections.defaultdict(list)
with open(sys.argv[1]) as f:
    for ln in f:
        if ln.startswith('#'): continue
        parts = ln.strip().split('\t')
        if len(parts) < 5 or parts[0] == 'j': continue
        j, _, _, cps, mem = parts
        d_bps[int(j)].append(float(cps))
        if mem != '-': d_mem[int(j)].append(float(mem))
for j in sorted(d_bps):
    mb = statistics.median(d_mem[j]) if d_mem.get(j) else float('nan')
    print(f"  j={j}  median {int(statistics.median(d_bps[j])):,} cmp/s   {mb:6.1f} MB")
PY

if [ -n "$BASELINE" ]; then
    if [ ! -s "$BASELINE" ]; then
        echo "error: baseline '$BASELINE' missing or empty" >&2; exit 2
    fi
    echo -e "\n# Delta vs $BASELINE (median cmp/s, median peak MB):" >&2
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
        print(f"  j={j}  {ba:>10,.0f} -> {bn:>10,.0f} cmp/s  ({dbps:+.1f}%)   "
              f"{ma:6.1f} -> {mn:6.1f} MB  ({dmem:+.1f}%)")
    else:
        print(f"  j={j}  {ba:>10,.0f} -> {bn:>10,.0f} cmp/s  ({dbps:+.1f}%)   (no mem data)")
PY
fi
