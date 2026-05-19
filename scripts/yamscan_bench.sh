#!/usr/bin/env bash
#
# yamscan_bench.sh — micro-benchmark for yamtk scan.
#
# Generates a deterministic synthetic FASTA + a small multi-motif MEME
# file (cached under /tmp), then runs `yamtk scan` at several thread
# counts and reports wall clock + bp/s. Optionally diffs against a
# prior run's TSV to surface per-thread deltas.
#
# Usage:
#   bash scripts/yamscan_bench.sh                 [ -- save current state ]
#   bash scripts/yamscan_bench.sh > baseline.tsv  [ -- capture baseline ]
#   bash scripts/yamscan_bench.sh --baseline baseline.tsv
#
# Run from the repo root.

set -euo pipefail

YAMTK="${YAMTK:-./yamtk}"
BENCH_DIR="${BENCH_DIR:-/tmp/yamscan_bench}"
FA="$BENCH_DIR/seqs.fa"
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
            sed -n '2,15p' "$0" >&2
            exit 0 ;;
        *) echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done

if [ ! -x "$YAMTK" ]; then
    echo "error: yamtk binary not found at '$YAMTK' (set YAMTK=path or build first)." >&2
    exit 1
fi

mkdir -p "$BENCH_DIR"

# ---- Generate fixture once, cache thereafter ----
if [ ! -s "$FA" ]; then
    echo "# generating $FA (one-time)..." >&2
    python3 - >"$FA" <<'PY'
import random
random.seed(0x5EAF00D)
b = 'ACGT'
# 50k seqs * 500 bp = 25 MB raw bases. Enough work for stable timing
# (~7-10 s/run at -j 1 on Apple silicon), but small enough that the
# full 12-trial sweep finishes in under 2 minutes.
for i in range(50000):
    s = ''.join(random.choice(b) for _ in range(500))
    print(f'>s{i}')
    print(s)
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

MOTIF wide  ACGTACGTACGTAC
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

# ---- Total bp for throughput ----
TOTAL_BP=$(awk '/^>/{next} {n+=length($0)} END{print n}' "$FA")

# ---- Warm up disk cache ----
cat "$FA" >/dev/null

# ---- Run ----
echo -e "j\ttrial\twall_s\tbp_per_s"
for j in "${JS[@]}"; do
    for ((t=1; t<=TRIALS; t++)); do
        start=$(python3 -c 'import time; print(time.time())')
        "$YAMTK" scan -m "$MM" -s "$FA" -j "$j" -o /dev/null 2>/dev/null
        end=$(python3 -c 'import time; print(time.time())')
        wall=$(python3 -c "print(f'{$end - $start:.4f}')")
        bps=$(python3 -c "print(int($TOTAL_BP / ($end - $start)))")
        echo -e "$j\t$t\t$wall\t$bps"
    done
done | tee "$BENCH_DIR/last.tsv"

# ---- Summary: median per -j ----
echo -e "\n# Median bp/s per -j:" >&2
awk -F'\t' '
    NR>1 && $1!="j" { vals[$1] = vals[$1] " " $4 }
    END {
        for (j in vals) {
            n=split(vals[j], a, " ");
            # a may have empty leading element; reindex
            m=0; for(i=1;i<=n;i++) if(a[i] != "") b[++m]=a[i];
            # sort b[]
            for(i=1;i<=m;i++) for(k=i+1;k<=m;k++) if(b[i]+0>b[k]+0) {tmp=b[i]; b[i]=b[k]; b[k]=tmp}
            mid=int((m+1)/2);
            printf "  j=%s  median %d bp/s\n", j, b[mid];
        }
    }' "$BENCH_DIR/last.tsv" | sort >&2

# ---- Compare against baseline ----
if [ -n "$BASELINE" ]; then
    if [ ! -s "$BASELINE" ]; then
        echo "error: baseline '$BASELINE' missing or empty" >&2; exit 2
    fi
    echo -e "\n# Delta vs $BASELINE (per j, median bp/s):" >&2
    python3 - "$BASELINE" "$BENCH_DIR/last.tsv" <<'PY' >&2
import sys, statistics, collections
def medians(path):
    d = collections.defaultdict(list)
    with open(path) as f:
        for ln in f:
            if ln.startswith('#'): continue
            parts = ln.strip().split('\t')
            if len(parts)<4 or parts[0]=='j': continue
            j, _, _, bps = parts
            d[int(j)].append(float(bps))
    return {j: statistics.median(v) for j, v in d.items()}
a = medians(sys.argv[1])
b = medians(sys.argv[2])
for j in sorted(set(a) | set(b)):
    base = a.get(j); new = b.get(j)
    if base and new:
        delta = (new - base) / base * 100
        print(f"  j={j}  {base:>12,.0f} -> {new:>12,.0f} bp/s  ({delta:+.1f}%)")
PY
fi
