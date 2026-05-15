#!/usr/bin/env python3
"""Choose a sensible default for yamme's -q flag.

Reads results/q_sweep/*.tsv (yamme at current defaults, no q-filter),
classifies every output motif as 'planted' (matches a ground-truth IUPAC
pattern) or 'noise', then scans q-value thresholds to find the one that
keeps 100% planted recall with minimum noise survival.

The scoring is IUPAC-aware (handles both specific motifs like D1's TGACTCAG
and degenerate motifs like D11's RGGNNYTYCC).
"""
import csv
import math
import os
from collections import defaultdict
from score import parse_yamme_tsv, RES
from score_degenerate import iupac_align, revcomp_iupac

Q_DIR = os.path.join(RES, "q_sweep")

DATASETS = [
    "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8",
    "D9",  "D9_v2",  "D9_v3",  "D9_v4",  "D9_v5",
    "D10", "D10_v2", "D10_v3", "D10_v4", "D10_v5",
    "D11", "D11_v2", "D11_v3", "D11_v4", "D11_v5",
    "D12", "D13", "D14",
    "D15", "D16",
]


def load_planted_map():
    """Return {dataset: [planted_iupac_str, ...]}."""
    m = defaultdict(list)
    with open(os.path.join(RES, "..", "fixtures", "planted.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            m[r["dataset"]].append(r["motif"])
    return m


def matches_any(consensus, planted_list):
    """Does discovered consensus match any planted IUPAC pattern within
    Hamming-2 over the overlap?  Tries forward and RC of each planted."""
    for p in planted_list:
        for query in (p, revcomp_iupac(p)):
            d, _ = iupac_align(consensus, query, max_dist=2)
            if d is not None:
                return True
    return False


def main():
    planted = load_planted_map()

    # Per-(dataset, planted-motif) collect q-values of discovered motifs that
    # match it.  Each entry is a SET of q-values (so each planted has a "best"
    # = smallest q-value).
    planted_keys = []
    for d in DATASETS:
        for p in planted[d]:
            planted_keys.append((d, p))

    planted_best_q = {k: None for k in planted_keys}
    noise_motifs = []  # list of (dataset, consensus, qvalue)

    for d in DATASETS:
        tsv = os.path.join(Q_DIR, f"{d}.tsv")
        if not os.path.exists(tsv):
            print(f"WARN: missing {tsv}")
            continue
        motifs = parse_yamme_tsv(tsv)
        for m in motifs:
            matched_any = False
            for p in planted[d]:
                for query in (p, revcomp_iupac(p)):
                    dist, _ = iupac_align(m["consensus"], query, max_dist=2)
                    if dist is None:
                        continue
                    matched_any = True
                    key = (d, p)
                    if planted_best_q[key] is None or m["qvalue"] < planted_best_q[key]:
                        planted_best_q[key] = m["qvalue"]
                    break
            if not matched_any:
                noise_motifs.append((d, m["consensus"], m["qvalue"]))

    total_planted = len(planted_keys)
    total_noise = len(noise_motifs)
    print(f"Total planted (dataset, motif) pairs: {total_planted}")
    print(f"Total noise motifs across all runs:   {total_noise}")
    print()

    # Sweep thresholds
    print(f"{'q-thresh':>12}  {'planted_kept':>15}  {'noise_kept':>11}  {'noise_drop':>11}  {'total_out':>10}")
    print("-" * 70)
    thresholds = [1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5,
                  1e-6, 1e-10, 1e-20]
    for tau in thresholds:
        n_planted = sum(1 for q in planted_best_q.values() if q is not None and q <= tau)
        n_noise = sum(1 for (_, _, q) in noise_motifs if q <= tau)
        noise_dropped = total_noise - n_noise
        out = n_planted + n_noise  # rough; planted ≤ total_planted (multiple matches collapsed)
        flag = "  ← lossless" if n_planted == total_planted else ""
        print(f"{tau:>12.0e}  {n_planted:>6}/{total_planted:<7}  {n_noise:>6}/{total_noise:<5}  "
              f"{noise_dropped:>6} ({100*noise_dropped/total_noise:.0f}%)  "
              f"{n_planted + n_noise:>10}{flag}")

    # Find the strictest threshold that still keeps all planted motifs
    candidate = None
    for tau in sorted(thresholds):
        n_planted = sum(1 for q in planted_best_q.values() if q is not None and q <= tau)
        if n_planted == total_planted:
            candidate = tau
            break
    print()
    print(f"Strictest threshold keeping ALL planted motifs: {candidate}")

    # Where do planted vs noise q-values cluster?
    print()
    print("Distribution of best q-value per planted motif:")
    qs = sorted(q for q in planted_best_q.values() if q is not None)
    if qs:
        print(f"  min={qs[0]:.2e}  q1={qs[len(qs)//4]:.2e}  median={qs[len(qs)//2]:.2e}  "
              f"q3={qs[3*len(qs)//4]:.2e}  max={qs[-1]:.2e}")
        missing = [k for k, v in planted_best_q.items() if v is None]
        if missing:
            print(f"  ... but {len(missing)} planted motif(s) were NOT recovered at all: {missing}")

    print()
    print("Distribution of noise q-values:")
    ns = sorted(q for (_, _, q) in noise_motifs)
    if ns:
        print(f"  min={ns[0]:.2e}  q1={ns[len(ns)//4]:.2e}  median={ns[len(ns)//2]:.2e}  "
              f"q3={ns[3*len(ns)//4]:.2e}  max={ns[-1]:.2e}")

    # If there's a clean gap, recommend the lowest noise q-value as the
    # threshold (rounded up to a clean number).
    if qs and ns:
        worst_planted_q = max(qs)
        best_noise_q = min(ns)
        print()
        print(f"Worst planted-motif q-value: {worst_planted_q:.2e}")
        print(f"Best  noise-motif    q-value: {best_noise_q:.2e}")
        if worst_planted_q < best_noise_q:
            # Pick something between them, rounded to a clean log-decade boundary.
            mid = (math.log10(worst_planted_q) + math.log10(best_noise_q)) / 2
            print(f"Clean gap; midpoint (log10) = 10^{mid:.2f} = {10**mid:.2e}")


if __name__ == "__main__":
    main()
