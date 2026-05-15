#!/usr/bin/env python3
"""Score the yamme parameter sweep.

Reads planted-motif ground truth (fixtures/planted.tsv) and every yamme
output TSV under results/tsv/, computes per-(dataset, k, r) metrics,
then writes results/summary.md with:
  - a 4x4 (TOP_K_SEEDS x REFINE_PASSES) summary table
  - per-dataset detail
  - a recommendation derived from the selection rule in the plan
"""
import csv
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
RES = os.path.join(HERE, "results")
TSV_DIR = os.path.join(RES, "tsv")


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def trim_consensus_with_pwm(consensus, pwm_probs, modal_threshold=0.5):
    """Trim flank columns whose modal-base probability is below threshold.
    consensus and pwm_probs are aligned column-by-column.
    """
    n = len(consensus)
    L = 0
    while L < n and max(pwm_probs[L]) < modal_threshold:
        L += 1
    R = n
    while R > L and max(pwm_probs[R - 1]) < modal_threshold:
        R -= 1
    return consensus[L:R]


def best_match(candidate, planted_motifs, max_dist=1, min_overlap=6):
    """Find any planted motif that aligns to the candidate (forward or RC)
    with at least `min_overlap` matching positions and at most `max_dist`
    mismatches over the overlap. Handles both:
      - planted embeds in candidate  (candidate has flanking context)
      - candidate embeds in planted  (candidate is a substring of the truth)
    Returns (motif, distance_over_overlap) or (None, None).
    """
    best = (None, None)
    for m in planted_motifs:
        for query in (m, revcomp(m)):
            # Try every shift such that the overlap is >= min_overlap.
            # Negative shift = query overhangs on the left of candidate.
            for shift in range(-(len(query) - min_overlap), len(candidate) - min_overlap + 1):
                c_start = max(0, shift)
                q_start = max(0, -shift)
                overlap = min(len(candidate) - c_start, len(query) - q_start)
                if overlap < min_overlap:
                    continue
                d = hamming(candidate[c_start:c_start + overlap],
                            query[q_start:q_start + overlap])
                if best[1] is None or d < best[1]:
                    best = (m, d)
                    if d == 0:
                        return best
    if best[1] is not None and best[1] <= max_dist:
        return best
    return (None, None)


def load_planted():
    planted = defaultdict(list)
    with open(os.path.join(FIX, "planted.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            planted[r["dataset"]].append(r["motif"])
    return planted


def parse_yamme_tsv(path):
    """Return list of motif dicts: {rank, width, consensus, pvalue, ...}."""
    motifs = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            parts = line.split("\t")
            # cols: motif, rank, width, consensus, nsites, seqs_pos, seqs_neg,
            #       sites_pos, sites_neg, n_pos, n_neg, pvalue, qvalue
            motifs.append({
                "name": parts[0],
                "rank": int(parts[1]),
                "width": int(parts[2]),
                "consensus": parts[3],
                "nsites": int(parts[4]),
                "seqs_pos": int(parts[5]),
                "seqs_neg": int(parts[6]),
                "pvalue": float(parts[11]),
                "qvalue": float(parts[12]),
            })
    return motifs


def load_timings():
    t = {}
    path = os.path.join(RES, "timings.tsv")
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            key = (row["dataset"], int(row["top_k"]), int(row["refine"]))
            t[key] = float(row["wall_s"])
    return t


def score_one(motifs, planted_motifs):
    """Score one (dataset, k, r) run."""
    # Without pwm_probs in the TSV we can only trim by surrounding context.
    # Use a heuristic: try matching against the raw consensus; the matcher
    # already searches all embeddings.
    recovered_any = False
    best_planted_rank = math.inf
    planted_p = None
    matched_motifs = set()
    n_match = 0
    n_noise = 0
    for m in motifs:
        match, dist = best_match(m["consensus"], planted_motifs)
        if match is not None:
            n_match += 1
            matched_motifs.add(match)
            if m["rank"] < best_planted_rank:
                best_planted_rank = m["rank"]
                planted_p = m["pvalue"]
        else:
            n_noise += 1
    recovery_top1 = motifs and best_match(motifs[0]["consensus"], planted_motifs)[0] is not None
    recall = len(matched_motifs) / len(planted_motifs) if planted_motifs else 0
    return {
        "recovery_top1": bool(recovery_top1),
        "recall": recall,
        "best_planted_rank": best_planted_rank if math.isfinite(best_planted_rank) else None,
        "planted_log10p": -math.log10(planted_p) if planted_p and planted_p > 0 else None,
        "n_motifs": len(motifs),
        "n_noise": n_noise,
        "n_match": n_match,
    }


def main():
    planted = load_planted()
    timings = load_timings()

    KS = [4, 10, 20, 40]
    RS = [2, 5, 10, 20]
    DATASETS = ["D1", "D2", "D3", "D4", "D5", "D6"]

    # Score grid: scores[(k, r)] = {dataset: metrics}
    scores = defaultdict(dict)
    for k in KS:
        for r in RS:
            for d in DATASETS:
                tsv = os.path.join(TSV_DIR, f"{d}_k{k}_r{r}.tsv")
                if not os.path.exists(tsv):
                    continue
                motifs = parse_yamme_tsv(tsv)
                s = score_one(motifs, planted[d])
                s["runtime_s"] = timings.get((d, k, r), float("nan"))
                scores[(k, r)][d] = s

    # Aggregate per grid point
    agg = {}
    for (k, r), per_d in scores.items():
        total_noise = sum(s["n_noise"] for s in per_d.values())
        total_time  = sum(s["runtime_s"] for s in per_d.values())
        log_ps = [s["planted_log10p"] for s in per_d.values() if s["planted_log10p"] is not None]
        mean_log_p = sum(log_ps) / len(log_ps) if log_ps else 0.0
        feasible_strict = all(per_d[d]["recovery_top1"] for d in DATASETS if d in per_d)
        n_recov = sum(1 for d in DATASETS if d in per_d and per_d[d]["recovery_top1"])
        agg[(k, r)] = {
            "total_noise": total_noise,
            "total_time": total_time,
            "mean_log10p": mean_log_p,
            "feasible_strict": feasible_strict and len(per_d) == len(DATASETS),
            "n_recov_top1": n_recov,
        }

    # Rank: feasible first, then min(total_noise), then min(total_time), then max(mean_log10p)
    ranked = sorted(agg.items(), key=lambda kv: (
        not kv[1]["feasible_strict"],
        -kv[1]["n_recov_top1"],
        kv[1]["total_noise"],
        kv[1]["total_time"],
        -kv[1]["mean_log10p"],
    ))

    out = os.path.join(RES, "summary.md")
    with open(out, "w") as f:
        f.write("# yamme parameter sweep — summary\n\n")
        f.write(f"Datasets: {', '.join(DATASETS)}\n\n")
        f.write("## Recommendation\n\n")
        if ranked:
            (best_k, best_r), best_agg = ranked[0]
            f.write(f"**Best `(TOP_K_SEEDS, REFINE_PASSES) = ({best_k}, {best_r})`**  ")
            f.write(f"(top-1 recovery on {best_agg['n_recov_top1']}/{len(DATASETS)} datasets, ")
            f.write(f"total noise motifs = {best_agg['total_noise']}, ")
            f.write(f"total time = {best_agg['total_time']:.1f}s, ")
            f.write(f"mean -log10p = {best_agg['mean_log10p']:.1f})\n\n")

        f.write("## Summary table\n\n")
        f.write("Cells show: recov/6, noise, time(s), mean(-log10 p)\n\n")
        f.write("| TOP_K \\ REFINE | " + " | ".join(str(r) for r in RS) + " |\n")
        f.write("|---|" + "|".join(["---"] * len(RS)) + "|\n")
        for k in KS:
            cells = []
            for r in RS:
                a = agg.get((k, r))
                if a is None:
                    cells.append(" ")
                    continue
                cells.append(f"{a['n_recov_top1']}/6 noise={a['total_noise']} {a['total_time']:.1f}s p={a['mean_log10p']:.1f}")
            f.write(f"| **k={k}** | " + " | ".join(cells) + " |\n")

        f.write("\n## Per-dataset detail\n\n")
        for d in DATASETS:
            f.write(f"### {d}  (planted: {', '.join(planted[d])})\n\n")
            f.write("| k | r | top1? | recall | best_rank | -log10p | n_motifs | n_noise | runtime_s |\n")
            f.write("|---|---|---|---|---|---|---|---|---|\n")
            for k in KS:
                for r in RS:
                    s = scores.get((k, r), {}).get(d)
                    if s is None:
                        continue
                    rec = "Y" if s["recovery_top1"] else "."
                    rk  = "-" if s["best_planted_rank"] is None else s["best_planted_rank"]
                    lp  = "-" if s["planted_log10p"]    is None else f"{s['planted_log10p']:.1f}"
                    f.write(f"| {k} | {r} | {rec} | {s['recall']:.2f} | {rk} | {lp} | {s['n_motifs']} | {s['n_noise']} | {s['runtime_s']:.1f} |\n")
            f.write("\n")

    print(f"Wrote {out}")
    if ranked:
        (best_k, best_r), _ = ranked[0]
        print(f"Recommended: TOP_K_SEEDS={best_k}  REFINE_PASSES={best_r}")


if __name__ == "__main__":
    main()
