#!/usr/bin/env python3
"""Score the -t sweep: for each (-t value, q-filter threshold) combination,
report planted recall, count of surviving noise motifs, and q-value stats of
recovered planted motifs.
"""
import csv
import math
import os
from statistics import mean
from score import parse_yamme_tsv, RES
from score_degenerate import iupac_align, revcomp_iupac
from score_q_sweep import load_planted_map, DATASETS

T_DIR = os.path.join(RES, "t_sweep")
TS = [1.0, 0.5, 0.05, 0.01, 0.001, 0.0001]
Q_THRESHES = [1.0, 1e-3, 1e-4]  # apply post-hoc at scoring time


def matches_any(consensus, planted_list):
    for p in planted_list:
        for query in (p, revcomp_iupac(p)):
            d, _ = iupac_align(consensus, query, max_dist=2)
            if d is not None:
                return True
    return False


def best_q_per_planted(motifs, planted_list, q_filter):
    """For each planted motif, find smallest q among matching discoveries
    that also pass q_filter. Returns dict {planted: best_q or None}."""
    out = {p: None for p in planted_list}
    for m in motifs:
        if m["qvalue"] > q_filter:
            continue
        for p in planted_list:
            for query in (p, revcomp_iupac(p)):
                dist, _ = iupac_align(m["consensus"], query, max_dist=2)
                if dist is None:
                    continue
                if out[p] is None or m["qvalue"] < out[p]:
                    out[p] = m["qvalue"]
                break
    return out


def load_timings():
    t = {}
    p = os.path.join(RES, "timings_t_sweep.tsv")
    with open(p) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            t[(r["dataset"], r["t"])] = float(r["wall_s"])
    return t


def main():
    planted = load_planted_map()
    timings = load_timings()
    # Total planted across all datasets
    total_planted = sum(len(planted[d]) for d in DATASETS)

    summary = os.path.join(RES, "summary.md")
    with open(summary, "a") as f:
        f.write("\n## -t sweep (with -q 1e-3 default applied at scoring time)\n\n")
        f.write(f"Total planted (dataset, motif) pairs: **{total_planted}**.  "
                "27 fixtures.\n\n")

        f.write("| -t | q-filter | planted recall | noise kept | "
                "median planted q | max planted q | total time (s) |\n")
        f.write("|---|---|---|---|---|---|---|\n")

        # For each (t, q_filter), compute aggregate metrics
        # Track best across the rows for picking a default
        candidates = []  # (t, q_filter, n_recov, n_noise, max_q, total_t)
        for t_val in TS:
            t_str = f"{t_val}"
            for q_filter in Q_THRESHES:
                n_recov = 0
                n_noise = 0
                planted_qs = []
                total_time = 0.0
                missing = []
                for d in DATASETS:
                    tsv = os.path.join(T_DIR, f"{d}_t{t_str}.tsv")
                    if not os.path.exists(tsv):
                        continue
                    motifs = parse_yamme_tsv(tsv)
                    best_qs = best_q_per_planted(motifs, planted[d], q_filter)
                    for p, q in best_qs.items():
                        if q is not None:
                            n_recov += 1
                            planted_qs.append(q)
                        else:
                            missing.append((d, p))
                    # Count noise = surviving motifs (q <= q_filter) that match no planted
                    for m in motifs:
                        if m["qvalue"] > q_filter:
                            continue
                        if not matches_any(m["consensus"], planted[d]):
                            n_noise += 1
                    total_time += timings.get((d, t_str), 0.0)
                median_q = sorted(planted_qs)[len(planted_qs)//2] if planted_qs else None
                max_q = max(planted_qs) if planted_qs else None
                median_str = f"{median_q:.2e}" if median_q else "—"
                max_str = f"{max_q:.2e}" if max_q else "—"
                q_str = "1.0 (off)" if q_filter == 1.0 else f"{q_filter:.0e}"
                f.write(f"| {t_val} | {q_str} | {n_recov}/{total_planted} | "
                        f"{n_noise} | {median_str} | {max_str} | "
                        f"{total_time:.1f} |\n")
                candidates.append((t_val, q_filter, n_recov, n_noise, max_q, total_time, planted_qs))

        # Recommendation: best -t for the default -q 1e-3
        f.write("\n### Recommendation\n\n")
        # filter to candidates with q_filter == 1e-3 and full recall
        full_recall = [c for c in candidates if c[1] == 1e-3 and c[2] == total_planted]
        if full_recall:
            # Among full-recall rows, minimize max_q (lowest worst-case q-value
            # for real motifs).
            best = min(full_recall, key=lambda c: (c[4], c[5]))
            f.write(f"With `-q 1e-3` applied: `-t {best[0]}` keeps full {total_planted}/{total_planted} "
                    f"recall, noise = {best[3]}, worst-case planted q-value = {best[4]:.2e}, "
                    f"total time = {best[5]:.1f}s.\n")
        else:
            f.write("No `-t` value achieves full recall under `-q 1e-3` — needs investigation.\n")

    print(f"Appended -t sweep to {summary}")


if __name__ == "__main__":
    main()
