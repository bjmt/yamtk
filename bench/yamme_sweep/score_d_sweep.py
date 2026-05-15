#!/usr/bin/env python3
"""Score the -D sweep: for each (-D value, post-hoc q-filter) combination,
report planted recall, noise count, worst-case planted q-value, and
specifically how multi-motif datasets (D3, D7, D8) fare under aggressive -D.
"""
import csv
import os
from score import parse_yamme_tsv, RES
from score_degenerate import iupac_align, revcomp_iupac
from score_q_sweep import load_planted_map, DATASETS

D_DIR = os.path.join(RES, "d_sweep")
DS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
Q_FILTERS = [1e-3]
MULTI_MOTIF_DS = ["D3", "D7", "D8"]


def matches_any(consensus, planted_list):
    for p in planted_list:
        for query in (p, revcomp_iupac(p)):
            d, _ = iupac_align(consensus, query, max_dist=2)
            if d is not None:
                return True
    return False


def best_q_per_planted(motifs, planted_list, q_filter):
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
    p = os.path.join(RES, "timings_d_sweep.tsv")
    with open(p) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            t[(r["dataset"], r["d"])] = float(r["wall_s"])
    return t


def main():
    planted = load_planted_map()
    timings = load_timings()
    total_planted = sum(len(planted[d]) for d in DATASETS)
    multi_planted = sum(len(planted[d]) for d in MULTI_MOTIF_DS)

    summary = os.path.join(RES, "summary.md")
    with open(summary, "a") as f:
        f.write("\n## -D sweep (cross-width dedup overlap threshold)\n\n")
        f.write(f"Total planted: **{total_planted}** across 27 fixtures.  "
                f"Multi-motif fixtures (D3, D7, D8): **{multi_planted}** planted.\n")
        f.write("All runs at `-t 1e-3, -q 1`; q-filter applied at scoring time.\n\n")

        f.write("| -D | overall recall | multi-motif recall (D3/D7/D8) | "
                "noise (q≤1e-3) | max planted q | total time (s) |\n")
        f.write("|---|---|---|---|---|---|\n")

        rows = []
        for d_val in DS:
            d_str = f"{d_val}"
            n_recov = 0
            multi_recov = 0
            n_noise = 0
            planted_qs = []
            total_time = 0.0
            for ds in DATASETS:
                tsv = os.path.join(D_DIR, f"{ds}_D{d_str}.tsv")
                if not os.path.exists(tsv):
                    continue
                motifs = parse_yamme_tsv(tsv)
                best_qs = best_q_per_planted(motifs, planted[ds], 1e-3)
                for p, q in best_qs.items():
                    if q is not None:
                        n_recov += 1
                        planted_qs.append(q)
                        if ds in MULTI_MOTIF_DS:
                            multi_recov += 1
                # Noise = motifs passing q≤1e-3 that don't match any planted
                for m in motifs:
                    if m["qvalue"] > 1e-3:
                        continue
                    if not matches_any(m["consensus"], planted[ds]):
                        n_noise += 1
                total_time += timings.get((ds, d_str), 0.0)

            max_q = max(planted_qs) if planted_qs else None
            max_str = f"{max_q:.2e}" if max_q else "—"
            rows.append((d_val, n_recov, multi_recov, n_noise, max_q, total_time))
            f.write(f"| {d_val} | {n_recov}/{total_planted} | "
                    f"{multi_recov}/{multi_planted} | "
                    f"{n_noise} | {max_str} | {total_time:.1f} |\n")

        # Recommendation
        f.write("\n### Recommendation\n\n")
        full_recall = [r for r in rows
                       if r[1] == total_planted and r[2] == multi_planted]
        if full_recall:
            # Among full-recall rows, pick the one with smallest max planted q
            # (least BH inflation).  Tiebreak: smaller noise count.
            best = min(full_recall, key=lambda r: (r[4], r[3]))
            f.write(f"`-D {best[0]}` gives full {best[1]}/{total_planted} recall, "
                    f"multi-motif recall {best[2]}/{multi_planted}, "
                    f"noise={best[3]}, worst-case planted q={best[4]:.2e}, "
                    f"time={best[5]:.1f}s.\n")
        else:
            f.write("No `-D` value achieves full recall — needs investigation.\n")

    print(f"Appended -D sweep to {summary}")


if __name__ == "__main__":
    main()
