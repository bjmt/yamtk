#!/usr/bin/env python3
"""Score the follow-up one-axis sweeps (HAMMING_MISMATCH, MIN_REFINE_HITS)
against the (k=4, r=2) baseline.  Appends a section to results/summary.md.
"""
import csv
import os
from score import (parse_yamme_tsv, score_one, load_planted, HERE, RES, TSV_DIR)

DATASETS = ["D1", "D2", "D3", "D4", "D5", "D6"]
VARIANTS = [
    ("baseline (k=4,r=2)",   "k4_r2"),       # we re-use main-sweep TSVs
    ("HAMMING_MISMATCH=2",   "hm2_k4_r2"),
    ("MIN_REFINE_HITS=5",    "mrh5_k4_r2"),
    ("MIN_REFINE_HITS=20",   "mrh20_k4_r2"),
]


def load_followup_timings():
    t = {}
    p = os.path.join(RES, "timings_followup.tsv")
    if not os.path.exists(p):
        return t
    with open(p) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            t[(row["dataset"], row["variant"])] = float(row["wall_s"])
    return t


def load_baseline_timings():
    t = {}
    with open(os.path.join(RES, "timings.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            if row["top_k"] == "4" and row["refine"] == "2":
                t[row["dataset"]] = float(row["wall_s"])
    return t


def main():
    planted = load_planted()
    fu_times = load_followup_timings()
    bl_times = load_baseline_timings()

    rows = []
    for label, key in VARIANTS:
        for d in DATASETS:
            if key == "k4_r2":
                tsv = os.path.join(TSV_DIR, f"{d}_k4_r2.tsv")
                wall = bl_times.get(d, float("nan"))
            else:
                tsv = os.path.join(TSV_DIR, f"{d}_{key}.tsv")
                wall = fu_times.get((d, key), float("nan"))
            if not os.path.exists(tsv):
                continue
            motifs = parse_yamme_tsv(tsv)
            s = score_one(motifs, planted[d])
            s["dataset"] = d
            s["variant"] = label
            s["runtime_s"] = wall
            rows.append(s)

    summary = os.path.join(RES, "summary.md")
    with open(summary, "a") as f:
        f.write("\n## Follow-up one-axis sweeps (at k=4, r=2)\n\n")
        for label, key in VARIANTS:
            v_rows = [r for r in rows if r["variant"] == label]
            n_recov = sum(1 for r in v_rows if r["recovery_top1"])
            tot_noise = sum(r["n_noise"] for r in v_rows)
            tot_time = sum(r["runtime_s"] for r in v_rows)
            avg_lp = ([r["planted_log10p"] for r in v_rows if r["planted_log10p"] is not None])
            mean_lp = sum(avg_lp)/len(avg_lp) if avg_lp else 0
            f.write(f"**{label}**: top-1 recovery {n_recov}/6, total noise {tot_noise}, "
                    f"total time {tot_time:.1f}s, mean -log10p {mean_lp:.1f}\n\n")
            f.write("| dataset | top1? | recall | best_rank | -log10p | n_noise | runtime_s |\n")
            f.write("|---|---|---|---|---|---|---|\n")
            for r in v_rows:
                rec = "Y" if r["recovery_top1"] else "."
                rk = "-" if r["best_planted_rank"] is None else r["best_planted_rank"]
                lp = "-" if r["planted_log10p"] is None else f"{r['planted_log10p']:.1f}"
                f.write(f"| {r['dataset']} | {rec} | {r['recall']:.2f} | {rk} | {lp} | "
                        f"{r['n_noise']} | {r['runtime_s']:.1f} |\n")
            f.write("\n")
    print(f"Appended follow-up section to {summary}")


if __name__ == "__main__":
    main()
