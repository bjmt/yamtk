#!/usr/bin/env python3
"""Score multi-motif fixtures (D7, D8). For each variant report which
planted motifs were found, at what rank, with what p-value."""
import csv
import os
from score import parse_yamme_tsv, best_match, revcomp, load_planted, RES, TSV_DIR

DATASETS = ["D7", "D8"]
VARIANTS = [
    ("default (k=20, r=2, MRH=10)",  "default"),
    ("recommended (k=4, r=2, MRH=20)","recommended"),
    ("hedge (k=10, r=2, MRH=20)",     "hedge_k10"),
]


def main():
    planted = load_planted()

    out = os.path.join(RES, "summary.md")
    with open(out, "a") as f:
        f.write("\n## Multi-motif fixtures (D7, D8)\n\n")
        f.write("D7: 3 motifs at decreasing strength (70%, 40%, 20%) in 400 seqs × 200 bp.\n")
        f.write("D8: 4 motifs at similar moderate strength (40-50%) in 400 seqs × 250 bp.\n\n")

        for label, key in VARIANTS:
            f.write(f"### {label}\n\n")
            for d in DATASETS:
                tsv = os.path.join(TSV_DIR, f"{d}_{key}.tsv")
                motifs = parse_yamme_tsv(tsv)
                truth = planted[d]
                # For each planted motif, find earliest matching rank
                found = {}
                for m in motifs:
                    match, dist = best_match(m["consensus"], truth)
                    if match and match not in found:
                        found[match] = (m["rank"], m["consensus"], m["pvalue"], dist)
                recall = len(found) / len(truth)
                f.write(f"**{d}** (recall = {len(found)}/{len(truth)}, "
                        f"motifs reported = {len(motifs)})\n\n")
                f.write("| planted | rank | consensus | dist | p-value |\n|---|---|---|---|---|\n")
                for t in truth:
                    if t in found:
                        rk, cons, p, d_ = found[t]
                        f.write(f"| {t} | {rk} | {cons} | {d_} | {p:.3e} |\n")
                    else:
                        f.write(f"| {t} | — | — | — | — (not found) |\n")
                f.write("\n")


if __name__ == "__main__":
    main()
