#!/usr/bin/env python3
"""Score the degenerate-preset verification sweep. Compares PWM fidelity on
D9–D11 across candidate "ambiguous mode" presets.
"""
import os
from score import parse_yamme_tsv, load_planted, RES, TSV_DIR
from score_degenerate import (
    parse_meme_pwms, best_planted_match, degenerate_columns_match,
    revcomp_iupac,
)

DATASETS = ["D9", "D10", "D11"]
VARIANTS = [
    ("recommended (k=4, r=2, MRH=20)",      "degen_baseline"),
    ("k=20, r=2, MRH=10  [old default]",    "degen_k20_mrh10"),
    ("k=20, r=2, MRH=5",                    "degen_k20_mrh5"),
    ("k=40, r=2, MRH=10",                   "degen_k40_mrh10"),
    ("k=40, r=5, MRH=5",                    "degen_k40_r5_mrh5"),
]


def main():
    planted = load_planted()

    summary = os.path.join(RES, "summary.md")
    with open(summary, "a") as f:
        f.write("\n## Degenerate-preset verification\n\n")
        f.write("Goal: confirm which preset gives best PWM fidelity on degenerate motifs (D9–D11).\n\n")
        f.write("| variant | D9 PWM | D10 PWM | D11 PWM | D9 p | D10 p | D11 p | total time |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for vlabel, vkey in VARIANTS:
            row = [vlabel]
            ps = {}
            pwm_fields = []
            total_time = 0.0
            for d in DATASETS:
                tsv = os.path.join(TSV_DIR, f"{d}_{vkey}.tsv")
                meme = os.path.join(TSV_DIR, f"{d}_{vkey}.meme")
                motifs = parse_yamme_tsv(tsv)
                pwms = parse_meme_pwms(meme)
                truth = planted[d][0]
                pwm_score = "—"
                p = "—"
                for m in motifs:
                    dist, align, strand = best_planted_match(m["consensus"], truth)
                    if dist is None:
                        continue
                    pwm = pwms.get(m["name"])
                    if pwm and align:
                        c0, c1, p0, p1 = align
                        pat = truth if strand == "+" else revcomp_iupac(truth)
                        tot, corr = degenerate_columns_match(
                            pwm["matrix"], (c0, c1), (p0, p1), pat
                        )
                        if tot > 0:
                            pwm_score = f"{corr}/{tot}"
                        p = f"{m['pvalue']:.1e}"
                    break
                row.append(pwm_score)
                ps[d] = p
            row.append(ps.get("D9", "—"))
            row.append(ps.get("D10", "—"))
            row.append(ps.get("D11", "—"))
            # Total time from timings file
            tt_path = os.path.join(RES, "timings_degen_verify.tsv")
            if os.path.exists(tt_path):
                import csv
                with open(tt_path) as ft:
                    rd = csv.DictReader(ft, delimiter="\t")
                    for r in rd:
                        if r["variant"] == vkey:
                            total_time += float(r["wall_s"])
            row.append(f"{total_time:.1f}s")
            f.write("| " + " | ".join(row) + " |\n")
    print(f"Appended degen-verify section to {summary}")


if __name__ == "__main__":
    main()
