#!/usr/bin/env python3
"""Score the full verification sweep:
 - yamme presets (default, proposed_A, hedge) on all degenerate fixtures
 - STREME on every fixture for external comparison
Aggregates across replicates by fixture class.
"""
import csv
import os
import re
import xml.etree.ElementTree as ET
from statistics import mean
from collections import defaultdict
from score import parse_yamme_tsv, load_planted, RES, TSV_DIR
from score_degenerate import (
    parse_meme_pwms, best_planted_match, degenerate_columns_match,
    revcomp_iupac, iupac_align, IUPAC, RC_TABLE,
)

STREME_DIR = os.path.join(RES, "streme")

# Replicate groups: rows in the aggregate table.
GROUPS = [
    # group label, planted motif, list of datasets in the group
    ("D9 mild-degen (TGASTCAG, 300×200bp, 60%)",
        "TGASTCAG", ["D9", "D9_v2", "D9_v3", "D9_v4", "D9_v5"]),
    ("D10 moderate-degen (WGATAAR, 300×200bp, 55%)",
        "WGATAAR", ["D10", "D10_v2", "D10_v3", "D10_v4", "D10_v5"]),
    ("D11 heavy-degen (RGGNNYTYCC, 400×200bp, 60%)",
        "RGGNNYTYCC", ["D11", "D11_v2", "D11_v3", "D11_v4", "D11_v5"]),
    ("D12 big-mild (TGASTCAG, 1500×200bp, 60%)",
        "TGASTCAG", ["D12"]),
    ("D13 big-moderate (WGATAAR, 1500×200bp, 55%)",
        "WGATAAR", ["D13"]),
    ("D14 big-heavy (RGGNNYTYCC, 2000×200bp, 60%)",
        "RGGNNYTYCC", ["D14"]),
    ("D15 stress-7IUPAC (WRYGATAARYW, 500×200bp, 60%)",
        "WRYGATAARYW", ["D15"]),
    ("D16 stress-IUPAC-core (ACGTSWKMACGT, 500×200bp, 55%)",
        "ACGTSWKMACGT", ["D16"]),
]

PRESETS = [
    ("default (k=4, r=2, MRH=20)",     "default"),
    ("proposed_A (k=20, r=2, MRH=10)", "proposed_A"),
    ("hedge (k=10, r=2, MRH=10)",      "hedge"),
]


def load_full_timings():
    t = {}
    p = os.path.join(RES, "timings_full_verify.tsv")
    if not os.path.exists(p):
        return t
    with open(p) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            t[(row["dataset"], row["variant"])] = float(row["wall_s"])
    return t


def score_yamme_run(dataset, preset_key, planted_iupac):
    """Score one yamme run. Returns dict with recovery, dist, pwm_frac, p, time."""
    tsv = os.path.join(TSV_DIR, f"{dataset}_{preset_key}.tsv")
    meme = os.path.join(TSV_DIR, f"{dataset}_{preset_key}.meme")
    if not os.path.exists(tsv):
        return None
    motifs = parse_yamme_tsv(tsv)
    pwms = parse_meme_pwms(meme)
    best = None
    for m in motifs:
        dist, align, strand = best_planted_match(m["consensus"], planted_iupac)
        if dist is None:
            continue
        if best is None or dist < best["dist"] or \
           (dist == best["dist"] and m["rank"] < best["rank"]):
            pwm_frac = None
            pwm = pwms.get(m["name"])
            if pwm and align:
                c0, c1, p0, p1 = align
                pat = planted_iupac if strand == "+" else revcomp_iupac(planted_iupac)
                tot, corr = degenerate_columns_match(
                    pwm["matrix"], (c0, c1), (p0, p1), pat
                )
                pwm_frac = corr / tot if tot > 0 else None
            best = {
                "rank": m["rank"],
                "consensus": m["consensus"],
                "dist": dist,
                "pwm_frac": pwm_frac,
                "pvalue": m["pvalue"],
            }
    return best


def score_streme_run(dataset, planted_iupac):
    """Score one STREME run. Returns dict similar to yamme."""
    xml_path = os.path.join(STREME_DIR, dataset, "streme.xml")
    if not os.path.exists(xml_path):
        return None
    tree = ET.parse(xml_path)
    root = tree.getroot()
    motifs_node = root.find("motifs")
    best = None
    for i, m in enumerate(motifs_node.findall("motif"), 1):
        name = m.get("id")
        consensus = name.split("-", 1)[1] if "-" in name else name
        # STREME consensus may contain IUPAC codes itself; treat it as such.
        # For matching, we compare two IUPAC strings: STREME consensus vs planted.
        # Strict: only equal codes match. Loose: STREME's degenerate code subsumes planted.
        # We'll be loose — if STREME's column "covers" the planted column's allowed set,
        # count as match.
        # Implementation: convert IUPAC→allowed-set per column, then count overlap.
        forward = streme_iupac_match(consensus, planted_iupac)
        rc = streme_iupac_match(consensus, revcomp_iupac(planted_iupac))
        m_dist, m_align, m_strand = max([forward, rc], key=lambda x: -x[0] if x[0] is not None else float("inf"))
        # We want the *smaller* distance. Re-pick:
        candidates = [c for c in [forward, rc] if c[0] is not None]
        if not candidates:
            continue
        candidates.sort(key=lambda x: x[0])
        dist, align, strand = candidates[0]
        if best is None or dist < best["dist"] or \
           (dist == best["dist"] and i < best["rank"]):
            # Parse STREME PWM for this motif
            pwm_matrix = []
            for pos in m.find("pos") if False else []:
                pass
            # STREME XML stores pwm under <pos>?? Easier to read text output.
            # Use streme.txt and find the matrix.
            pwm_matrix = parse_streme_matrix(dataset, name)
            pwm_frac = None
            if pwm_matrix and align:
                c0, c1, p0, p1 = align
                pat = planted_iupac if strand == "+" else revcomp_iupac(planted_iupac)
                tot, corr = degenerate_columns_match(
                    pwm_matrix, (c0, c1), (p0, p1), pat
                )
                pwm_frac = corr / tot if tot > 0 else None
            test_p = m.get("test_pvalue")
            best = {
                "rank": i,
                "consensus": consensus,
                "dist": dist,
                "pwm_frac": pwm_frac,
                "pvalue": float(test_p) if test_p else float("nan"),
                "strand": strand,
            }
    return best


def streme_iupac_match(candidate_iupac, planted_iupac, min_overlap=None):
    """Align candidate (possibly IUPAC) against planted (IUPAC), counting
    column mismatch as 0 if candidate's allowed set intersects planted's,
    1 otherwise."""
    if min_overlap is None:
        min_overlap = max(4, len(planted_iupac) - 4)
    pL, cL = len(planted_iupac), len(candidate_iupac)
    max_dist = 2  # allow some slack
    best = (None, None, None)
    for shift in range(-(pL - min_overlap), cL - min_overlap + 1):
        c_start = max(0, shift)
        p_start = max(0, -shift)
        overlap = min(cL - c_start, pL - p_start)
        if overlap < min_overlap:
            continue
        d = 0
        for i in range(overlap):
            ca = IUPAC.get(candidate_iupac[c_start + i].upper(), set())
            pa = IUPAC.get(planted_iupac[p_start + i].upper(), set())
            if not (ca & pa):
                d += 1
        if best[0] is None or d < best[0]:
            best = (d, (c_start, c_start + overlap, p_start, p_start + overlap), "+")
            if d == 0:
                return best
    if best[0] is not None and best[0] <= max_dist:
        return best
    return (None, None, None)


def parse_streme_matrix(dataset, motif_full_name):
    """Parse the PWM matrix for a given motif name from STREME text output."""
    txt = os.path.join(STREME_DIR, dataset, "streme.txt")
    if not os.path.exists(txt):
        return None
    with open(txt) as f:
        text = f.read()
    # Pattern: "MOTIF <full_name>" followed by lines, then matrix
    pat = re.compile(
        rf"MOTIF\s+{re.escape(motif_full_name)}\b.*?letter-probability matrix:[^\n]*\n"
        r"((?:[ \t]*[0-9.eE+\-]+[ \t]+[0-9.eE+\-]+[ \t]+[0-9.eE+\-]+[ \t]+[0-9.eE+\-]+[ \t]*\n)+)",
        re.DOTALL
    )
    m = pat.search(text)
    if not m:
        return None
    rows = []
    for line in m.group(1).strip().splitlines():
        parts = line.split()
        if len(parts) >= 4:
            rows.append([float(x) for x in parts[:4]])
    return rows


def main():
    planted = load_planted()
    timings = load_full_timings()

    summary = os.path.join(RES, "summary.md")
    with open(summary, "a") as f:
        f.write("\n## Full verification: yamme presets vs STREME on degenerate fixtures\n\n")
        f.write("**Setup**: 20 datasets — 15 replicates of D9–D11, 3 big variants D12–D14, "
                "2 stress tests D15–D16. yamme runs at 3 presets, STREME at default settings, "
                "all with seed=42 where applicable.\n\n")
        f.write("Per-row aggregation across replicates: `recall = fraction of replicates with rank-1 hit`, "
                "`mean PWM fidelity = mean(corr/tot)` across replicates that recovered, "
                "`time` = mean wall-clock seconds.\n\n")

        f.write("| group | tool/preset | recall@r=1 | mean PWM fidelity | mean -log10 p | mean time (s) |\n")
        f.write("|---|---|---|---|---|---|\n")

        for group_label, planted_iupac, datasets in GROUPS:
            # yamme presets
            for preset_label, preset_key in PRESETS:
                results = []
                for d in datasets:
                    r = score_yamme_run(d, preset_key, planted_iupac)
                    if r is None:
                        continue
                    r["time"] = timings.get((d, preset_key), float("nan"))
                    results.append(r)
                if not results:
                    continue
                recall = sum(1 for r in results if r["rank"] == 1) / len(datasets)
                pwm_fracs = [r["pwm_frac"] for r in results if r["pwm_frac"] is not None]
                mean_pwm = mean(pwm_fracs) if pwm_fracs else float("nan")
                logps = [-__import__("math").log10(r["pvalue"]) for r in results if r["pvalue"] > 0]
                mean_lp = mean(logps) if logps else float("nan")
                mean_t = mean(r["time"] for r in results if r["time"] == r["time"])
                f.write(f"| {group_label} | yamme {preset_label} | {recall:.2f} | "
                        f"{mean_pwm:.2f} | {mean_lp:.1f} | {mean_t:.2f} |\n")

            # STREME
            results = []
            for d in datasets:
                r = score_streme_run(d, planted_iupac)
                if r is None:
                    continue
                r["time"] = timings.get((d, "streme"), float("nan"))
                results.append(r)
            if results:
                recall = sum(1 for r in results if r["rank"] == 1) / len(datasets)
                pwm_fracs = [r["pwm_frac"] for r in results if r["pwm_frac"] is not None]
                mean_pwm = mean(pwm_fracs) if pwm_fracs else float("nan")
                logps = [-__import__("math").log10(r["pvalue"]) for r in results if r["pvalue"] > 0]
                mean_lp = mean(logps) if logps else float("nan")
                mean_t = mean(r["time"] for r in results if r["time"] == r["time"])
                f.write(f"| {group_label} | STREME | {recall:.2f} | "
                        f"{mean_pwm:.2f} | {mean_lp:.1f} | {mean_t:.2f} |\n")
            f.write("|  |  |  |  |  |  |\n")

        # Aggregate summary
        f.write("\n### Aggregate across all degenerate groups\n\n")
        f.write("| preset | mean recall | mean PWM fidelity | mean time (s) |\n")
        f.write("|---|---|---|---|\n")
        all_tools = [(p[1], p[0], "yamme") for p in PRESETS] + [("streme", "STREME", "streme")]
        for tool_key, tool_label, kind in all_tools:
            recalls = []
            pwms = []
            times = []
            for group_label, planted_iupac, datasets in GROUPS:
                for d in datasets:
                    if kind == "yamme":
                        r = score_yamme_run(d, tool_key, planted_iupac)
                        tt = timings.get((d, tool_key), float("nan"))
                    else:
                        r = score_streme_run(d, planted_iupac)
                        tt = timings.get((d, "streme"), float("nan"))
                    if r is None:
                        continue
                    recalls.append(1 if r["rank"] == 1 else 0)
                    if r["pwm_frac"] is not None:
                        pwms.append(r["pwm_frac"])
                    if tt == tt:
                        times.append(tt)
            mr = mean(recalls) if recalls else 0
            mp = mean(pwms) if pwms else float("nan")
            mt = mean(times) if times else float("nan")
            f.write(f"| {tool_label} | {mr:.2f} | {mp:.2f} | {mt:.2f} |\n")
    print(f"Appended full-verify section to {summary}")


if __name__ == "__main__":
    main()
