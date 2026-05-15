#!/usr/bin/env python3
"""Score yamme runs against IUPAC-degenerate planted motifs.

For each (variant, dataset) pair, parse the yamme TSV and (where present)
the MEME PWM, then report:
  - which planted motifs were found (consensus matches IUPAC pattern)
  - PWM fidelity at degenerate positions (do the discovered probabilities
    actually reflect the degeneracy?)
"""
import csv
import os
import re
from collections import OrderedDict
from score import parse_yamme_tsv, load_planted, RES, TSV_DIR

IUPAC = {
    'A': set('A'), 'C': set('C'), 'G': set('G'), 'T': set('T'),
    'R': set('AG'), 'Y': set('CT'), 'S': set('CG'), 'W': set('AT'),
    'K': set('GT'), 'M': set('AC'),
    'B': set('CGT'), 'D': set('AGT'), 'H': set('ACT'), 'V': set('ACG'),
    'N': set('ACGT'),
}
RC_TABLE = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def iupac_match_pos(obs, exp):
    return obs in IUPAC.get(exp, set())


def revcomp_iupac(s):
    return s.translate(RC_TABLE)[::-1]


def iupac_align(candidate, pattern, min_overlap=None, max_dist=1):
    """Best alignment of `pattern` (IUPAC) against `candidate` (ACGT).
    Allows partial overlap (either string can overhang). Returns
    (best_dist, best_aligned_pos_range_in_candidate) or (None, None).
    """
    if min_overlap is None:
        # Require enough overlap that we're really matching the motif's core,
        # but allow short candidates against long IUPAC patterns (yamme often
        # converges on a shorter informative core of a long degenerate motif).
        min_overlap = max(5, min(len(pattern), len(candidate)) - 1)
    best = (None, None)
    pL, cL = len(pattern), len(candidate)
    for shift in range(-(pL - min_overlap), cL - min_overlap + 1):
        c_start = max(0, shift)
        p_start = max(0, -shift)
        overlap = min(cL - c_start, pL - p_start)
        if overlap < min_overlap:
            continue
        d = sum(
            0 if iupac_match_pos(candidate[c_start + i], pattern[p_start + i]) else 1
            for i in range(overlap)
        )
        if best[0] is None or d < best[0]:
            best = (d, (c_start, c_start + overlap, p_start, p_start + overlap))
            if d == 0:
                return best
    if best[0] is not None and best[0] <= max_dist:
        return best
    return (None, None)


def best_planted_match(candidate, pattern, min_overlap=None, max_dist=1):
    """Try forward and RC of the pattern, return (dist, alignment, strand)."""
    fwd = iupac_align(candidate, pattern, min_overlap, max_dist)
    rcp = revcomp_iupac(pattern)
    rcv = iupac_align(candidate, rcp, min_overlap, max_dist)
    if fwd[0] is not None and (rcv[0] is None or fwd[0] <= rcv[0]):
        return fwd[0], fwd[1], "+"
    if rcv[0] is not None:
        return rcv[0], rcv[1], "-"
    return None, None, None


def parse_meme_pwms(path):
    """Return {motif_name: {'consensus': str, 'matrix': [[a,c,g,t], ...]}}."""
    out = OrderedDict()
    if not os.path.exists(path):
        return out
    with open(path) as f:
        text = f.read()
    blocks = re.findall(
        r"MOTIF\s+(\S+)\s+(\S+).*?letter-probability matrix:[^\n]*\n((?:[ \t]*[0-9.eE+\-]+[ \t]+[0-9.eE+\-]+[ \t]+[0-9.eE+\-]+[ \t]+[0-9.eE+\-]+[ \t]*\n)+)",
        text, flags=re.DOTALL,
    )
    for name, cons, mat in blocks:
        rows = []
        for line in mat.strip().splitlines():
            parts = line.split()
            if len(parts) >= 4:
                rows.append([float(x) for x in parts[:4]])
        out[name] = {"consensus": cons, "matrix": rows}
    return out


def column_information_content(row):
    """Per-column IC in bits (background = uniform)."""
    import math
    bits = 2.0
    for p in row:
        if p > 0:
            bits += p * math.log2(p)
    return max(0.0, bits)


def degenerate_columns_match(matrix, c_range, p_range, pattern):
    """For an aligned planted-IUPAC region [p_range], check whether discovered
    PWM columns put significant probability on the right bases.

    Returns (n_pos, n_pos_correct) where pos_correct means: for an IUPAC code
    spanning k bases, the discovered PWM should sum >= 0.7 over those bases.
    """
    c0, c1, p0, p1 = c_range[0], c_range[1], p_range[0], p_range[1]
    n_correct = 0
    total = 0
    base_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(c1 - c0):
        if c0 + i >= len(matrix):
            break
        code = pattern[p0 + i]
        allowed = IUPAC[code]
        if len(allowed) == 4:
            continue  # 'N' — no constraint
        row = matrix[c0 + i]
        psum = sum(row[base_idx[b]] for b in allowed)
        total += 1
        if psum >= 0.7:
            n_correct += 1
    return total, n_correct


def main():
    planted = load_planted()
    variants = [
        ("default (k=20, r=2, MRH=10)",   "default"),
        ("recommended (k=4, r=2, MRH=20)","recommended"),
    ]
    datasets = ["D9", "D10", "D11"]

    summary = os.path.join(RES, "summary.md")
    with open(summary, "a") as f:
        f.write("\n## Degenerate motifs (D9–D11)\n\n")
        f.write("D9: TGASTCAG  (1 degenerate position, AP-1-like) @ 60% in 300×200bp.\n")
        f.write("D10: WGATAAR   (2 degenerate positions, GATA-like) @ 55% in 300×200bp.\n")
        f.write("D11: RGGNNYTYCC (5 degenerate positions, NF-kB-like) @ 60% in 400×200bp.\n\n")

        for vlabel, vkey in variants:
            f.write(f"### {vlabel}\n\n")
            for d in datasets:
                tsv = os.path.join(TSV_DIR, f"{d}_{vkey}.tsv")
                meme = os.path.join(TSV_DIR, f"{d}_{vkey}.meme")
                if not os.path.exists(tsv):
                    f.write(f"_{d}: results not found_\n\n")
                    continue
                motifs = parse_yamme_tsv(tsv)
                pwms = parse_meme_pwms(meme)
                truth_list = planted[d]
                f.write(f"**{d}** (motifs reported = {len(motifs)})\n\n")
                f.write("| planted (IUPAC) | rank | discovered consensus | dist | aligned | strand | PWM degen. cols correct | p-value |\n")
                f.write("|---|---|---|---|---|---|---|---|\n")
                for t in truth_list:
                    found_rank, found_cons, found_p = None, None, None
                    found_dist, found_align, found_strand = None, None, None
                    best_so_far = None
                    for m in motifs:
                        d_, align, strand = best_planted_match(m["consensus"], t)
                        if d_ is None:
                            continue
                        if best_so_far is None or d_ < best_so_far[0]:
                            best_so_far = (d_, m, align, strand)
                    if best_so_far is not None:
                        found_dist, m, align, strand = best_so_far
                        found_rank = m["rank"]
                        found_cons = m["consensus"]
                        found_p = m["pvalue"]
                        found_align = align
                        found_strand = strand
                    if found_rank is None:
                        f.write(f"| {t} | — | — | — | — | — | — | not found |\n")
                        continue
                    # PWM degenerate-column check
                    pwm_score = "n/a"
                    pwm_row = pwms.get(motifs[found_rank-1]["name"])
                    if pwm_row and found_align:
                        c0, c1, p0, p1 = found_align
                        pattern_to_use = t if found_strand == "+" else revcomp_iupac(t)
                        tot, corr = degenerate_columns_match(
                            pwm_row["matrix"], (c0, c1), (p0, p1), pattern_to_use
                        )
                        pwm_score = f"{corr}/{tot}" if tot > 0 else "—"
                    a_str = f"c[{found_align[0]}:{found_align[1]}]→p[{found_align[2]}:{found_align[3]}]" if found_align else "—"
                    f.write(f"| {t} | {found_rank} | {found_cons} | {found_dist} | {a_str} | {found_strand} | {pwm_score} | {found_p:.2e} |\n")
                f.write("\n")


if __name__ == "__main__":
    main()
