#!/usr/bin/env python3
"""Generate synthetic FASTA fixtures for the yamme parameter sweep.

Writes 6 datasets to fixtures/D1_pos.fa .. D6_pos.fa, spanning the realistic
operating envelope (motif strength, GC skew, sequence length, multi-motif).
All datasets use deterministic per-fixture seeds, so re-running this script
reproduces the same files byte-for-byte.

Also writes fixtures/planted.tsv: a TSV listing the planted motif(s) per
fixture, consumed by score.py for ground-truth matching.
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
FIX_DIR = os.path.join(HERE, "fixtures")
os.makedirs(FIX_DIR, exist_ok=True)


def gen_seq(rng, length, gc=0.5):
    at_p = (1.0 - gc) / 2.0
    gc_p = gc / 2.0
    return rng.choices("ACGT", weights=[at_p, gc_p, gc_p, at_p], k=length)


IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
    'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'N': 'ACGT',
}


def sample_iupac(motif, rng):
    """Concrete instance of an IUPAC pattern by uniform sampling per position."""
    return ''.join(rng.choice(IUPAC[c]) for c in motif)


def plant(rng, seqs, motif, frac, margin=10):
    """Plant `motif` (which may contain IUPAC degenerate codes) in `frac` of seqs.
    Each insertion samples a concrete instance from the IUPAC pattern."""
    n_plant = int(len(seqs) * frac)
    chosen = rng.sample(range(len(seqs)), n_plant)
    L = len(motif)
    is_degen = any(c not in "ACGT" for c in motif)
    for i in chosen:
        seq = seqs[i]
        slen = len(seq)
        if slen < L + 2 * margin:
            pos = (slen - L) // 2
        else:
            pos = rng.randint(margin, slen - L - margin)
        inst = sample_iupac(motif, rng) if is_degen else motif
        seq[pos:pos + L] = list(inst)
    return n_plant


def write_fa(path, seqs, name_prefix="seq"):
    with open(path, "w") as f:
        for i, s in enumerate(seqs, 1):
            f.write(f">{name_prefix}_{i}\n{''.join(s)}\n")


def build(seed, n_seqs, seq_len, gc=0.5):
    rng = random.Random(seed)
    return rng, [gen_seq(rng, seq_len, gc=gc) for _ in range(n_seqs)]


def main():
    planted = []

    # D1 — easy: 300x200bp, TGACTCAG in 75%
    rng, seqs = build(seed=1001, n_seqs=300, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "TGACTCAG", 0.75)
    write_fa(os.path.join(FIX_DIR, "D1_pos.fa"), seqs)
    planted.append(("D1", "TGACTCAG", n, len(seqs)))

    # D2 — weak: 300x200bp, TGACTCAG in 30%
    rng, seqs = build(seed=1002, n_seqs=300, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "TGACTCAG", 0.30)
    write_fa(os.path.join(FIX_DIR, "D2_pos.fa"), seqs)
    planted.append(("D2", "TGACTCAG", n, len(seqs)))

    # D3 — two motifs: CACGTGCACG in 50%, GCAAACCC in 40% (independent plantings)
    rng, seqs = build(seed=1003, n_seqs=300, seq_len=200, gc=0.5)
    n1 = plant(rng, seqs, "CACGTGCACG", 0.50)
    n2 = plant(rng, seqs, "GCAAACCC", 0.40)
    write_fa(os.path.join(FIX_DIR, "D3_pos.fa"), seqs)
    planted.append(("D3", "CACGTGCACG", n1, len(seqs)))
    planted.append(("D3", "GCAAACCC", n2, len(seqs)))

    # D4 — GC-skewed bg: 300x200bp, TGACTCAG in 60%, GC=0.70
    rng, seqs = build(seed=1004, n_seqs=300, seq_len=200, gc=0.70)
    n = plant(rng, seqs, "TGACTCAG", 0.60)
    write_fa(os.path.join(FIX_DIR, "D4_pos.fa"), seqs)
    planted.append(("D4", "TGACTCAG", n, len(seqs)))

    # D5 — short seqs: 600x60bp, TGACTCAG in 70%
    rng, seqs = build(seed=1005, n_seqs=600, seq_len=60, gc=0.5)
    n = plant(rng, seqs, "TGACTCAG", 0.70, margin=5)
    write_fa(os.path.join(FIX_DIR, "D5_pos.fa"), seqs)
    planted.append(("D5", "TGACTCAG", n, len(seqs)))

    # D6 — long seqs: 100x1000bp, CACGTGCACG in 80%
    rng, seqs = build(seed=1006, n_seqs=100, seq_len=1000, gc=0.5)
    n = plant(rng, seqs, "CACGTGCACG", 0.80)
    write_fa(os.path.join(FIX_DIR, "D6_pos.fa"), seqs)
    planted.append(("D6", "CACGTGCACG", n, len(seqs)))

    # D7 — three motifs of decreasing strength (recall test):
    #   strong CACGTGCACG @ 70%, medium GCAAACCC @ 40%, weak TTGACGTCA @ 20%
    rng, seqs = build(seed=1007, n_seqs=400, seq_len=200, gc=0.5)
    n1 = plant(rng, seqs, "CACGTGCACG", 0.70)
    n2 = plant(rng, seqs, "GCAAACCC",   0.40)
    n3 = plant(rng, seqs, "TTGACGTCA",  0.20)
    write_fa(os.path.join(FIX_DIR, "D7_pos.fa"), seqs)
    planted.append(("D7", "CACGTGCACG", n1, len(seqs)))
    planted.append(("D7", "GCAAACCC",   n2, len(seqs)))
    planted.append(("D7", "TTGACGTCA",  n3, len(seqs)))

    # D9 — single mildly degenerate motif: TGASTCAG (S = C/G, AP-1 style) @ 60%
    rng, seqs = build(seed=1009, n_seqs=300, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "TGASTCAG", 0.60)
    write_fa(os.path.join(FIX_DIR, "D9_pos.fa"), seqs)
    planted.append(("D9", "TGASTCAG", n, len(seqs)))

    # D10 — moderately degenerate motif: WGATAAR (W=A/T, R=A/G, GATA-like) @ 55%
    rng, seqs = build(seed=1010, n_seqs=300, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "WGATAAR", 0.55)
    write_fa(os.path.join(FIX_DIR, "D10_pos.fa"), seqs)
    planted.append(("D10", "WGATAAR", n, len(seqs)))

    # D11 — heavily degenerate motif: RGGNNYTYCC (NF-kB-like) @ 60%
    rng, seqs = build(seed=1011, n_seqs=400, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "RGGNNYTYCC", 0.60)
    write_fa(os.path.join(FIX_DIR, "D11_pos.fa"), seqs)
    planted.append(("D11", "RGGNNYTYCC", n, len(seqs)))

    # ---- Replicates of degenerate fixtures (different RNG seeds) ----
    # Pattern: <D9|D10|D11>_v<2..5>  with seeds 2000+idx
    # Replicate config matches the D9/D10/D11 base config.
    for v in range(2, 6):
        rng, seqs = build(seed=2000 + (v - 2) * 10 + 9, n_seqs=300, seq_len=200, gc=0.5)
        n = plant(rng, seqs, "TGASTCAG", 0.60)
        path = os.path.join(FIX_DIR, f"D9_v{v}_pos.fa")
        write_fa(path, seqs)
        planted.append((f"D9_v{v}", "TGASTCAG", n, len(seqs)))

        rng, seqs = build(seed=2000 + (v - 2) * 10 + 10, n_seqs=300, seq_len=200, gc=0.5)
        n = plant(rng, seqs, "WGATAAR", 0.55)
        write_fa(os.path.join(FIX_DIR, f"D10_v{v}_pos.fa"), seqs)
        planted.append((f"D10_v{v}", "WGATAAR", n, len(seqs)))

        rng, seqs = build(seed=2000 + (v - 2) * 10 + 11, n_seqs=400, seq_len=200, gc=0.5)
        n = plant(rng, seqs, "RGGNNYTYCC", 0.60)
        write_fa(os.path.join(FIX_DIR, f"D11_v{v}_pos.fa"), seqs)
        planted.append((f"D11_v{v}", "RGGNNYTYCC", n, len(seqs)))

    # ---- Bigger sequence sets (5x larger than D9-D11) ----
    rng, seqs = build(seed=3009, n_seqs=1500, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "TGASTCAG", 0.60)
    write_fa(os.path.join(FIX_DIR, "D12_pos.fa"), seqs)
    planted.append(("D12", "TGASTCAG", n, len(seqs)))

    rng, seqs = build(seed=3010, n_seqs=1500, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "WGATAAR", 0.55)
    write_fa(os.path.join(FIX_DIR, "D13_pos.fa"), seqs)
    planted.append(("D13", "WGATAAR", n, len(seqs)))

    rng, seqs = build(seed=3011, n_seqs=2000, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "RGGNNYTYCC", 0.60)
    write_fa(os.path.join(FIX_DIR, "D14_pos.fa"), seqs)
    planted.append(("D14", "RGGNNYTYCC", n, len(seqs)))

    # ---- Heavily degenerate, stress tests ----
    # D15: 7 IUPAC positions out of 10 (very degenerate)
    rng, seqs = build(seed=4015, n_seqs=500, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "WRYGATAARYW", 0.60)
    write_fa(os.path.join(FIX_DIR, "D15_pos.fa"), seqs)
    planted.append(("D15", "WRYGATAARYW", n, len(seqs)))

    # D16: longer motif with central IUPAC core
    rng, seqs = build(seed=4016, n_seqs=500, seq_len=200, gc=0.5)
    n = plant(rng, seqs, "ACGTSWKMACGT", 0.55)
    write_fa(os.path.join(FIX_DIR, "D16_pos.fa"), seqs)
    planted.append(("D16", "ACGTSWKMACGT", n, len(seqs)))

    # D8 — four motifs at similar moderate strength (40-50%):
    rng, seqs = build(seed=1008, n_seqs=400, seq_len=250, gc=0.5)
    n1 = plant(rng, seqs, "TGACTCAG",    0.50)
    n2 = plant(rng, seqs, "CACGTG",      0.45)
    n3 = plant(rng, seqs, "GATAAG",      0.45)
    n4 = plant(rng, seqs, "AACCGGTT",    0.40)
    write_fa(os.path.join(FIX_DIR, "D8_pos.fa"), seqs)
    planted.append(("D8", "TGACTCAG", n1, len(seqs)))
    planted.append(("D8", "CACGTG",   n2, len(seqs)))
    planted.append(("D8", "GATAAG",   n3, len(seqs)))
    planted.append(("D8", "AACCGGTT", n4, len(seqs)))

    with open(os.path.join(FIX_DIR, "planted.tsv"), "w") as f:
        f.write("dataset\tmotif\tn_planted\tn_total\n")
        for row in planted:
            f.write("\t".join(map(str, row)) + "\n")

    print(f"Wrote {len(planted)} planted-motif records to {FIX_DIR}/planted.tsv")
    for d, m, n, t in planted:
        print(f"  {d}: {m}  {n}/{t} ({100*n/t:.0f}%)")


if __name__ == "__main__":
    main()
