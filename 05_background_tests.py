"""
Stage 5: Background enrichment tests (Fisher's exact test, two-sided).

1. cCRE class distribution: 45 loci vs 500 random autosomal positions, one test
   per V4 class (8 tests) with Benjamini-Hochberg FDR across classes. CTCF-bound
   status is tested separately (single test, uncorrected).
2. Gene-level PheWAS Moderate-or-above rate: 45 loci vs 500 random
   protein-coding genes classified with the same rank-based rule (single test).

Prints results for transcription into the manuscript; writes no file. Runs
offline on the files in data/.

Usage:  python 05_background_tests.py
Input:  data/table1_final.csv
        data/background_ccre_n500_v4.json     (Stage 6)
        data/phewas_background_n500.json      (Stages 4, 9, 10)
"""

import json
import os

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from config import load_loci, DATA_DIR, PHEWAS_STRONG_RANK_CUTOFF


FOOTNOTE_MARKERS = ["*", "\u2020", "\u2021", "\u00a7", "#", "\u00b6", "\u2016"]


def strip_footnote_markers(value: str) -> str:
    """Table 1's PheWAS column carries footnote markers appended to the tier
    (e.g. 'Moderate**', 'None\u2021'). Strip them to get the bare tier name."""
    for marker in FOOTNOTE_MARKERS:
        value = value.split(marker)[0]
    return value.strip()


def ccre_class_enrichment():
    print("=" * 60)
    print("TEST 1: cCRE class distribution vs. random genomic background (V4)")
    print("=" * 60)

    loci = load_loci()
    real_classes = pd.Series([l["cCRE class"] for l in loci]).value_counts()
    real_total = len(loci)

    bg_path = os.path.join(DATA_DIR, "background_ccre_n500_v4.json")
    if not os.path.exists(bg_path):
        print(f"  [SKIP] {bg_path} not found — run 06_generate_ccre_background.py first.")
        return
    with open(bg_path) as f:
        bg_data = json.load(f)
    bg_classified = [d for d in bg_data if d.get("classification")]
    bg_total = len(bg_data)  # denominator includes positions with no cCRE within window
    bg_classes = pd.Series([d["classification"] for d in bg_classified]).value_counts()

    classes = sorted(set(real_classes.index) | set(bg_classes.index))
    pvals, rows = [], []
    for c in classes:
        r = int(real_classes.get(c, 0))
        b = int(bg_classes.get(c, 0))
        table = [[r, real_total - r], [b, bg_total - b]]
        odds, p = fisher_exact(table)
        pvals.append(p)
        rows.append((c, r, real_total, b, bg_total, odds, p))

    _, p_adj, _, _ = multipletests(pvals, method="fdr_bh")

    print(f"\n{'Class':<18}{'Real':<10}{'Background':<14}{'OR':<8}{'p':<10}{'p_FDR':<8}")
    for (c, r, rt, b, bt, odds, p), padj in zip(rows, p_adj):
        print(f"{c:<18}{r}/{rt:<7}{b}/{bt:<11}{odds:<8.3f}{p:<10.4f}{padj:<8.4f}")

    # CTCF-bound composite comparison (separate axis from class; see Methods)
    real_ctcf = sum(1 for l in loci if l["CTCF-bound"] == "Yes")
    bg_ctcf = sum(1 for d in bg_data if d.get("ctcf_bound") is True)
    odds_c, p_c = fisher_exact([[real_ctcf, real_total - real_ctcf], [bg_ctcf, bg_total - bg_ctcf]])
    print(f"\nCTCF-bound: real {real_ctcf}/{real_total} ({real_ctcf/real_total*100:.1f}%) vs "
          f"background {bg_ctcf}/{bg_total} ({bg_ctcf/bg_total*100:.1f}%)  OR={odds_c:.3f} p={p_c:.4f}")


def phewas_background_enrichment():
    print("\n" + "=" * 60)
    print("TEST 2: Gene-level PheWAS Moderate-or-above vs. random gene background")
    print("=" * 60)

    loci = load_loci()
    real_mod_plus = sum(1 for l in loci if strip_footnote_markers(l["PheWAS"]) in ("Strong", "Moderate"))
    real_total = len(loci)

    bg_path = os.path.join(DATA_DIR, "phewas_background_n500.json")
    if not os.path.exists(bg_path):
        print(f"  [SKIP] {bg_path} not found — run "
              f"a concatenation of the 04_phewas_query.py results for the 300- and 200-gene background sets (see README).")
        return
    with open(bg_path) as f:
        bg_data = json.load(f)
    bg_total = len(bg_data)
    # Stage 4's output nests the classification under result.tier (rank-based,
    # same rule as the real loci) — NOT a top-level pain_flag field (that was
    # the older, score-based file's schema; see module docstring).
    bg_mod_plus = sum(1 for g in bg_data
                       if g.get("result") and g["result"].get("tier") in ("Strong", "Moderate"))
    bg_unresolved = sum(1 for g in bg_data if g.get("result") is None)
    if bg_unresolved:
        print(f"  [NOTE] {bg_unresolved}/{bg_total} background genes failed to resolve "
              f"(e.g. symbol not found in Open Targets) — excluded from tier counts but "
              f"still counted in the denominator below unless you decide otherwise.")

    table = [[real_mod_plus, real_total - real_mod_plus], [bg_mod_plus, bg_total - bg_mod_plus]]
    odds, p = fisher_exact(table)

    print(f"\nReal:       {real_mod_plus}/{real_total} = {real_mod_plus/real_total*100:.1f}%")
    print(f"Background: {bg_mod_plus}/{bg_total} = {bg_mod_plus/bg_total*100:.1f}%")
    print(f"OR = {odds:.3f}, p = {p:.4f}")


if __name__ == "__main__":
    ccre_class_enrichment()
    phewas_background_enrichment()
