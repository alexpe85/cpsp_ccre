"""
Stage 1: Pairwise LD resolution.

Computes r^2 and D' for all candidate-locus pairs within 2Mb on the same
chromosome, using the Ensembl REST API's /ld/ endpoint against 1000 Genomes
phase 3 reference panels. Collapses non-independent clusters to a single
representative variant (strongest reported GWAS p-value), and confirms
independence (r^2 < 0.01) or flags LD (r^2 above threshold) for pairs found
close together.

NOT covered by this script: UKBB-LD. No public API exposes UK
Biobank-derived pairwise LD for arbitrary variants. See the README.
The MTF2/FNBP1L cross-check reported in the paper (r^2=0.780) required the
authors' own UK Biobank RAP access and cannot be reproduced by this script.

Usage:
    python 01_ld_resolution.py

Output:
    data/ld_resolution_results.json  (raw pairwise LD results)
    Prints a summary of clusters found and the retained locus list.
"""

import json
import os
import time
from itertools import combinations

import requests

from config import CANDIDATE_LOCI, PANEL_BY_SOURCE, DATA_DIR

ENSEMBL_LD_URL = "https://rest.ensembl.org/ld/human/pairwise/{rsid_a}/{rsid_b}"
LD_WINDOW_BP = 2_000_000
REQUEST_DELAY_S = 0.5  # be polite to the free Ensembl REST endpoint


def panel_for(source: str) -> str | None:
    return PANEL_BY_SOURCE.get(source)


def query_ld(rsid_a: str, rsid_b: str, population: str) -> dict | None:
    """
    Query Ensembl REST for pairwise LD between two rsIDs in a given 1000
    Genomes population. Returns the raw JSON response (a list; Ensembl
    sometimes returns multiple windows/phasing results for a pair), or None
    on failure.
    """
    url = ENSEMBL_LD_URL.format(rsid_a=rsid_a, rsid_b=rsid_b)
    params = {"population_name": population}
    headers = {"Content-Type": "application/json"}
    try:
        resp = requests.get(url, params=params, headers=headers, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException as e:
        print(f"  [WARN] LD query failed for {rsid_a}/{rsid_b} ({population}): {e}")
        return None


def find_close_pairs(loci: list[dict], window_bp: int = LD_WINDOW_BP) -> list[tuple[dict, dict]]:
    """All same-chromosome pairs within window_bp of each other, excluding
    variants with unresolved positions (pos=None, i.e. not yet lifted over)."""
    resolvable = [l for l in loci if l.get("pos") is not None]
    pairs = []
    for a, b in combinations(resolvable, 2):
        if a["chr"] == b["chr"] and abs(a["pos"] - b["pos"]) <= window_bp:
            pairs.append((a, b))
    return pairs


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    missing_pos = [l["rsid"] for l in CANDIDATE_LOCI if l.get("pos") is None]
    if missing_pos:
        print(
            "[NOTICE] The following variants have no confirmed GRCh38 position in "
            "config.py and will be SKIPPED (their source paper reported GRCh37 "
            "coordinates that need independent liftover confirmation before use):"
        )
        for rsid in missing_pos:
            print(f"    {rsid}")
        print(
            "This affects LD checks involving these variants; confirm coordinates "
            "via Ensembl's /map/ endpoint (GRCh37->GRCh38) before treating results "
            "as final.\n"
        )

    pairs = find_close_pairs(CANDIDATE_LOCI)
    print(f"Found {len(pairs)} same-chromosome pairs within {LD_WINDOW_BP:,} bp.\n")

    results = []
    for a, b in pairs:
        panel_a = panel_for(a["source"])
        panel_b = panel_for(b["source"])
        # Use the more specific / smaller-ancestry panel if sources differ;
        # otherwise use the shared panel. Cross-source pairs should be
        # spot-checked manually — this default is a reasonable starting point,
        # not a substitute for judgment (see the UKB-2/GLR-2 cross-check in
        # the paper, which used 1000G EUR despite one locus being a UK
        # Biobank/own-GWAS signal, because that's what was queryable).
        panel = panel_a or panel_b
        if panel is None:
            print(f"  [SKIP] {a['rsid']}/{b['rsid']}: no panel defined for source(s) "
                  f"{a['source']}/{b['source']}")
            continue

        dist = abs(a["pos"] - b["pos"])
        print(f"Querying {a['rsid']} vs {b['rsid']} ({dist:,} bp apart, panel={panel})...")
        raw = query_ld(a["rsid"], b["rsid"], panel)
        time.sleep(REQUEST_DELAY_S)

        results.append({
            "rsid_a": a["rsid"], "rsid_b": b["rsid"],
            "chr": a["chr"], "distance_bp": dist,
            "panel": panel, "raw_response": raw,
        })

    out_path = os.path.join(DATA_DIR, "ld_resolution_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {len(results)} pairwise LD results to {out_path}")
    print(
        "\nNEXT STEP (manual): review ld_resolution_results.json, decide cluster "
        "collapses using strongest-p-value-representative rule, and confirm "
        "independence (r^2 < 0.01) for any pair not collapsed. This step is "
        "left manual deliberately — the collapse decision also depends on each "
        "locus's reported GWAS p-value, which lives in the source papers' own "
        "tables, not in this script."
    )


if __name__ == "__main__":
    main()
