"""
Stage 9: Extend the PheWAS background gene set from 300 to 500.

Draws 200 NEW random protein-coding genes (excluding the 45 study genes
and the existing 300 background genes already used) so the PheWAS
background can be expanded to match the cCRE background's n=500, without
discarding the 300 genes already queried via 04_phewas_query.py.

Builds a real, verified protein-coding gene catalog by querying Ensembl's
overlap/region endpoint (the same endpoint already used in
08_ld_proxy_annotation.py's cCRE-overlap logic) across each autosome, then
samples from it -- rather than trusting a hardcoded or half-remembered
gene list.

Usage:
    python 09_generate_phewas_background_extension.py

Input:
    data/background_genes_n300.txt   (existing 300, excluded from the draw)
    study_genes_exclude.txt          (45 loci's assigned genes, excluded)
      -- if this file isn't present, pass genes directly with --exclude-file

Output:
    data/background_genes_extension_n200.txt   (200 new gene symbols)
    data/background_genes_n500.txt             (existing 300 + new 200, for
                                                  convenience, so
                                                  --genes-file can target
                                                  either the extension alone
                                                  or the full combined set)
"""

import argparse
import os
import random
import time

import requests

from config import DATA_DIR

ENSEMBL_OVERLAP_URL = "https://rest.ensembl.org/overlap/region/human/{region}"
# Autosomes only, consistent with the cCRE background positions.
AUTOSOMES = [str(i) for i in range(1, 23)]
# Approximate GRCh38 chromosome lengths (bp), used only to define the query
# region per chromosome -- doesn't need to be exact, just needs to cover
# the whole chromosome so no genes are missed.
CHROM_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468,
}


def fetch_protein_coding_genes_for_chrom(chrom: str, max_retries: int = 3) -> list[str]:
    """Fetch all protein-coding gene symbols on one chromosome. Ensembl's
    overlap/region has a maximum query span (5Mb per the documented
    /overlap/region limits for the 'gene' feature), so this chunks the
    chromosome into 5Mb windows rather than querying the whole thing at
    once."""
    length = CHROM_LENGTHS[chrom]
    window = 5_000_000
    symbols = []
    start = 1
    while start <= length:
        end = min(start + window - 1, length)
        region = f"{chrom}:{start}-{end}"
        url = ENSEMBL_OVERLAP_URL.format(region=region)
        for attempt in range(1, max_retries + 1):
            try:
                resp = requests.get(
                    url, params={"feature": "gene", "biotype": "protein_coding"},
                    headers={"Content-Type": "application/json"}, timeout=30,
                )
                resp.raise_for_status()
                genes = resp.json()
                symbols.extend(g.get("external_name") for g in genes
                                if g.get("biotype") == "protein_coding" and g.get("external_name"))
                break
            except requests.RequestException as e:
                if attempt < max_retries:
                    wait = 2 ** attempt
                    print(f"    [RETRY {attempt}/{max_retries}] chr{chrom}:{start}-{end} "
                          f"failed ({e}); waiting {wait}s...")
                    time.sleep(wait)
                else:
                    print(f"    [WARN] chr{chrom}:{start}-{end} failed after {max_retries} "
                          f"attempts, skipping this window: {e}")
        start = end + 1
        time.sleep(0.3)
    return symbols


def build_genome_wide_catalog() -> set[str]:
    all_symbols = set()
    for chrom in AUTOSOMES:
        print(f"  chr{chrom}...")
        symbols = fetch_protein_coding_genes_for_chrom(chrom)
        all_symbols.update(symbols)
        print(f"    {len(symbols)} genes on chr{chrom} (running total: {len(all_symbols)})")
    return all_symbols


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=200, help="Number of new genes to draw (default 200)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    exclude = set()
    existing_bg_path = os.path.join(DATA_DIR, "background_genes_n300.txt")
    if os.path.exists(existing_bg_path):
        with open(existing_bg_path) as f:
            exclude.update(line.strip() for line in f if line.strip())
    study_genes_path = os.path.join(DATA_DIR, "study_genes_exclude.txt")
    if os.path.exists(study_genes_path):
        with open(study_genes_path) as f:
            exclude.update(line.strip() for line in f if line.strip())
    print(f"Excluding {len(exclude)} genes (45 study loci + existing 300 background).\n")

    print("Building genome-wide protein-coding gene catalog via Ensembl "
          "(this will take a while -- ~22 chromosomes x several 5Mb windows each)...")
    catalog = build_genome_wide_catalog()
    print(f"\nTotal unique protein-coding genes found: {len(catalog)}")

    candidates = sorted(catalog - exclude)
    print(f"Candidates after exclusion: {len(candidates)}")

    random.seed(args.seed)
    drawn = random.sample(candidates, args.n)

    os.makedirs(DATA_DIR, exist_ok=True)
    ext_path = os.path.join(DATA_DIR, "background_genes_extension_n200.txt")
    with open(ext_path, "w") as f:
        for g in drawn:
            f.write(g + "\n")
    print(f"\nWrote {len(drawn)} new gene symbols to {ext_path}")

    if os.path.exists(existing_bg_path):
        combined_path = os.path.join(DATA_DIR, "background_genes_n500.txt")
        with open(existing_bg_path) as f:
            existing = [line.strip() for line in f if line.strip()]
        with open(combined_path, "w") as f:
            for g in existing + drawn:
                f.write(g + "\n")
        print(f"Wrote combined 500-gene list to {combined_path}")

    print(
        "\nNEXT STEP: run\n"
        f"  python 04_phewas_query.py --genes-file {os.path.basename(ext_path)}\n"
        "to query just the 200 new genes (faster than redoing all 500 -- the "
        "existing 300 already have results in phewas_background_n300.json). "
        "Then merge the two result files before re-running 05_background_tests.py."
    )


if __name__ == "__main__":
    main()
