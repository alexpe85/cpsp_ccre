"""
Stage 8: LD-proxy regulatory overlap.

For each lead SNP, takes all biallelic SNVs within +/-500 kb in the
ancestry-matched 1000 Genomes panel (GRCh38 liftover), computes r^2 to the lead
SNP from phased haplotypes (compute_ld from Stage 1b), keeps proxies with
r^2 >= 0.8, and tests whether the lead SNP or any proxy lies inside a V4 cCRE
(variant position within element boundaries; not the nearest-within-50 kb
criterion of Stage 2).

Panels by source study: PANEL_BY_SOURCE in config.py (CEU, EUR or GBR).
Requires `tabix` on PATH.

Usage:
    python 08_ld_proxy_annotation.py --test GLR-4
    python 08_ld_proxy_annotation.py --all
Input:  data/table1_final.csv
Output: data/ld_proxy_results.json (written after every locus)
"""

import argparse
import csv
import json
import os
import subprocess
import time

import requests

from config import DATA_DIR, load_loci, PANEL_BY_SOURCE

# --- Reuse Stage 1b's VCF access + haplotype/r^2 machinery ---
import importlib.util
_spec1b = importlib.util.spec_from_file_location(
    "stage1b", os.path.join(os.path.dirname(__file__), "01b_ld_longrange_direct.py")
)
stage1b = importlib.util.module_from_spec(_spec1b)
_spec1b.loader.exec_module(stage1b)

# --- Reuse Stage 2's V4 cCRE index-building machinery ---
_spec2 = importlib.util.spec_from_file_location(
    "stage2", os.path.join(os.path.dirname(__file__), "02_ccre_annotation.py")
)
stage2 = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(stage2)

R2_THRESHOLD = 0.8
WINDOW_BP = 500_000

# Maps table1_final.csv's free-text Surgery/Source field to a PANEL_BY_SOURCE key.
SOURCE_TO_PANEL_KEY = [
    ("UK Biobank (own)", "own_GWAS"),
    ("Li 2025 (Br J Anaesth", "Li_BrJAnaesth"),
    ("Li 2025 (suggestive)", "Li_Anaesthesia"),
    ("Li 2025 (GWS)", "Li_Anaesthesia"),
    ("Warner", "Warner"),
    ("Parisien", "Parisien"),
]


def panel_for_source(source_text: str) -> str:
    for prefix, key in SOURCE_TO_PANEL_KEY:
        if source_text.startswith(prefix) or prefix in source_text:
            panel = PANEL_BY_SOURCE.get(key)
            if panel:
                return panel
    raise ValueError(f"No panel mapping found for source text: {source_text!r}")


def _run_tabix_with_retry(args: list[str], timeout: int, max_retries: int = 3):
    """Shared retry helper for tabix subprocess calls -- same transient-failure
    pattern already seen twice in this pipeline (VCF header fetch, panel file
    fetch) hardened here too, proactively, rather than waiting for a third
    crash report on yet another unretried network call in this same script."""
    last_result = None
    for attempt in range(1, max_retries + 1):
        result = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        if result.returncode == 0 and result.stdout.strip():
            return result
        last_result = result
        if attempt < max_retries:
            wait = 2 ** attempt
            print(f"    [RETRY {attempt}/{max_retries}] {' '.join(args[:2])}... returned empty/error "
                  f"(rc={result.returncode}); waiting {wait}s before retry...")
            time.sleep(wait)
    return last_result


def fetch_region_records(chrom: str, center_pos: int, window_bp: int):
    """
    Fetch all biallelic SNV records in [center_pos - window_bp, center_pos + window_bp]
    via tabix, reusing Stage 1b's contig-naming auto-detection. Returns a list of
    dicts with pos/rsid/raw fields -- actual genotype extraction happens per-sample-
    subset separately, since the sample list can differ by population.
    """
    vcf_url = stage1b.VCF_URL_TEMPLATE.format(chrom=chrom)
    list_result = _run_tabix_with_retry(["tabix", "-l", vcf_url], timeout=60)
    available_contigs = set(list_result.stdout.strip().split("\n"))
    contig = chrom if chrom in available_contigs else (f"chr{chrom}" if f"chr{chrom}" in available_contigs else None)
    if contig is None:
        raise RuntimeError(f"Could not find contig '{chrom}' or 'chr{chrom}' in {vcf_url}")

    region = f"{contig}:{max(1, center_pos - window_bp)}-{center_pos + window_bp}"
    result = _run_tabix_with_retry(["tabix", vcf_url, region], timeout=300)
    if result.returncode != 0:
        raise RuntimeError(f"tabix failed for {region} after retries: {result.stderr}")

    records = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        pos, rsid, ref, alt = int(fields[1]), fields[2], fields[3], fields[4]
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue  # biallelic SNVs only, same convention as Stage 1b
        records.append({"pos": pos, "rsid": rsid, "line_fields": fields})
    return vcf_url, records


def extract_haplotypes_for_samples(fields, sample_indices):
    """Same phased-genotype extraction logic as Stage 1b's get_haplotypes,
    applied to one already-fetched record's raw fields."""
    format_field = fields[8].split(":")
    gt_idx = format_field.index("GT")
    sample_fields = fields[9:]
    haplotypes = []
    for idx in sample_indices:
        gt_str = sample_fields[idx].split(":")[gt_idx]
        sep = "|" if "|" in gt_str else "/"
        alleles = gt_str.split(sep)
        haplotypes.extend(int(a) for a in alleles if a in ("0", "1"))
    return haplotypes


def find_ld_proxies(rsid: str, chrom: str, pos: int, panel: str, window_bp: int = WINDOW_BP,
                     r2_threshold: float = R2_THRESHOLD):
    """
    Fetch all SNVs in the window, compute r^2 between the lead SNP and each
    other SNP using Stage 1b's verified formula, return those >= r2_threshold.
    """
    population = panel.split(":")[-1]  # e.g. "1000GENOMES:phase_3:EUR" -> "EUR"
    samples = stage1b.fetch_population_samples(population)

    vcf_url, records = fetch_region_records(chrom, pos, window_bp)
    print(f"  {len(records)} biallelic SNVs in {window_bp*2/1000:.0f}kb window around {rsid}")

    all_samples = stage1b._get_vcf_sample_order(vcf_url)
    sample_indices = [all_samples.index(s) for s in samples if s in all_samples]

    lead_record = next((r for r in records if r["rsid"] == rsid or r["pos"] == pos), None)
    if lead_record is None:
        print(f"  [WARN] Lead SNP {rsid} not found in fetched records -- cannot compute proxies.")
        return []
    lead_hap = extract_haplotypes_for_samples(lead_record["line_fields"], sample_indices)

    proxies = []
    for rec in records:
        if rec["pos"] == lead_record["pos"]:
            continue
        other_hap = extract_haplotypes_for_samples(rec["line_fields"], sample_indices)
        if len(other_hap) != len(lead_hap) or len(lead_hap) == 0:
            continue
        ld = stage1b.compute_ld(lead_hap, other_hap)
        if ld["r2"] is not None and ld["r2"] >= r2_threshold:
            proxies.append({"rsid": rec["rsid"], "pos": rec["pos"], "r2": ld["r2"], "d_prime": ld["d_prime"]})

    proxies.sort(key=lambda p: -p["r2"])
    print(f"  {len(proxies)} proxies at r2>={r2_threshold}")
    return proxies


def check_ccre_overlap(chrom: str, pos: int, indexes: dict, ctcf_bound_accessions: set):
    """Direct overlap (not nearest-within-window) -- a variant either sits
    inside a cCRE's [start, end] or it doesn't."""
    target_chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom
    for cls, index in indexes.items():
        for mid, start, end, accession, cls_label in index.get(target_chrom, []):
            if start <= pos <= end:
                return {"overlaps": True, "class": cls_label, "accession": accession,
                         "ctcf_bound": accession in ctcf_bound_accessions}
    return {"overlaps": False}


def annotate_locus_with_proxies(locus_row: dict, indexes: dict, ctcf_bound_accessions: set):
    rsid = locus_row["rsID"]
    chrom = locus_row["Chr"]
    pos = int(locus_row["Pos(GRCh38)"])
    panel = panel_for_source(locus_row["Surgery/Source"])

    print(f"\n{locus_row['Locus']} / {rsid} (chr{chrom}:{pos:,}, panel={panel})...")
    proxies = find_ld_proxies(rsid, chrom, pos, panel)

    lead_overlap = check_ccre_overlap(chrom, pos, indexes, ctcf_bound_accessions)
    proxy_overlaps = []
    for p in proxies:
        ov = check_ccre_overlap(chrom, p["pos"], indexes, ctcf_bound_accessions)
        if ov["overlaps"]:
            proxy_overlaps.append({**p, **ov})

    any_overlap = lead_overlap["overlaps"] or bool(proxy_overlaps)
    print(f"  Lead SNP overlaps cCRE: {lead_overlap['overlaps']}. "
          f"{len(proxy_overlaps)}/{len(proxies)} proxies overlap a cCRE. "
          f"LD-expanded overlap: {any_overlap}")

    return {
        "locus": locus_row["Locus"], "rsid": rsid, "chr": chrom, "pos": pos, "panel": panel,
        "lead_snp_ccre_overlap": lead_overlap,
        "n_proxies": len(proxies),
        "proxies_overlapping_ccre": proxy_overlaps,
        "any_overlap_lead_or_proxy": any_overlap,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", metavar="LOCUS_ID", help="Run a single locus first, e.g. --test GLR-4")
    parser.add_argument("--all", action="store_true", help="Run all 45 loci")
    args = parser.parse_args()

    loci = load_loci()
    if args.test:
        loci = [l for l in loci if l["Locus"] == args.test]
        if not loci:
            parser.error(f"Locus '{args.test}' not found in table1_final.csv")
    elif not args.all:
        parser.error("Specify --test <LOCUS_ID> or --all")

    print("Downloading/loading V4 cCRE class files and building indexes "
          "(reuses cache from earlier stages if present)...")
    class_paths = stage2.ensure_all_files_downloaded()
    indexes = stage2.build_all_indexes(class_paths)
    ctcf_bound_accessions = stage2.load_ctcf_bound_accessions(class_paths["_CTCF_BOUND_COMPOSITE"])
    print()

    results = []
    out_path = os.path.join(DATA_DIR, "ld_proxy_results.json")
    for row in loci:
        result = annotate_locus_with_proxies(row, indexes, ctcf_bound_accessions)
        results.append(result)
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)

    print(f"\nWrote {len(results)} results to {out_path}")

    if len(results) > 1:
        lead_only = sum(1 for r in results if r["lead_snp_ccre_overlap"]["overlaps"])
        expanded = sum(1 for r in results if r["any_overlap_lead_or_proxy"])
        print(f"\nSummary: lead-SNP-only overlap {lead_only}/{len(results)}; "
              f"LD-expanded overlap {expanded}/{len(results)}")


if __name__ == "__main__":
    main()
