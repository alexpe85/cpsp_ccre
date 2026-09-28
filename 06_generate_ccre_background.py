"""
Stage 6: Annotate the 500 background positions against the V4 registry.

Annotates a fixed set of 500 random autosomal positions
(data/background_positions_n500.txt; GRCh38, >= 5 Mb from chromosome ends,
ENCODE-blacklist-excluded, >= 1 Mb from any study locus) with the same indexed
nearest-cCRE lookup used for the study loci (02_ccre_annotation.py), so both
sets are annotated identically.

Usage:  python 06_generate_ccre_background.py
Input:  data/background_positions_n500.txt   (one chr:pos per line)
Output: data/background_ccre_n500_v4.json
"""

import json
import os

# Reuse Stage 2's download/parse/nearest-match logic directly rather than
# duplicating it, so background and real-loci annotation are guaranteed
# methodologically identical.
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "stage2", os.path.join(os.path.dirname(__file__), "02_ccre_annotation.py")
)
stage2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stage2)

from config import DATA_DIR


def main():
    pos_path = os.path.join(DATA_DIR, "background_positions_n500.txt")
    if not os.path.exists(pos_path):
        raise FileNotFoundError(f"{pos_path} not found.")
    with open(pos_path) as f:
        old_bg = [{"position": line.strip()} for line in f if line.strip()]
    print(f"Loaded {len(old_bg)} background positions.")

    print("Ensuring V4 cCRE class files are downloaded (reuses cache from Stage 2 if already present)...")
    class_paths = stage2.ensure_all_files_downloaded()
    ctcf_bound_accessions = stage2.load_ctcf_bound_accessions(class_paths["_CTCF_BOUND_COMPOSITE"])
    print(f"  {len(ctcf_bound_accessions):,} CTCF-bound accessions loaded.\n")

    print("Building per-chromosome indexes (one-time cost, enables fast lookups for all 500 positions)...")
    indexes = stage2.build_all_indexes(class_paths)
    print()

    results = []
    for i, entry in enumerate(old_bg, 1):
        pos_str = entry["position"]  # e.g. "chr11:6713774"
        chrom_part, pos_part = pos_str.split(":")
        chrom = chrom_part.replace("chr", "")
        pos = int(pos_part)

        if i % 50 == 0 or i == 1:
            print(f"  [{i}/{len(old_bg)}] {pos_str}...")

        result = stage2.annotate_locus_indexed(chrom, pos, indexes, ctcf_bound_accessions)
        results.append({
            "position": pos_str,
            "nearest_cCRE_accession": result["accession"] if result else None,
            "classification": result["class"] if result else None,
            "distance_bp": result["distance_bp"] if result else None,
            "ctcf_bound": result["ctcf_bound"] if result else None,
        })

    out_path = os.path.join(DATA_DIR, "background_ccre_n500_v4.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    n_found = sum(1 for r in results if r["classification"] is not None)
    print(f"\nWrote {len(results)} results to {out_path}")
    print(f"  {n_found}/{len(results)} positions had a cCRE within the window.")
    if n_found == 0:
        print("  [NOTICE] Zero positions found any cCRE — check the same access risks noted in "
              "02_ccre_annotation.py's docstring (auth/download issues) before trusting this file.")


if __name__ == "__main__":
    main()
