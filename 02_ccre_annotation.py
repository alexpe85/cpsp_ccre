"""
Stage 2: ENCODE4 cCRE annotation (V4 registry; Moore et al. 2026, Nature).

For each locus in data/table1_final.csv, finds the nearest candidate
cis-regulatory element within +/-50 kb across all eight V4 classes (PLS, pELS,
dELS, CA-CTCF, CA-H3K4me3, CA-TF, CA, TF) and reports its class, accession and
distance to the element midpoint. CTCF-bound status is membership in SCREEN's
CTCF-bound cCRE composite file.

Data: SCREEN's static per-class BED files (GRCh38; 2,348,854 cCREs) and the
CTCF-bound composite file, downloaded once into data/ccre_v4/ (~130 MB,
git-ignored). SCREEN's public GraphQL API required authentication when this
analysis was run, so the bulk files are used. They serve the V4 registry only;
class labels are not comparable with annotations based on earlier registry
versions.

Format: parse_bed_line() takes the accession from column 5 and the class from
column 6 when present, otherwise from the file name. Run --peek to print the
first lines of each downloaded file and confirm the layout.

Usage:
    python 02_ccre_annotation.py --peek
    python 02_ccre_annotation.py --test 4:174713747      # single position
    python 02_ccre_annotation.py                         # all loci
Output: data/ccre_annotation_results.json
"""

import argparse
import json
import os

import requests

from config import load_loci, DATA_DIR

CCRE_V4_DIR = os.path.join(DATA_DIR, "ccre_v4")
BASE_URL = "https://downloads.wenglab.org/Registry-V4"
CLASS_FILES = {
    "PLS": f"{BASE_URL}/GRCh38-cCREs.PLS.bed",
    "pELS": f"{BASE_URL}/GRCh38-cCREs.pELS.bed",
    "dELS": f"{BASE_URL}/GRCh38-cCREs.dELS.bed",
    "CA-CTCF": f"{BASE_URL}/GRCh38-cCREs.CA-CTCF.bed",
    "CA-H3K4me3": f"{BASE_URL}/GRCh38-cCREs.CA-H3K4me3.bed",
    "CA-TF": f"{BASE_URL}/GRCh38-cCREs.CA-TF.bed",
    "CA": f"{BASE_URL}/GRCh38-cCREs.CA.bed",
    "TF": f"{BASE_URL}/GRCh38-cCREs.TF.bed",
}
CTCF_BOUND_URL = "https://downloads.wenglab.org/GRCh38-cCREs.CTCF-bound.bed"
CCRE_WINDOW_BP = 50_000


def download_if_needed(url: str, dest_path: str):
    if os.path.exists(dest_path):
        print(f"  (cached) {os.path.basename(dest_path)}")
        return
    print(f"  Downloading {url} -> {dest_path} ...")
    resp = requests.get(url, stream=True, timeout=300)
    resp.raise_for_status()
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1 << 20):
            f.write(chunk)
    print(f"  Done ({os.path.getsize(dest_path):,} bytes).")


def ensure_all_files_downloaded():
    os.makedirs(CCRE_V4_DIR, exist_ok=True)
    paths = {}
    for cls, url in CLASS_FILES.items():
        dest = os.path.join(CCRE_V4_DIR, os.path.basename(url))
        download_if_needed(url, dest)
        paths[cls] = dest
    ctcf_dest = os.path.join(CCRE_V4_DIR, os.path.basename(CTCF_BOUND_URL))
    download_if_needed(CTCF_BOUND_URL, ctcf_dest)
    paths["_CTCF_BOUND_COMPOSITE"] = ctcf_dest
    return paths


def peek_format(path: str, n: int = 3):
    """Print the first n lines of a downloaded file so its actual column
    structure can be visually confirmed before trusting parse_bed_line."""
    print(f"  First {n} lines of {os.path.basename(path)}:")
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            print(f"    {line.rstrip()}")


def parse_bed_line(line: str, fallback_class: str) -> dict | None:
    """
    Parse one BED line. Handles both the 6-column format confirmed from a
    comparable Weng-lab file (chrom, start, end, d_accession, e_accession,
    class_label) and a possible simpler per-class-file format lacking the
    class column, using fallback_class (the source filename's class) in
    that case. Returns None for malformed/short lines rather than crashing
    the whole run on one bad line.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 5:
        return None
    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
    e_accession = fields[4] if len(fields) >= 5 else fields[-1]
    label = fields[5] if len(fields) >= 6 else fallback_class
    return {"chrom": chrom, "start": start, "end": end,
            "accession": e_accession, "class": label.split(",")[0]}


def find_nearest_in_file(path: str, chrom: str, pos: int, window_bp: int,
                          fallback_class: str) -> dict | None:
    """Stream the file (do not load fully into memory — some class files
    are 80MB+) and keep only the nearest candidate within the window for
    this one position. O(file size) per call; fine for 45 loci x 8 files
    x one-time run, not fine for repeated interactive use — see README if
    this needs to scale to hundreds of positions (build an index instead)."""
    target_chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom
    best = None
    best_dist = None
    with open(path) as f:
        for line in f:
            if not line.startswith(target_chrom + "\t"):
                continue
            rec = parse_bed_line(line, fallback_class)
            if rec is None:
                continue
            mid = (rec["start"] + rec["end"]) / 2
            if abs(mid - pos) > window_bp:
                continue
            dist = abs(mid - pos)
            if best is None or dist < best_dist:
                best, best_dist = rec, dist
    if best is not None:
        best["distance_bp"] = round(best_dist)
    return best


def load_ctcf_bound_accessions(path: str) -> set[str]:
    accs = set()
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                accs.add(fields[4])
    return accs


def build_index(path: str, fallback_class: str) -> dict:
    """
    Load one class BED file fully into memory, once, as a per-chromosome
    list of (midpoint, start, end, accession, class) tuples sorted by
    midpoint. Enables O(log n) nearest-neighbor lookup via bisect instead
    of an O(file size) linear stream per query — necessary once query
    counts go beyond a handful (e.g. 500 background positions; see
    06_generate_ccre_background.py's runtime warning about the
    per-query-stream approach at that scale).
    """
    import bisect
    index: dict[str, list] = {}
    with open(path) as f:
        for line in f:
            rec = parse_bed_line(line, fallback_class)
            if rec is None:
                continue
            mid = (rec["start"] + rec["end"]) / 2
            index.setdefault(rec["chrom"], []).append((mid, rec["start"], rec["end"], rec["accession"], rec["class"]))
    for chrom in index:
        index[chrom].sort(key=lambda t: t[0])
    return index


def build_all_indexes(class_paths: dict) -> dict:
    """Build indexes for all 8 class files. One-time cost; reused across
    every query in a run, unlike find_nearest_in_file's per-call file scan."""
    import time
    indexes = {}
    for cls, path in class_paths.items():
        if cls == "_CTCF_BOUND_COMPOSITE":
            continue
        t0 = time.time()
        indexes[cls] = build_index(path, fallback_class=cls)
        print(f"  Indexed {cls}: {sum(len(v) for v in indexes[cls].values()):,} records "
              f"({time.time()-t0:.1f}s)")
    return indexes


def nearest_in_index(index: dict, chrom: str, pos: int, window_bp: int) -> dict | None:
    import bisect
    target_chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom
    chrom_list = index.get(target_chrom)
    if not chrom_list:
        return None
    midpoints = [t[0] for t in chrom_list]
    i = bisect.bisect_left(midpoints, pos)
    candidates = []
    for j in (i - 1, i):
        if 0 <= j < len(chrom_list):
            candidates.append(chrom_list[j])
    if not candidates:
        return None
    best = min(candidates, key=lambda t: abs(t[0] - pos))
    dist = abs(best[0] - pos)
    if dist > window_bp:
        return None
    return {"chrom": target_chrom, "start": best[1], "end": best[2],
            "accession": best[3], "class": best[4], "distance_bp": round(dist)}


def annotate_locus_indexed(chrom: str, pos: int, indexes: dict, ctcf_bound_accessions: set[str]) -> dict | None:
    """Indexed equivalent of annotate_locus() — same output shape, much
    faster for many queries since indexes are built once and reused."""
    candidates = []
    for cls, index in indexes.items():
        hit = nearest_in_index(index, chrom, pos, CCRE_WINDOW_BP)
        if hit is not None:
            candidates.append(hit)
    if not candidates:
        return None
    nearest = min(candidates, key=lambda r: r["distance_bp"])
    nearest["ctcf_bound"] = nearest["accession"] in ctcf_bound_accessions
    return nearest


def annotate_locus(chrom: str, pos: int, class_paths: dict, ctcf_bound_accessions: set[str]) -> dict | None:
    candidates = []
    for cls, path in class_paths.items():
        if cls == "_CTCF_BOUND_COMPOSITE":
            continue
        hit = find_nearest_in_file(path, chrom, pos, CCRE_WINDOW_BP, fallback_class=cls)
        if hit is not None:
            candidates.append(hit)
    if not candidates:
        return None
    nearest = min(candidates, key=lambda r: r["distance_bp"])
    nearest["ctcf_bound"] = nearest["accession"] in ctcf_bound_accessions
    return nearest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", metavar="CHR:POS",
                         help="Single position, e.g. --test 4:174713747 (GLRA3/rs17298280)")
    parser.add_argument("--peek", action="store_true",
                         help="Print first lines of each downloaded file (format check) and exit")
    args = parser.parse_args()

    print("Ensuring V4 cCRE class files are downloaded (one-time, ~130MB total)...")
    class_paths = ensure_all_files_downloaded()

    if args.peek:
        for cls, path in class_paths.items():
            peek_format(path)
        print("\nCompare the columns above against this script's docstring assumptions "
              "before trusting any annotation output.")
        return

    print("Loading CTCF-bound composite accession set...")
    ctcf_bound_accessions = load_ctcf_bound_accessions(class_paths["_CTCF_BOUND_COMPOSITE"])
    print(f"  {len(ctcf_bound_accessions):,} CTCF-bound accessions loaded.\n")

    if args.test:
        chrom, pos = args.test.split(":")
        pos = int(pos)
        print(f"Single-locus test: chr{chrom}:{pos:,}...")
        result = annotate_locus(chrom, pos, class_paths, ctcf_bound_accessions)
        print(json.dumps(result, indent=2) if result else "No cCRE found in window.")
        return

    loci = load_loci()
    print(f"Loaded {len(loci)} loci.\n")

    results = []
    for locus in loci:
        chrom = locus["Chr"]
        pos = int(locus["Pos(GRCh38)"])
        rsid = locus["rsID"]
        print(f"Annotating {rsid} (chr{chrom}:{pos:,})...")
        result = annotate_locus(chrom, pos, class_paths, ctcf_bound_accessions)
        if result is None:
            print("  No cCRE found within window.")
        else:
            print(f"  -> {result['accession']}, class={result['class']}, "
                  f"dist={result['distance_bp']}bp, CTCF-bound={result['ctcf_bound']}")
        results.append({"locus": locus["Locus"], "rsid": rsid, "chr": chrom, "pos": pos, "ccre": result})

    out_path = os.path.join(DATA_DIR, "ccre_annotation_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {len(results)} annotation results to {out_path}")


if __name__ == "__main__":
    main()
