"""
Stage 7: Nearest protein-coding gene (Ensembl REST overlap/region).

Finds the nearest protein-coding gene to a variant within a window (default
500 kb) and reports the distance to the nearest gene edge (0 if the variant lies
inside the gene). The protein_coding filter is applied in the request and again
client-side.

Usage:
    python 07_nearest_gene.py --test rs572547799:6:103805760:500000
    python 07_nearest_gene.py --variants rsid:chrom:pos:window ...
    python 07_nearest_gene.py --table1
--table1 re-derives the nearest gene for every locus in data/table1_final.csv,
compares it with the Gene column, and lists assignments beyond 50 kb.

Output: data/nearest_gene_results.json (--test / --variants; overwritten on each
        run) or data/gene_assignment_check.csv (--table1)
"""

import argparse
import csv
import json
import os
import re
import time

import requests

from config import DATA_DIR

ENSEMBL_OVERLAP_URL = "https://rest.ensembl.org/overlap/region/human/{region}"


def query_nearest_gene(rsid: str, chrom: str, pos: int, window_bp: int):
    region = f"{chrom}:{max(1, pos - window_bp)}-{pos + window_bp}"
    url = ENSEMBL_OVERLAP_URL.format(region=region)
    # Retry on transient failures (5xx errors such as 503, and connection errors);
    # Ensembl's public REST server is intermittently unavailable.
    genes = None
    last_error = None
    max_retries = 4
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(
                url,
                params={"feature": "gene", "biotype": "protein_coding"},
                headers={"Content-Type": "application/json"},
                timeout=30,
            )
            resp.raise_for_status()
            genes = resp.json()
            break
        except requests.RequestException as e:
            last_error = e
            status = getattr(getattr(e, "response", None), "status_code", None)
            transient = status is None or status >= 500 or status == 429
            if transient and attempt < max_retries:
                wait = 3 * (2 ** (attempt - 1))  # 3s, 6s, 12s
                print(f"    [RETRY {attempt}/{max_retries}] {rsid}: {e}; waiting {wait}s...")
                time.sleep(wait)
                continue
            break
    if genes is None:
        print(f"  [WARN] Query failed for {rsid} ({region}) after retries: {last_error}. "
              f"If this is a 5xx error, Ensembl is temporarily unavailable; try again later.")
        return {"rsid": rsid, "chrom": chrom, "pos": pos, "window_bp": window_bp,
                "error": str(last_error)}

    # Client-side filter: do not rely on the server-side biotype parameter alone.
    protein_coding = [g for g in genes if g.get("biotype") == "protein_coding"]
    if len(protein_coding) != len(genes):
        print(f"  [NOTE] Server returned {len(genes)} gene(s), {len(protein_coding)} "
              f"after client-side protein_coding filter -- biotype param may not be "
              f"fully honored server-side; filtering locally regardless.")

    if not protein_coding:
        print(f"  {rsid}: no protein-coding gene within {window_bp:,}bp")
        return {"rsid": rsid, "chrom": chrom, "pos": pos, "window_bp": window_bp,
                "nearest_gene": None, "distance_bp": None}

    def distance(g):
        start, end = g["start"], g["end"]
        if start <= pos <= end:
            return 0
        return min(abs(pos - start), abs(pos - end))

    nearest = min(protein_coding, key=distance)
    dist = distance(nearest)
    print(f"  {rsid}: nearest protein-coding gene = {nearest.get('external_name', nearest.get('id'))} "
          f"({dist:,}bp)")
    return {"rsid": rsid, "chrom": chrom, "pos": pos, "window_bp": window_bp,
            "nearest_gene": nearest.get("external_name", nearest.get("id")),
            "ensembl_id": nearest.get("id"), "distance_bp": dist,
            "gene_start": nearest["start"], "gene_end": nearest["end"]}


def parse_variant_arg(arg: str):
    """Parse 'rsid:chrom:pos:window_bp' into its parts."""
    parts = arg.split(":")
    if len(parts) != 4:
        raise ValueError(f"Expected rsid:chrom:pos:window_bp, got: {arg}")
    rsid, chrom, pos, window_bp = parts
    return rsid, chrom, int(pos), int(window_bp)


def parse_assigned_genes(gene_cell: str):
    """Gene symbol(s) named in a Table 1 Gene cell, or None for a desert cell."""
    cell = gene_cell.strip()
    if cell.startswith("("):
        return None
    m = re.match(r"^[A-Za-z0-9.\-]+(?:/[A-Za-z0-9.\-]+)*", cell)
    return set(m.group(0).split("/")) if m else None


def check_table1(window_bp: int = 500_000):
    """Re-derive the nearest protein-coding gene for all loci in
    data/table1_final.csv and compare with the Gene column."""
    from config import load_loci
    rows = []
    for locus in load_loci():
        res = query_nearest_gene(locus["rsID"], locus["Chr"], int(locus["Pos(GRCh38)"]), window_bp)
        time.sleep(0.4)
        assigned = parse_assigned_genes(locus["Gene"])
        nearest, dist = res.get("nearest_gene"), res.get("distance_bp")
        if "error" in res:
            status = "ERROR (query failed; rerun)"
        elif assigned is None and nearest is None:
            status = "desert confirmed"
        elif assigned is None:
            status = f"CHECK: Table 1 lists no gene, but {nearest} is {dist:,} bp away"
        elif nearest is None:
            status = "CHECK: no protein-coding gene found within window"
        elif nearest in assigned:
            status = "match"
        else:
            status = f"CHECK: code finds {nearest} ({dist:,} bp)"
        rows.append({"locus": locus["Locus"], "rsID": locus["rsID"], "table1_gene": locus["Gene"],
                     "code_nearest_gene": nearest or "", "distance_bp": dist if dist is not None else "",
                     "status": status})
    out_path = os.path.join(DATA_DIR, "gene_assignment_check.csv")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\n=== Summary ===")
    print(f"match: {sum(r['status']=='match' for r in rows)}, "
          f"desert confirmed: {sum(r['status']=='desert confirmed' for r in rows)}, "
          f"to review: {sum(r['status'].startswith(('CHECK','ERROR')) for r in rows)}")
    for r in rows:
        if r["status"].startswith(("CHECK", "ERROR")):
            print(f"  {r['locus']}: {r['status']}")
    far = [r for r in rows if r["status"] == "match" and r["distance_bp"] != "" and r["distance_bp"] > 50_000]
    print(f"\nMatched assignments beyond 50 kb ({len(far)}):")
    for r in sorted(far, key=lambda r: r["distance_bp"]):
        print(f"  {r['locus']}: {r['code_nearest_gene']} at {r['distance_bp']/1000:.1f} kb")
    print(f"\nWrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", help="Single variant: rsid:chrom:pos:window_bp")
    parser.add_argument("--variants", nargs="+",
                         help="Multiple variants: rsid:chrom:pos:window_bp rsid:chrom:pos:window_bp ...")
    parser.add_argument("--table1", action="store_true",
                         help="Check every locus in data/table1_final.csv (500 kb window)")
    args = parser.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)

    if args.table1:
        check_table1()
        return

    variant_args = []
    if args.test:
        variant_args = [args.test]
    elif args.variants:
        variant_args = args.variants
    else:
        parser.error("Specify --test <variant> or --variants <variant> [<variant> ...]")

    results = []
    for v in variant_args:
        rsid, chrom, pos, window_bp = parse_variant_arg(v)
        results.append(query_nearest_gene(rsid, chrom, pos, window_bp))

    out_path = os.path.join(DATA_DIR, "nearest_gene_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {len(results)} results to {out_path}")


if __name__ == "__main__":
    main()
