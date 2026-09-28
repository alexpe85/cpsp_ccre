"""
Stage 4: Gene-level PheWAS classification (Open Targets Platform).

For each gene, pages through ALL disease associations ranked by association
score and classifies:
    Strong   = a pain-related association at rank <= 5
    Moderate = pain-related association(s) present, none at rank <= 5
    None     = no pain-related association
"Indirect" (tangentially pain-related terms) is assigned by manual review of the
returned hits, not by this script. Pain-related terms are listed in
config.PAIN_TERMS. A full scan is used because a top-5 lookup misses lower-ranked
associations.

Note: the API serves the current Open Targets release; record the release date
when running.

Usage:
    python 04_phewas_query.py --gene GLRA3
    python 04_phewas_query.py --all                    # genes in data/table1_final.csv
    python 04_phewas_query.py --genes-file FILE        # one symbol per line
Output: data/phewas_results.json (--gene / --all), or
        data/phewas_results_<input file name>.json (--genes-file)
"""

import argparse
import json
import os
import time

import requests

from config import load_loci, DATA_DIR, PAIN_TERMS, PHEWAS_STRONG_RANK_CUTOFF

OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
PAGE_SIZE = 100
REQUEST_DELAY_S = 0.3

SEARCH_QUERY = """
query searchGene($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) {
    hits { id name }
  }
}
"""

ASSOCIATIONS_QUERY = """
query geneAssociations($ensemblId: String!, $index: Int!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    associatedDiseases(page: {index: $index, size: $size}) {
      count
      rows {
        disease { id name }
        score
      }
    }
  }
}
"""


def resolve_ensembl_id(gene_symbol: str) -> str | None:
    resp = requests.post(OT_GRAPHQL_URL, json={"query": SEARCH_QUERY, "variables": {"q": gene_symbol}}, timeout=30)
    resp.raise_for_status()
    hits = resp.json().get("data", {}).get("search", {}).get("hits", [])
    for hit in hits:
        if hit.get("name", "").upper() == gene_symbol.upper():
            return hit["id"]
    return hits[0]["id"] if hits else None


def full_scan_associations(ensembl_id: str) -> list[dict]:
    """Page through ALL disease associations for a gene, ranked by score
    (Open Targets returns associatedDiseases pre-sorted by score descending)."""
    all_rows = []
    index = 0
    while True:
        variables = {"ensemblId": ensembl_id, "index": index, "size": PAGE_SIZE}
        resp = requests.post(OT_GRAPHQL_URL, json={"query": ASSOCIATIONS_QUERY, "variables": variables}, timeout=30)
        resp.raise_for_status()
        data = resp.json().get("data", {}).get("target", {})
        if data is None:
            break
        assoc = data.get("associatedDiseases", {})
        rows = assoc.get("rows", [])
        all_rows.extend(rows)
        total = assoc.get("count", 0)
        if len(all_rows) >= total or not rows:
            break
        index += 1
        time.sleep(REQUEST_DELAY_S)
    return all_rows


def classify(rows: list[dict]) -> dict:
    """Apply the rank-based classification rule. rows must already be in
    rank order (Open Targets returns them score-sorted)."""
    pain_hits = []
    for rank, row in enumerate(rows, start=1):
        name = row["disease"]["name"].lower()
        if any(term in name for term in PAIN_TERMS):
            pain_hits.append({
                "rank": rank,
                "disease_name": row["disease"]["name"],
                "disease_id": row["disease"]["id"],
                "score": row["score"],
            })

    if not pain_hits:
        tier = "None"
    elif any(h["rank"] <= PHEWAS_STRONG_RANK_CUTOFF for h in pain_hits):
        tier = "Strong"
    else:
        tier = "Moderate"
    # NOTE: "Indirect" (tangential/comorbidity-adjacent terms) is a judgment
    # call this script does not automate — it requires reading the actual
    # disease names for borderline cases (see e.g. LI-9/FUT8, classified
    # Indirect in the paper). Review pain_hits manually for tier=="Moderate"
    # or "Strong" cases where the matched term is only loosely pain-related.

    return {"tier": tier, "pain_hits": pain_hits, "total_associations": len(rows)}


def run_gene(gene_symbol: str) -> dict:
    print(f"Resolving {gene_symbol}...")
    ensembl_id = resolve_ensembl_id(gene_symbol)
    if ensembl_id is None:
        print(f"  [WARN] Could not resolve {gene_symbol} to an Ensembl ID.")
        return {"gene": gene_symbol, "ensembl_id": None, "result": None}

    print(f"  -> {ensembl_id}. Fetching full association scan...")
    rows = full_scan_associations(ensembl_id)
    result = classify(rows)
    print(f"  -> {result['total_associations']} total associations, "
          f"{len(result['pain_hits'])} pain-related, tier={result['tier']}")
    return {"gene": gene_symbol, "ensembl_id": ensembl_id, "result": result}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", help="Single gene symbol to query")
    parser.add_argument("--all", action="store_true", help="Query every gene in table1_final.csv")
    parser.add_argument("--genes-file", help="Query every gene symbol listed in this file (one per "
                         "line) — use for the 300-gene background set")
    args = parser.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)

    if args.gene:
        results = [run_gene(args.gene)]
        out_path = os.path.join(DATA_DIR, "phewas_results.json")
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {len(results)} PheWAS results to {out_path}")

    elif args.all:
        loci = load_loci()
        genes = sorted({l["Gene"].split(" ")[0].split("/")[0]
                         for l in loci if l["Gene"] and "none within" not in l["Gene"]})
        print(f"Querying {len(genes)} unique genes.\n")
        results = []
        out_path = os.path.join(DATA_DIR, "phewas_results.json")
        for gene in genes:
            results.append(run_gene(gene))
            time.sleep(REQUEST_DELAY_S)
            # Incremental save — a 45-gene run is short enough that this
            # matters less than for --genes-file's 300, but cheap to do
            # consistently and avoids a repeat of Stage 3's all-work-lost bug.
            with open(out_path, "w") as f:
                json.dump(results, f, indent=2)
        print(f"\nWrote {len(results)} PheWAS results to {out_path}")

    elif args.genes_file:
        # Resolve the path robustly: try as given (relative to wherever the
        # user is standing when they run this), then fall back to
        # DATA_DIR-relative — so `--genes-file background_genes_n300.txt`
        # and `--genes-file data/background_genes_n300.txt` both work
        # regardless of whether you're running from scripts/ or the repo
        # root, instead of silently depending on cwd like a plain open()
        # call would.
        genes_file_path = args.genes_file
        if not os.path.exists(genes_file_path):
            candidate = os.path.join(DATA_DIR, os.path.basename(args.genes_file))
            if os.path.exists(candidate):
                genes_file_path = candidate
            else:
                raise FileNotFoundError(
                    f"Could not find '{args.genes_file}' as given, or as "
                    f"'{candidate}'. Pass either a full/relative path from "
                    f"your current directory, or just the filename (e.g. "
                    f"background_genes_n300.txt) to have it looked up in {DATA_DIR}."
                )
        with open(genes_file_path) as f:
            genes = [line.strip() for line in f if line.strip()]
        print(f"Querying {len(genes)} genes from {genes_file_path}.\n")
        results = []
        # Output filename is derived from the INPUT filename, not
        # hardcoded — a hardcoded "phewas_background_n300.json" here
        # silently overwrote a real 300-gene result set with a later
        # 200-gene extension run that used a different --genes-file, since
        # both runs wrote to the same fixed name regardless of input.
        input_stem = os.path.splitext(os.path.basename(genes_file_path))[0]
        out_path = os.path.join(DATA_DIR, f"phewas_results_{input_stem}.json")
        for i, gene in enumerate(genes, 1):
            if i % 25 == 0 or i == 1:
                print(f"[{i}/{len(genes)}]")
            results.append(run_gene(gene))
            time.sleep(REQUEST_DELAY_S)
            # Incremental save is important here — a 300-gene run takes
            # long enough that losing it all to a late failure (as
            # happened once already with Stage 3's much shorter 15-locus
            # run) would be a real waste of an hour-plus of API calls.
            with open(out_path, "w") as f:
                json.dump(results, f, indent=2)
        print(f"\nWrote {len(results)} PheWAS results to {out_path}")

    else:
        parser.error("Specify --gene <SYMBOL>, --all, or --genes-file <FILE>")


if __name__ == "__main__":
    main()
