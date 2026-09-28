"""
Stage 10: Build the 500-gene PheWAS background result file.

Concatenates the two 04_phewas_query.py --genes-file result sets (the original
300 background genes and the 200-gene extension drawn by Stage 9) into
data/phewas_background_n500.json, the file read by 05_background_tests.py.
Stops if the two sets share any gene.

Usage:
    python 10_merge_phewas_background.py

Input:
    data/phewas_results_background_genes_n300.json
    data/phewas_results_background_genes_extension_n200.json
Output:
    data/phewas_background_n500.json
"""
import json
import os

from config import DATA_DIR

parts = ["phewas_results_background_genes_n300.json",
         "phewas_results_background_genes_extension_n200.json"]
merged = []
for name in parts:
    with open(os.path.join(DATA_DIR, name)) as f:
        merged.extend(json.load(f))
genes = [d["gene"] for d in merged]
if len(genes) != len(set(genes)):
    raise SystemExit("Duplicate genes across the two background sets; refusing to merge.")
out = os.path.join(DATA_DIR, "phewas_background_n500.json")
with open(out, "w") as f:
    json.dump(merged, f, indent=2)
print(f"Wrote {len(merged)} genes to {out}")
