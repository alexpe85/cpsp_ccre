"""
Stage 11: Compare Stage 4 gene-level results with the PheWAS column of Table 1.

For every locus, the tier returned by 04_phewas_query.py for the assigned gene is
compared with the tier in data/table1_final.csv (footnote markers removed).
"Indirect" is assigned by manual review of the returned hits and cannot be
reproduced by the script, so those loci are listed separately. Gene deserts have
no gene to query.

Usage:
    python 04_phewas_query.py --all
    python 11_check_phewas_table1.py
Input:  data/table1_final.csv, data/phewas_results.json
Output: data/phewas_table1_check.csv
"""
import csv
import json
import os

from config import DATA_DIR, load_loci

MARKERS = "*\u2020\u2021\u00a7\u2016\u00b6#"


def strip_markers(value: str) -> str:
    return value.rstrip(MARKERS).strip()


def main():
    with open(os.path.join(DATA_DIR, "phewas_results.json")) as f:
        results = {d["gene"]: (d["result"] or {}).get("tier") for d in json.load(f)}

    rows = []
    for locus in load_loci():
        cell = locus["Gene"]
        table_tier = strip_markers(locus["PheWAS"])
        if "none within" in cell:
            rows.append((locus["Locus"], "", table_tier, "", "gene desert (no gene to query)"))
            continue
        gene = cell.split(" ")[0].split("/")[0]
        script_tier = results.get(gene)
        if script_tier is None:
            status = "not in phewas_results.json (run 04 --all)"
        elif table_tier == "Indirect":
            status = f"manual class 'Indirect'; script tier = {script_tier}"
        elif script_tier == table_tier:
            status = "match"
        else:
            status = "MISMATCH"
        rows.append((locus["Locus"], gene, table_tier, script_tier or "", status))

    out = os.path.join(DATA_DIR, "phewas_table1_check.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["locus", "gene", "table1_tier", "script_tier", "status"])
        w.writerows(rows)

    n_match = sum(r[4] == "match" for r in rows)
    print(f"match: {n_match}, gene deserts: {sum(r[4].startswith('gene desert') for r in rows)}, "
          f"manual Indirect: {sum(r[4].startswith('manual') for r in rows)}, "
          f"not queried: {sum(r[4].startswith('not in') for r in rows)}, "
          f"MISMATCH: {sum(r[4] == 'MISMATCH' for r in rows)}")
    for r in rows:
        if r[4] == "MISMATCH" or r[4].startswith("not in"):
            print(f"  {r[0]} ({r[1]}): Table 1 = {r[2]}, script = {r[3] or 'n/a'}  [{r[4]}]")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
