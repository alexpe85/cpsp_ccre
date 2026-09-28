"""
Stage 3: CTCF motif scan.

Scans the sequence of the nearest cCRE at each CTCF-bound locus for the
canonical CTCF motif (JASPAR MA0139.1) using a position weight matrix
(Biopython; both strands; log-odds score in bits; threshold 8.45 bits,
p < 1e-4) and reports the best-scoring window. This is a presence/absence check
against one known motif, not de novo motif discovery.

Sources: JASPAR REST API (matrix); Ensembl REST /sequence/region (GRCh38
sequence). Requests are retried on 5xx errors; results are written after every
locus.

Usage:  python 03_ctcf_motif_scan.py
Input:  data/ccre_annotation_results.json (Stage 2)
Output: data/ctcf_motif_results.json
"""

import json
import os
import time

import requests
from Bio import motifs
from Bio.Seq import Seq

from config import DATA_DIR, CTCF_MOTIF_BIT_THRESHOLD

JASPAR_MATRIX_ID = "MA0139.1"
JASPAR_API_URL = f"https://jaspar.elixir.no/api/v1/matrix/{JASPAR_MATRIX_ID}/"
ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/region/human/{region}"
REQUEST_DELAY_S = 0.5


def fetch_jaspar_pwm() -> motifs.Motif:
    """Fetch the MA0139.1 position-frequency matrix from JASPAR and build a
    Biopython motif object with log-odds scoring."""
    resp = requests.get(JASPAR_API_URL, timeout=30)
    resp.raise_for_status()
    pfm_data = resp.json()["pfm"]  # {"A": [...], "C": [...], "G": [...], "T": [...]}

    counts = {base: [int(round(v)) for v in vals] for base, vals in pfm_data.items()}
    motif = motifs.Motif(alphabet="ACGT", counts=counts)
    motif.pseudocounts = 0.8  # standard JASPAR-recommended pseudocount
    motif.background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    return motif


def fetch_sequence(chrom: str, start: int, end: int, max_retries: int = 3) -> str | None:
    region = f"{chrom}:{start}-{end}"
    url = ENSEMBL_SEQUENCE_URL.format(region=region)
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(url, headers={"Content-Type": "application/json"},
                                 params={"coord_system_version": "GRCh38"}, timeout=30)
            resp.raise_for_status()
            return resp.json().get("seq")
        except requests.RequestException as e:
            is_last = attempt == max_retries
            # Retry on transient server errors (5xx) and connection issues;
            # a real 4xx (bad request) won't be fixed by retrying, so don't
            # waste time on those.
            status = getattr(e.response, "status_code", None)
            transient = status is None or status >= 500
            if transient and not is_last:
                wait = 2 ** attempt  # 2s, 4s, 8s
                print(f"  [RETRY {attempt}/{max_retries}] {region} failed ({e}); "
                      f"waiting {wait}s before retry...")
                time.sleep(wait)
                continue
            print(f"  [WARN] Sequence fetch failed for {region}: {e}")
            return None


def best_match_score(motif_obj: motifs.Motif, sequence: str) -> tuple[float, str, int, str]:
    """
    Scan both strands, return (best_score_bits, matched_subsequence, position, strand).
    """
    pssm = motif_obj.pssm
    seq_obj = Seq(sequence.upper())

    best_score = float("-inf")
    best_pos = None
    best_strand = "+"

    for strand, s in [("+", seq_obj), ("-", seq_obj.reverse_complement())]:
        for pos, score in pssm.search(s, threshold=float("-inf"), both=False):
            # Biopython's PSSM search can return numpy scalar types; cast to native
            # Python types so results are JSON-serializable (numpy.bool_ is not).
            score = float(score)
            pos = int(pos)
            if score > best_score:
                best_score = score
                best_pos = pos if pos >= 0 else len(s) + pos
                best_strand = strand

    motif_len = pssm.length
    if best_strand == "+":
        matched = str(seq_obj[best_pos:best_pos + motif_len])
    else:
        matched = str(seq_obj.reverse_complement()[best_pos:best_pos + motif_len])

    return best_score, matched, best_pos, best_strand


def main():
    ccre_path = os.path.join(DATA_DIR, "ccre_annotation_results.json")
    if not os.path.exists(ccre_path):
        raise FileNotFoundError(f"{ccre_path} not found — run 02_ccre_annotation.py first.")

    with open(ccre_path) as f:
        annotations = json.load(f)

    ctcf_bound = [a for a in annotations if a.get("ccre") and a["ccre"].get("ctcf_bound")]
    print(f"{len(ctcf_bound)} CTCF-bound loci to scan.\n")

    print("Fetching MA0139.1 PWM from JASPAR...")
    motif_obj = fetch_jaspar_pwm()
    print(f"  Motif length: {motif_obj.length}bp\n")

    out_path = os.path.join(DATA_DIR, "ctcf_motif_results.json")
    results = []
    for a in ctcf_bound:
        ccre = a["ccre"]
        # ccre dict now comes from the V4 bulk-BED annotation (start/end
        # directly), not the old V3 GraphQL output (start/len) — support
        # both so this script works against either an old or new
        # ccre_annotation_results.json without silently misreading one.
        chrom, start = a["chr"], ccre["start"]
        end = ccre["end"] if "end" in ccre else start + ccre["len"]
        print(f"Scanning {a['rsid']} ({a['locus']}) cCRE {ccre['accession']} "
              f"(chr{chrom}:{start}-{end})...")

        seq = fetch_sequence(chrom, start, end)
        time.sleep(REQUEST_DELAY_S)
        if seq is None:
            results.append({**a, "motif_scan": None})
            with open(out_path, "w") as f:
                json.dump(results, f, indent=2)
            continue

        score, matched_seq, pos, strand = best_match_score(motif_obj, seq)
        meets_threshold = bool(score >= CTCF_MOTIF_BIT_THRESHOLD)  # native bool, not numpy.bool_
        print(f"  -> best score {score:.2f} bits, threshold_met={meets_threshold}")

        results.append({
            **a,
            "motif_scan": {
                "best_score_bits": round(float(score), 4),
                "threshold_met": meets_threshold,
                "matched_sequence": matched_seq,
                "position": int(pos),
                "strand": strand,
            },
        })

        # Write after every locus, not just once at the end — an earlier
        # version of this script lost all 15 completed scans to a
        # serialization bug that only surfaced at the final json.dump().
        # Incremental writes mean a late failure loses at most one locus's
        # result, not the whole run.
        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)

    print(f"\nWrote {len(results)} motif scan results to {out_path}")


if __name__ == "__main__":
    main()
