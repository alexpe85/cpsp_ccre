"""
Stage 1b: Direct genotype-level LD (r^2, D') from phased 1000 Genomes data.

Computes r^2 and D' between two variants from phased haplotypes in the GRCh38
liftover of 1000 Genomes phase 3 (EBI-hosted, tabix-indexed; only the
requested regions are read). Works at any distance, unlike Ensembl's pairwise
LD endpoint (01_ld_resolution.py), which returned no result for the long-range
pairs analysed here.

    D = p_AB - p_A*p_B;  r^2 = D^2 / [p_A(1-p_A) p_B(1-p_B)];  D' = |D| / Dmax

Validation: the run starts with a self-test on a pair whose r^2 Ensembl
reports (rs6792369 / rs1042779, KHV; 0.9516) and stops if the direct
computation differs by 0.05 or more (observed difference: 0.0002).

Pairs computed: the two long-range independence checks, and the within-cluster
pairs quoted in Results and Figure 2 (positions resolved from rsIDs via the
Ensembl variation endpoint). Each recomputed r^2 is printed next to the value
reported in the manuscript. The UK Biobank LD value for MTF2/FNBP1L (r^2 = 0.780)
was obtained with the authors' own UK Biobank access and is not reproduced here
(see README).

Requires: `tabix` (htslib) on PATH; access to ftp.1000genomes.ebi.ac.uk and
rest.ensembl.org.

Usage:  python 01b_ld_longrange_direct.py
Output: data/ld_longrange_results.json
"""

import json
import os
import subprocess
import time

import requests

from config import DATA_DIR

VCF_URL_TEMPLATE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/"
    "phase3.chr{chrom}.GRCh38.GT.crossmap.vcf.gz"
)
PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

# Long-range independence checks (Ensembl's pairwise endpoint returned no result
# for either pair).
LONGRANGE_PAIRS = [
    {"rsid_a": "rs76991866", "chr": "17", "pos_a": 5076637,
     "rsid_b": "rs11653414", "pos_b": 5511718, "population": "CEU", "reported_r2": 0.006},
    {"rsid_a": "rs186819635", "chr": "1", "pos_a": 94605623,
     "rsid_b": "rs151142747", "pos_b": 93517560, "population": "EUR", "reported_r2": 0.000006},
]

# Within-cluster pairs quoted in Results and Figure 2. Positions are resolved from
# rsIDs at run time; reported_r2 is the value given in the manuscript, printed
# next to the recomputed value.
REPORTED_PAIRS = [
    {"rsid_a": "rs76991866", "rsid_b": "rs114604537", "chr": "17", "population": "CEU", "reported_r2": 0.65},
    {"rsid_a": "rs76991866", "rsid_b": "rs3026120", "chr": "17", "population": "CEU", "reported_r2": 0.12},
    {"rsid_a": "rs115024240", "rsid_b": "rs151142747", "chr": "1", "population": "CEU", "reported_r2": 1.0},
    {"rsid_a": "rs17298280", "rsid_b": "rs10032594", "chr": "4", "population": "EUR", "reported_r2": 0.487},
    {"rsid_a": "rs186819635", "rsid_b": "rs145636748", "chr": "1", "population": "EUR", "reported_r2": 0.371},
]

# Ensembl's own documented example pair, short-range, KNOWN to return
# r2=0.951626 via the REST API. Used to validate this script before it is
# applied to the pairs above.
# Positions are resolved automatically at runtime via Ensembl's variation
# lookup (rather than hardcoded here) to avoid trusting a manually-copied
# coordinate for the one pair this whole approach is validated against.
SELF_TEST_RSIDS = ("rs6792369", "rs1042779")
SELF_TEST_POPULATION = "KHV"
SELF_TEST_EXPECTED_R2 = 0.9516


def resolve_variant_position(rsid: str) -> tuple[str, int]:
    """Look up an rsID's GRCh38 chr:pos via Ensembl's variation endpoint."""
    url = f"https://rest.ensembl.org/variation/human/{rsid}"
    for attempt in range(1, 4):
        try:
            resp = requests.get(url, headers={"Content-Type": "application/json"},
                                 params={"genotypes": 0}, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            break
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    # mappings[] can contain multiple entries (e.g. patches); take the first
    # primary-assembly GRCh38 mapping.
    for mapping in data.get("mappings", []):
        if mapping.get("assembly_name", "").startswith("GRCh38") or "assembly_name" not in mapping:
            return mapping["seq_region_name"], int(mapping["start"])
    raise ValueError(f"No GRCh38 mapping found for {rsid}")


_PANEL_FILE_CACHE: str | None = None
_POPULATION_SAMPLES_CACHE: dict[str, list[str]] = {}


def fetch_population_samples(population_or_superpop: str, max_retries: int = 3) -> list[str]:
    """
    Return the list of 1000G sample IDs belonging to a population (e.g. CEU)
    or super-population (e.g. EUR). The panel file has columns:
    sample  pop  super_pop  gender

    Retries on transient network failures. Caches the raw panel file and the
    per-population sample lists, so a multi-locus run fetches the panel once.
    """
    global _PANEL_FILE_CACHE

    if population_or_superpop in _POPULATION_SAMPLES_CACHE:
        return _POPULATION_SAMPLES_CACHE[population_or_superpop]

    if _PANEL_FILE_CACHE is None:
        last_error = None
        for attempt in range(1, max_retries + 1):
            try:
                resp = requests.get(PANEL_URL, timeout=30)
                resp.raise_for_status()
                _PANEL_FILE_CACHE = resp.text
                break
            except requests.RequestException as e:
                last_error = e
                if attempt < max_retries:
                    wait = 2 ** attempt
                    print(f"    [RETRY {attempt}/{max_retries}] panel file fetch failed ({e}); "
                          f"waiting {wait}s before retry...")
                    time.sleep(wait)
        else:
            raise RuntimeError(f"Could not fetch {PANEL_URL} after {max_retries} attempts: {last_error}")

    samples = []
    for line in _PANEL_FILE_CACHE.strip().split("\n")[1:]:  # skip header
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        sample, pop, super_pop = parts[0], parts[1], parts[2]
        if population_or_superpop in (pop, super_pop):
            samples.append(sample)

    _POPULATION_SAMPLES_CACHE[population_or_superpop] = samples
    return samples


def get_haplotypes(chrom: str, pos: int, rsid: str, samples: list[str]) -> list[int] | None:
    """
    Query the remote tabix-indexed VCF for a small window around pos via the
    standalone `tabix` CLI, find the record matching
    rsid (or pos, as fallback), and return a flat list of 0/1 alleles across
    all haplotypes (2 per sample, phased) for the given sample subset.
    Returns None if the variant isn't found.

    Uses the standalone tabix binary for remote region queries.
    """
    # First, determine correct contig naming by listing this file's contigs
    # via a lightweight tabix -l call (chromosome names only, no records).
    vcf_url = VCF_URL_TEMPLATE.format(chrom=chrom)
    try:
        list_result = subprocess.run(
            ["tabix", "-l", vcf_url], capture_output=True, text=True, timeout=60
        )
        available_contigs = set(list_result.stdout.strip().split("\n"))
    except FileNotFoundError:
        raise RuntimeError(
            "`tabix` command not found. Install htslib first (see README): "
            "conda install -c bioconda -c conda-forge htslib; "
            "apt install tabix; or brew install htslib."
        )
    except subprocess.TimeoutExpired:
        print(f"  [WARN] tabix -l timed out listing contigs for {vcf_url}")
        return None

    if chrom in available_contigs:
        contig = chrom
    elif f"chr{chrom}" in available_contigs:
        contig = f"chr{chrom}"
    else:
        print(f"  [WARN] Neither '{chrom}' nor 'chr{chrom}' found in this file's "
              f"contigs ({sorted(available_contigs)[:5]}...). Cannot query.")
        return None

    window = 5000
    region = f"{contig}:{max(1, pos - window)}-{pos + window}"

    try:
        result = subprocess.run(
            ["tabix", vcf_url, region], capture_output=True, text=True, timeout=120
        )
    except subprocess.TimeoutExpired:
        print(f"  [WARN] tabix query timed out for region {region}")
        return None

    if result.returncode != 0:
        print(f"  [WARN] tabix exited with error for region {region}: {result.stderr.strip()}")
        return None

    lines = [l for l in result.stdout.strip().split("\n") if l]
    if not lines:
        print(f"  [WARN] Zero records returned for region {region} — check "
              f"contig naming ('{contig}') and that pos={pos} is genuinely "
              f"GRCh38 (not accidentally GRCh37).")
        return None

    # We need the sample column order from the file header to know which
    # columns in each data line correspond to `samples`. tabix -h would give
    # us the header, but a second network round-trip is wasteful — instead
    # fetch the header once via `tabix -H` (samtools >=1.11) or fall back to
    # `bcftools view -h` if that's unavailable. Simplest robust option:
    # query the header once per chromosome file and cache it.
    header_samples = _get_vcf_sample_order(vcf_url)
    sample_indices = [header_samples.index(s) for s in samples if s in header_samples]
    missing = set(samples) - set(header_samples)
    if missing:
        print(f"  [WARN] {len(missing)} requested samples not found in VCF header "
              f"(e.g. {sorted(missing)[:3]}) — proceeding with the {len(sample_indices)} found.")

    for line in lines:
        fields = line.split("\t")
        vcf_id, vcf_pos = fields[2], int(fields[1])
        if vcf_id != rsid and vcf_pos != pos:
            continue
        ref, alt = fields[3], fields[4]
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue  # skip multiallelic/indel for a clean biallelic r^2

        format_field = fields[8].split(":")
        gt_idx = format_field.index("GT")
        sample_fields = fields[9:]

        haplotypes = []
        for idx in sample_indices:
            gt_str = sample_fields[idx].split(":")[gt_idx]
            sep = "|" if "|" in gt_str else "/"
            if sep == "/":
                print(f"  [WARN] Unphased genotype encountered ({gt_str}) at "
                      f"{vcf_id or vcf_pos} — this file is expected to be fully "
                      f"phased; treating as phased anyway, but flag this if it recurs.")
            alleles = gt_str.split(sep)
            haplotypes.extend(int(a) for a in alleles if a in ("0", "1"))
        return haplotypes

    print(f"  [WARN] Records found in window, but none matched rsID={rsid} "
          f"or pos={pos} exactly — variant may not be in this call set, "
          f"or its exact position differs slightly from what was queried.")
    return None


_HEADER_CACHE: dict[str, list[str]] = {}


def _get_vcf_sample_order(vcf_url: str, max_retries: int = 3) -> list[str]:
    """Fetch and cache the sample column order from a remote VCF's header line.
    Retries on transient failures (empty or malformed output)."""
    if vcf_url in _HEADER_CACHE:
        return _HEADER_CACHE[vcf_url]

    last_error = None
    for attempt in range(1, max_retries + 1):
        result = subprocess.run(["tabix", "-H", vcf_url], capture_output=True, text=True, timeout=60)
        header_line = next(
            (l for l in result.stdout.split("\n") if l.startswith("#CHROM")), None
        )
        if header_line is not None:
            columns = header_line.lstrip("#").split("\t")
            samples = columns[9:]  # standard VCF: first 9 columns are fixed fields
            _HEADER_CACHE[vcf_url] = samples
            return samples
        last_error = result.stderr.strip() or "(empty stdout, no #CHROM line, no stderr)"
        if attempt < max_retries:
            wait = 2 ** attempt
            print(f"    [RETRY {attempt}/{max_retries}] tabix -H returned no #CHROM line for "
                  f"{vcf_url}; waiting {wait}s before retry...")
            time.sleep(wait)

    raise RuntimeError(
        f"Could not find #CHROM header line for {vcf_url} after {max_retries} attempts. "
        f"Last error/output: {last_error}. `tabix -H` output may differ by version — "
        f"inspect manually with `tabix -H {vcf_url} | tail -5`."
    )


def compute_ld(hap_a: list[int], hap_b: list[int]) -> dict:
    """
    Exact r^2 and D' from two equal-length phased haplotype allele lists
    (paired by position — hap_a[i] and hap_b[i] must come from the same
    physical haplotype).
    """
    n = len(hap_a)
    assert n == len(hap_b) and n > 0

    freq_a = sum(hap_a) / n
    freq_b = sum(hap_b) / n
    freq_ab = sum(1 for a, b in zip(hap_a, hap_b) if a == 1 and b == 1) / n

    D = freq_ab - freq_a * freq_b

    if D >= 0:
        d_max = min(freq_a * (1 - freq_b), (1 - freq_a) * freq_b)
    else:
        d_max = min(freq_a * freq_b, (1 - freq_a) * (1 - freq_b))
    d_prime = abs(D / d_max) if d_max != 0 else None

    denom = freq_a * (1 - freq_a) * freq_b * (1 - freq_b)
    r2 = (D ** 2) / denom if denom != 0 else None

    return {
        "n_haplotypes": n, "freq_a": round(freq_a, 6), "freq_b": round(freq_b, 6),
        "r2": round(r2, 6) if r2 is not None else None,
        "d_prime": round(d_prime, 6) if d_prime is not None else None,
    }


def run_pair(pair: dict) -> dict:
    print(f"\n{pair['rsid_a']} vs {pair['rsid_b']} (population={pair['population']})...")
    samples = fetch_population_samples(pair["population"])
    print(f"  {len(samples)} samples in {pair['population']}")

    hap_a = get_haplotypes(pair["chr"], pair["pos_a"], pair["rsid_a"], samples)
    hap_b = get_haplotypes(pair["chr"], pair["pos_b"], pair["rsid_b"], samples)

    if hap_a is None or hap_b is None:
        missing = pair["rsid_a"] if hap_a is None else pair["rsid_b"]
        print(f"  [WARN] Variant not found in VCF: {missing}")
        return {**pair, "ld": None, "error": f"{missing} not found"}

    ld = compute_ld(hap_a, hap_b)
    print(f"  -> r2={ld['r2']}, D'={ld['d_prime']}, n_haplotypes={ld['n_haplotypes']}")
    return {**pair, "ld": ld}


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    print("=" * 60)
    print("SELF-TEST: known short-range pair (should match Ensembl's r2~=0.9516)")
    print("=" * 60)
    chr_a, pos_a = resolve_variant_position(SELF_TEST_RSIDS[0])
    chr_b, pos_b = resolve_variant_position(SELF_TEST_RSIDS[1])
    assert chr_a == chr_b, f"Self-test pair unexpectedly on different chromosomes: {chr_a} vs {chr_b}"
    self_test_pair = {
        "rsid_a": SELF_TEST_RSIDS[0], "chr": chr_a, "pos_a": pos_a,
        "rsid_b": SELF_TEST_RSIDS[1], "pos_b": pos_b, "population": SELF_TEST_POPULATION,
    }
    test_result = run_pair(self_test_pair)
    if test_result["ld"] is None:
        raise RuntimeError(
            "Self-test FAILED — could not compute LD for Ensembl's own known "
            "pair. Do not trust results below until this is fixed."
        )
    diff = abs(test_result["ld"]["r2"] - SELF_TEST_EXPECTED_R2)
    if diff >= 0.05:
        raise AssertionError(
            f"Self-test FAILED — got r2={test_result['ld']['r2']}, expected "
            f"~{SELF_TEST_EXPECTED_R2}. Do not trust results below until this "
            f"discrepancy (diff={diff:.4f}) is understood — likely a phasing, "
            f"sample-subset, or allele-orientation bug in this script."
        )
    print(f"Self-test PASSED (diff from expected: {diff:.4f}) — "
          f"proceeding to long-range pairs.\n")

    print("=" * 60)
    print("Pairs reported in the manuscript")
    print("=" * 60)
    pairs = list(LONGRANGE_PAIRS)
    for p in REPORTED_PAIRS:
        chr_a, pos_a = resolve_variant_position(p["rsid_a"])
        chr_b, pos_b = resolve_variant_position(p["rsid_b"])
        if chr_a != p["chr"] or chr_b != p["chr"]:
            print(f"  [WARN] {p['rsid_a']}/{p['rsid_b']}: resolved on chr{chr_a}/chr{chr_b}, "
                  f"expected chr{p['chr']}; skipped")
            continue
        pairs.append({**p, "pos_a": pos_a, "pos_b": pos_b})
    results = [run_pair(p) for p in pairs]

    print("\n" + "=" * 60)
    print("Recomputed vs reported r2")
    print("=" * 60)
    for r in results:
        if r.get("ld") is None or "reported_r2" not in r:
            continue
        got, want = r["ld"]["r2"], r["reported_r2"]
        flag = "ok" if abs(got - want) <= 0.02 else "DIFFERS"
        print(f"  {r['rsid_a']} / {r['rsid_b']} ({r['population']}): "
              f"recomputed {got} vs reported {want}  [{flag}]")

    out_path = os.path.join(DATA_DIR, "ld_longrange_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {len(results)} results to {out_path}")


if __name__ == "__main__":
    main()
