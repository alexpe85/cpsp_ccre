"""
Shared configuration and locus list for the CPSP cCRE annotation pipeline.

Loads the working locus set from data/table1_final.csv so every stage
operates on the same canonical list. If that file is absent, falls back to
CANDIDATE_LOCI below (the pre-LD-resolution compilation) so 01_ld_resolution.py
can be run as the very first step of a from-scratch re-analysis.
"""

import os
import csv

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TABLE1_PATH = os.path.join(DATA_DIR, "table1_final.csv")

GENOME_BUILD = "GRCh38"

# Reference panel per source study, matched to each study's reported ancestry.
# See Methods, "Linkage Disequilibrium Resolution".
PANEL_BY_SOURCE = {
    "Parisien": "1000GENOMES:phase_3:CEU",
    "own_GWAS": "1000GENOMES:phase_3:CEU",  # cross-checked against UKBB-LD where available; see README
    "Li_BrJAnaesth": "1000GENOMES:phase_3:EUR",
    # Li_Anaesthesia uses the same UK Biobank source as Li_BrJAnaesth (EUR). Warner's
    # cohort is Nottinghamshire, UK; 1000G GBR is the closest population-matched panel.
    "Li_Anaesthesia": "1000GENOMES:phase_3:EUR",
    "Warner": "1000GENOMES:phase_3:GBR",
}

# cCRE search window: nearest element within +/-50 kb. Gene assignment uses a separate rule
# (nearest protein-coding gene within 50 kb, extended to 500 kb where none is found).
CCRE_WINDOW_BP = 50_000

# Gene-assignment window: strict for the 8-locus Li et al. Br J Anaesth batch,
# extended for the original 37-locus set where the initial window was empty.
GENE_WINDOW_STRICT_BP = 50_000
GENE_WINDOW_EXTENDED_BP = 500_000

# CTCF-bound threshold (Methods)
CTCF_MAXZ_THRESHOLD = 1.64

# CTCF motif (JASPAR MA0139.1) match threshold, bits (Methods: p<1e-4)
CTCF_MOTIF_BIT_THRESHOLD = 8.45

# PheWAS classification rule (rank-based; see Methods, "Gene-Level Phenotype
# Association Analysis" and the corrections discussed for WAR-1/FRMD4A/NXPH2):
#   Strong   = a pain-related hit at rank <= 5 anywhere in the full scan
#   Moderate = a pain-related hit present, but no hit at rank <= 5
#   Indirect = only tangential/comorbidity-adjacent pain terms present
#   None     = no pain-related hit found in the full scan
PHEWAS_STRONG_RANK_CUTOFF = 5
PAIN_TERMS = ["pain", "migraine", "fibromyalgia", "neuropathic", "headache"]


def load_loci(path=TABLE1_PATH):
    """Load the 45-locus table. Returns a list of dicts, one per locus."""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found. Run 01_ld_resolution.py first, or supply "
            "table1_final.csv manually — see README."
        )
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# Pre-LD-resolution candidate compilation (50 lead variants), for a genuine
# from-scratch re-run of 01_ld_resolution.py. chr/pos are GRCh38.
# Source label maps to PANEL_BY_SOURCE above.
CANDIDATE_LOCI = [
    # --- Parisien et al. 2024 (24 loci) ---
    {"rsid": "rs143994530", "chr": "2", "pos": 117579040, "source": "Parisien"},
    {"rsid": "rs78454748", "chr": "5", "pos": 30194366, "source": "Parisien"},
    {"rsid": "rs74896012", "chr": "6", "pos": 1016856, "source": "Parisien"},
    {"rsid": "rs117244643", "chr": "6", "pos": 96807356, "source": "Parisien"},
    {"rsid": "rs10488532", "chr": "7", "pos": 93135176, "source": "Parisien"},
    {"rsid": "rs10969869", "chr": "9", "pos": 30797793, "source": "Parisien"},
    {"rsid": "rs78326996", "chr": "9", "pos": 96697147, "source": "Parisien"},
    {"rsid": "rs17882261", "chr": "10", "pos": 79565219, "source": "Parisien"},
    {"rsid": "rs76985919", "chr": "14", "pos": 94409210, "source": "Parisien"},
    {"rsid": "rs4347600", "chr": "15", "pos": 90521232, "source": "Parisien"},
    {"rsid": "rs79423608", "chr": "16", "pos": 25066103, "source": "Parisien"},
    {"rsid": "rs140299794", "chr": "19", "pos": 24176470, "source": "Parisien"},
    {"rsid": "rs142578347", "chr": "17", "pos": 52937473, "source": "Parisien"},
    {"rsid": "rs112030385", "chr": "18", "pos": 30885936, "source": "Parisien"},
    {"rsid": "rs111290518", "chr": "7", "pos": 79167224, "source": "Parisien"},
    {"rsid": "rs7946537", "chr": "11", "pos": 68550719, "source": "Parisien"},
    {"rsid": "rs76991866", "chr": "17", "pos": 5076637, "source": "Parisien"},
    {"rsid": "rs114604537", "chr": "17", "pos": None, "source": "Parisien"},  # excluded via LD (17p13.1 cluster)
    {"rsid": "rs3026120", "chr": "17", "pos": None, "source": "Parisien"},    # excluded via LD (17p13.1 cluster)
    {"rsid": "rs11653414", "chr": "17", "pos": 5511718, "source": "Parisien"},
    {"rsid": "rs4142128", "chr": "20", "pos": 26269908, "source": "Parisien"},
    {"rsid": "rs1596479", "chr": "18", "pos": 65522255, "source": "Parisien"},
    {"rsid": "rs138190025", "chr": "2", "pos": 51168346, "source": "Parisien"},
    {"rsid": "rs114837251", "chr": "4", "pos": 155683957, "source": "Parisien"},
    # --- Li et al. 2025, Anaesthesia (11 loci) ---
    {"rsid": "rs17047504", "chr": "1", "pos": 218203154, "source": "Li_Anaesthesia"},
    {"rsid": "rs6531281", "chr": "2", "pos": 17093972, "source": "Li_Anaesthesia"},
    {"rsid": "rs13127505", "chr": "4", "pos": 140959616, "source": "Li_Anaesthesia"},
    {"rsid": "rs56052023", "chr": "5", "pos": 60108542, "source": "Li_Anaesthesia"},
    {"rsid": "rs182762077", "chr": "5", "pos": 112497128, "source": "Li_Anaesthesia"},
    {"rsid": "rs116169715", "chr": "6", "pos": 23158218, "source": "Li_Anaesthesia"},
    {"rsid": "rs146141654", "chr": "7", "pos": 105368442, "source": "Li_Anaesthesia"},
    {"rsid": "rs185545327", "chr": "8", "pos": 23734934, "source": "Li_Anaesthesia"},
    {"rsid": "rs78134813", "chr": "14", "pos": 65292034, "source": "Li_Anaesthesia"},
    {"rsid": "rs117920312", "chr": "14", "pos": 68596205, "source": "Li_Anaesthesia"},
    {"rsid": "rs4843341", "chr": "16", "pos": 86074561, "source": "Li_Anaesthesia"},
    # --- Li et al. 2025, Br J Anaesth (10 raw lead SNPs -> 8 after within-study LD) ---
    {"rsid": "rs17298280", "chr": "4", "pos": 174713747, "source": "Li_BrJAnaesth"},
    {"rsid": "rs10032594", "chr": "4", "pos": None, "source": "Li_BrJAnaesth"},   # excluded, LD with rs17298280
    {"rsid": "rs186819635", "chr": "1", "pos": 94605623, "source": "Li_BrJAnaesth"},
    {"rsid": "rs145636748", "chr": "1", "pos": None, "source": "Li_BrJAnaesth"},  # excluded, LD with rs186819635
    {"rsid": "rs184832856", "chr": "2", "pos": 138798261, "source": "Li_BrJAnaesth"},
    {"rsid": "rs140330443", "chr": "16", "pos": 1622581, "source": "Li_BrJAnaesth"},
    {"rsid": "rs12143186", "chr": "1", "pos": 68272150, "source": "Li_BrJAnaesth"},
    {"rsid": "rs117130005", "chr": "17", "pos": 19701301, "source": "Li_BrJAnaesth"},
    {"rsid": "rs2160419", "chr": "13", "pos": 71740214, "source": "Li_BrJAnaesth"},
    {"rsid": "rs138470454", "chr": "10", "pos": 14059742, "source": "Li_BrJAnaesth"},
    # --- Warner et al. 2017 (1 locus) ---
    {"rsid": "rs887797", "chr": "17", "pos": 66583327, "source": "Warner"},
    # --- Authors' own UK Biobank GWAS (4 loci) ---
    {"rsid": "rs151142747", "chr": "1", "pos": 93517560, "source": "own_GWAS"},
    {"rsid": "rs115024240", "chr": "1", "pos": None, "source": "own_GWAS"},  # excluded, LD with rs151142747
    {"rsid": "rs572547799", "chr": "6", "pos": 103805760, "source": "own_GWAS"},
    {"rsid": "rs35363701", "chr": "16", "pos": 83984133, "source": "own_GWAS"},
]
