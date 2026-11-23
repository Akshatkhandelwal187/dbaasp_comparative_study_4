"""Shared definitions for the survey scripts (paths, group definitions, minimum group size)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SURVEY = ROOT / "survey"
DELIV = ROOT / "deliverables"
TABLES = SURVEY / "tables"
FIGS = SURVEY / "figures"

MIN_N = 50          # a group needs at least this many non-redundant peptides to be tested
ITERATIONS = 10000  # permutations, same as the Composition Profiler default
ALPHA = 0.05 / 20   # Bonferroni over the 20 amino acids, same as deliverables/04_aa_enrichment.py
SEED = 1

TARGET_GROUPS = {   # file name in by_target_group/  ->  label used in tables and figures
    "antibacterial_gram_neg": "Anti-Gram-negative",
    "antibacterial_gram_pos": "Anti-Gram-positive",
    "antifungal": "Antifungal",
    "anticancer": "Anticancer",
    "antiviral": "Antiviral",
    "antibiofilm": "Antibiofilm",
    "antiparasitic": "Antiparasitic",
    "insecticidal": "Insecticidal",
    "antiprotozoal": "Antiprotozoal",
    "antimollicute": "Anti-Mollicute",
    "nematicidal": "Nematicidal",
    "mammalian_cytotoxic": "Mammalian-cytotoxic",
}
# the 20 amino acids ordered hydrophilic -> hydrophobic, same order as the Natural/Synthetic figure
AA_ORDER = list("RKDENQHSTGPACYMVWILF")
