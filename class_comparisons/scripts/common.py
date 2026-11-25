"""Shared definitions for the class comparisons: paths, the five comparisons and their classes, plotting style.

Every comparison uses the 11,716 peptides of deliverables/ (cleaned, X -> Ala, CD-HIT 90 %, Natural + Synthetic clustered together),
so all figures are comparable with the Natural vs Synthetic figures in deliverables/figures/.
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CC = ROOT / "class_comparisons"
DELIV = ROOT / "deliverables"
RAW_ZIP = ROOT / "raw_data" / "dbaasp_datasets_corrected.zip"
OVERVIEW = CC / "overview"
BACKGROUNDS = {"SwissProt": DELIV / "background" / "swissprot51_5k.fasta",
               "Non-AMP": DELIV / "background" / "pepbench_nonamp_background.fasta"}

MIN_N = 50            # a class needs at least this many peptides to get figures and tests
ITERATIONS = 10000    # Composition Profiler permutations (tool default)
ALPHA = 0.05          # before Bonferroni over the 20 amino acids

# DBAASP TARGET GROUP token -> label (same tokens as the notebook's activity columns)
TARGET_GROUPS = {"Gram-": "Anti-Gram-negative", "Gram+": "Anti-Gram-positive", "Fungus": "Antifungal", "Cancer": "Anticancer",
                 "Virus": "Antiviral", "Biofilm": "Antibiofilm", "Parasite": "Antiparasitic", "Insect": "Insecticidal",
                 "Protista": "Antiprotozoal", "Mollicute": "Anti-Mollicute", "Nematode": "Nematicidal",
                 "Mammalian Cell": "Mammalian cell"}
TARGET_OBJECTS = ["Lipid Bilayer", "Virus entry", "Cytoplasmic Protein", "Virus Integrase", "DNA / RNA", "Virus replication",
                  "Virus Protease", "Virus Reverse Transcriptase", "gp41", "Membrane Protein", "Immunomodulation", "EgpDV", "gp120",
                  "VIF", "APJ receptor", "Sgp", "Ribosome", "gpE", "Quorum Sensing", "Lysosomal membrane",
                  "Mitochondrial membrane", "gpHA", "gD",
                  "Lipid II"]   # Lipid II is not in the notebook's list of 23 (so there is no by_target_object/Lipid_II.csv), added here
KINGDOMS = ["Animalia", "Bacteria", "Plantae", "Fungi", "Archaea", "Protista", "Virus", "Multiple", "Unclassified"]

COMPARISONS = {   # folder -> title, population, how classes are defined
    "source_organism": {"title": "Source organism (Natural peptides)", "population": "Natural",
                        "note": "Kingdom is only given for Natural peptides in DBAASP."},
    "synthesis_type": {"title": "Synthesis type", "population": "All",
                       "note": "Ribosomal + Nonribosomal = Natural."},
    "function_all": {"title": "Function: target group (all peptides)", "population": "All",
                     "note": "A peptide can belong to several classes; 'tested against', not 'potent against'."},
    "function_natural": {"title": "Function: target group (Natural peptides)", "population": "Natural",
                         "note": "A peptide can belong to several classes; 'tested against', not 'potent against'."},
    "target_object": {"title": "Target object (mechanism, all peptides)", "population": "All",
                      "note": "A peptide can belong to several classes."},
}


def slug(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def load_table():
    return pd.read_csv(OVERVIEW / "tables" / "class_peptides.csv")


def classes(df, key):
    """[(class name, boolean mask)] for one comparison, every class of its population (also the small ones), largest first."""
    pop = df["ORIGIN"] == "Natural" if COMPARISONS[key]["population"] == "Natural" else pd.Series(True, index=df.index)
    if key == "source_organism":
        out = [(k, pop & (df["KINGDOM_CLASS"] == k)) for k in KINGDOMS]
    elif key == "synthesis_type":
        out = [(k, df["SYNTHESIS TYPE"] == k) for k in ["Ribosomal", "Nonribosomal", "Synthetic"]]
    elif key.startswith("function"):
        out = [(lab, pop & df[f"TG_{lab}"]) for lab in TARGET_GROUPS.values()]
    else:
        out = [(o, pop & df[f"TO_{o}"]) for o in TARGET_OBJECTS]
    out = [(n, m) for n, m in out if m.sum() > 0]
    if key != "synthesis_type":
        out.sort(key=lambda c: -int(c[1].sum()))
    return out, pop


def tested(df, key):
    cl, pop = classes(df, key)
    return [(n, m) for n, m in cl if m.sum() >= MIN_N], pop


# ---------------------------------------------------------------- plotting style (same tokens as deliverables/scripts/02 and 05)
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = "#2a78d6"                                    # categorical slot 1 (one class per panel, so one hue is enough)
LENGTH_RAMP = ["#b7d3f6", "#6da7ec", "#2a78d6", "#104281"]   # sequential blue, Ultrashort -> Long
AA_ORDER = "RKDENQHPYWSTGAMCFLVI"                     # Kyte-Doolittle order of the Composition Profiler (hydrophilic -> hydrophobic)
CATEGORIES = [("Ultrashort", "3–9 aa", 3, 9), ("Short", "10–24 aa", 10, 24), ("Medium", "25–50 aa", 25, 50), ("Long", "51–100 aa", 51, 100)]
