"""Rebuild raw_data/dbaasp_datasets_corrected.zip from the DBAASP export pages (script version of dbaasp_data_final.ipynb).

Reads   data collection/combined_df/peptides *.csv   13 export pages of the DBAASP peptide list
        data collection/kingdom/peptides *.csv       DBAASP peptide lists per kingdom
Writes  raw_data/dbaasp_datasets_corrected.zip       by_origin/ by_length_class/ by_domain/ by_kingdom/ by_target_group/
                                                     by_target_object/ + manifest.csv   (same layout and columns as the notebook)
        manifest.csv                                 rows per subset (repository root)

The steps are the notebook's, in the same order and with the same labels, so the subsets keep their columns and file names:
  * ORIGIN: Ribosomal / Nonribosomal -> Natural, Synthetic -> Synthetic
  * Length_of_Sequence = number of characters of SEQUENCE; Length_Class after the APD classification
  * tokenised TARGET GROUP / TARGET OBJECT (stored as text lists)
  * KINGDOM by matching the ID with the per-kingdom lists (several kingdoms -> joined with ', '; none -> 'Unclassified'); DOMAIN from KINGDOM
  * 12 activity columns from TARGET GROUP (added before the target-group and target-object subsets are written, as in the notebook;
    one difference: an empty TARGET GROUP counts as 0 activities here, the notebook counted it as one activity called "nan")
Usage:  python "data collection/build_subsets.py"  [--out <zip>]
"""
import argparse
import io
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

PAGES = ["peptides 0-2000.csv", "peptides 2000-4000.csv", "peptides 4000-6000.csv", "peptides 6000-8000.csv",
         "peptides 8000-10000.csv", "peptides 10000-12000.csv", "peptides 12000-14000.csv", "peptides 14000-16000.csv",
         "peptides 16000-18000.csv", "peptides 18000-20000.csv", "peptides 20000-22000.csv", "peptides 22000-24000.csv",
         "peptides 24000-25542.csv"]
KINGDOM_FILES = {   # same order as the notebook (decides the order of joined labels such as "Bacteria, Fungi")
    "Animalia": ["peptides Animalia 0-2000.csv", "peptides Animalia 2000-2645.csv"],
    "Archaea": ["peptides Archaea 0-10.csv"],
    "Bacteria": ["peptides bacteria 0-789.csv"],
    "Fungi": ["peptides fungi 0-307.csv"],
    "Protista": ["peptides protista 0-2.csv"],
    "Virus": ["peptides virus 0-2.csv"],
    "Plantae": ["peptides plantae 0-352.csv"],
}
DOMAIN_MAP = {"Bacteria": "Bacteria", "Archaea": "Archaea", "Animalia": "Eukarya", "Plantae": "Eukarya", "Fungi": "Eukarya",
              "Protista": "Eukarya", "Virus": "N/A (acellular; outside the three-domain system)"}
KNOWN_OBJECTS = ["Lipid Bilayer", "Virus entry", "Cytoplasmic Protein", "Virus Integrase", "DNA / RNA", "Virus replication",
                 "Virus Protease", "Virus Reverse Transcriptase", "gp41", "Membrane Protein", "Immunomodulation", "EgpDV", "gp120",
                 "VIF", "APJ receptor", "Sgp", "Ribosome", "gpE", "Quorum Sensing", "Lysosomal membrane",
                 "Mitochondrial membrane", "gpHA", "gD"]
ACTIVITY_LABELS = {"Gram+": "antibacterial_gram_pos", "Gram-": "antibacterial_gram_neg", "Fungus": "antifungal",
                   "Virus": "antiviral", "Cancer": "anticancer", "Biofilm": "antibiofilm", "Parasite": "antiparasitic",
                   "Insect": "insecticidal", "Protista": "antiprotozoal", "Mollicute": "antimollicute", "Nematode": "nematicidal",
                   "Mammalian Cell": "mammalian_cytotoxic"}

parser = argparse.ArgumentParser()
parser.add_argument("--out", default=str(ROOT / "raw_data" / "dbaasp_datasets_corrected.zip"))
args = parser.parse_args()

df = pd.concat([pd.read_csv(HERE / "combined_df" / p) for p in PAGES], ignore_index=True)
print(f"{len(df):,} rows, {df['ID'].nunique():,} distinct IDs, {int(df.duplicated().sum()):,} exact duplicate rows")

# tokenised target group (text list, as the notebook stores it)
tok = df["TARGET GROUP"].fillna("").str.replace("Mammalian Cell", "Mammalian_Cell", regex=False).str.split()
df["tokenised_TARGET GROUP"] = tok.apply(lambda items: [i.replace("_", " ") for i in items]).astype(str)

pattern = "|".join(re.escape(o) for o in sorted(KNOWN_OBJECTS, key=len, reverse=True))


def parse_target_objects(val):
    return [] if pd.isna(val) else re.findall(pattern, str(val))


df["tokenised_TARGET OBJECT"] = df["TARGET OBJECT"].apply(parse_target_objects).astype(str)
df["ORIGIN"] = df["SYNTHESIS TYPE"].map({"Ribosomal": "Natural", "Nonribosomal": "Natural", "Synthetic": "Synthetic"})
df["Length_of_Sequence"] = df["SEQUENCE"].str.len()
L = df["Length_of_Sequence"]
df["Length_Class"] = np.select([(L >= 2) & (L <= 9), (L >= 10) & (L <= 24), (L >= 25) & (L <= 50), (L > 50) & (L <= 100), L > 100],
                               ["Ultrashort (2–9 aa)", "Short (10–24 aa)", "Medium (25–50 aa)", "Long (50–100 aa)",
                                "Antimicrobial Protein (>100 aa)"], default="Out of range (<2 aa)")

id_to_kingdoms = {}
for kingdom, files in KINGDOM_FILES.items():
    for f in files:
        for pid in pd.read_csv(HERE / "kingdom" / f)["ID"]:
            id_to_kingdoms.setdefault(pid, []).append(kingdom)
df["KINGDOM"] = df["ID"].map(lambda pid: ", ".join(id_to_kingdoms.get(pid, [])))
df.loc[df["KINGDOM"] == "", "KINGDOM"] = "Unclassified"
df["DOMAIN"] = df["KINGDOM"].map(DOMAIN_MAP)

# ---------------------------------------------------------------- subsets
files, manifest = {}, []


def slugify(value):
    value = re.sub(r"[\\/]+", "-", str(value).strip())
    value = re.sub(r"\s+", "_", value)
    return re.sub(r"[^\w\-]", "", value)


def write_subset(subset, dimension, label):
    name = f"{dimension}/{slugify(label)}.csv"
    files[name] = subset.to_csv(index=False)
    manifest.append({"dimension": dimension, "value": label, "n_rows": len(subset), "file": name})


for col, dim in [("DOMAIN", "by_domain"), ("ORIGIN", "by_origin"), ("Length_Class", "by_length_class")]:
    for value, group in df.groupby(col):
        write_subset(group, dim, value)
for k in ["Animalia", "Archaea", "Bacteria", "Fungi", "Protista", "Virus", "Plantae"]:
    write_subset(df[df["KINGDOM"].str.contains(rf"\b{k}\b", na=False)], "by_kingdom", k)
write_subset(df[df["KINGDOM"].isin(["Unclassified"])], "by_kingdom", "Unresolved")

tokens = df["TARGET GROUP"].map(lambda v: [] if pd.isna(v) or str(v).strip() == "" else
                                [p.strip() for p in re.split(r"\s{2,}", str(v).strip()) if p.strip()])
for raw, label in ACTIVITY_LABELS.items():
    df[label] = tokens.map(lambda t, tok_=raw: tok_ in t)
df["n_activities"] = tokens.map(len)
df["activity_list"] = tokens.map(lambda t: "; ".join(ACTIVITY_LABELS.get(x, x) for x in t))
for label in ACTIVITY_LABELS.values():
    write_subset(df[df[label]], "by_target_group", label)
objs = df["TARGET OBJECT"].apply(parse_target_objects)
for obj in KNOWN_OBJECTS:
    write_subset(df[objs.apply(lambda lst, o=obj: o in lst)], "by_target_object", obj)

man = pd.DataFrame(manifest)
with zipfile.ZipFile(args.out, "w", zipfile.ZIP_DEFLATED) as z:
    for name, text in files.items():
        z.writestr(name, text)
    z.writestr("manifest.csv", man.to_csv(index=False))
man.to_csv(ROOT / "manifest.csv", index=False)
print(f"wrote {len(files)} subsets to {args.out}")
print(man[man["dimension"].isin(["by_origin", "by_kingdom"])].to_string(index=False))
