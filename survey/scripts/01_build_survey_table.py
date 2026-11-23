"""Survey step 1 - one table with every peptide of the final non-redundant set and its class labels.

Peptides:  deliverables/fasta/amp_cdhit90.fasta (10,824 peptides, cleaned, X -> Ala, CD-HIT 90 %)
Labels:    ORIGIN / SYNTHESIS TYPE / KINGDOM / DOMAIN from by_origin/*.csv (the cleaned DBAASP columns),
           function (TARGET GROUP) from by_target_group/*.csv. No new parameter is calculated.
Kingdom is only known for Natural peptides: 8,315 of the 8,316 Synthetic peptides are 'Unclassified' in DBAASP.
Writes survey/tables/survey_peptides.csv and survey/tables/group_sizes.csv
"""
import pandas as pd
from common import ROOT, DELIV, TABLES, TARGET_GROUPS, MIN_N

TABLES.mkdir(parents=True, exist_ok=True)

seqs = {}
for line in (DELIV / "fasta" / "amp_cdhit90.fasta").read_text().splitlines():
    if line.startswith(">"):
        key = int(line[1:].replace("DBAASP_", ""))
    else:
        seqs[key] = line.strip()
df = pd.DataFrame({"ID": list(seqs), "SEQUENCE": list(seqs.values())})
df["LENGTH"] = df["SEQUENCE"].str.len()

meta = pd.concat([pd.read_csv(ROOT / "by_origin" / f"{o}.csv") for o in ("Natural", "Synthetic")])
meta = meta[["ID", "ORIGIN", "SYNTHESIS TYPE", "KINGDOM", "DOMAIN"]].drop_duplicates("ID")
df = df.merge(meta, on="ID", how="left", validate="1:1")
assert df["ORIGIN"].notna().all()

# kingdom: peptides with two kingdoms in DBAASP are called 'Multiple'; 'Unclassified' stays as it is
df["KINGDOM"] = df["KINGDOM"].where(~df["KINGDOM"].str.contains(",", na=False), "Multiple")
df["SOURCE_CLASS"] = df["SYNTHESIS TYPE"].replace({"Synthetic": "Synthetic"})   # Ribosomal / Nonribosomal / Synthetic

for key, label in TARGET_GROUPS.items():
    ids = set(pd.read_csv(ROOT / "by_target_group" / f"{key}.csv")["ID"])
    df[label] = df["ID"].isin(ids)
df["N_FUNCTIONS"] = df[list(TARGET_GROUPS.values())].sum(axis=1)

df.to_csv(TABLES / "survey_peptides.csv", index=False)

rows = []
nat = df[df["ORIGIN"] == "Natural"]
for k, g in nat.groupby("KINGDOM"):
    rows.append(("Source organism (Natural peptides only)", k, len(g)))
for k, g in df.groupby("SOURCE_CLASS"):
    rows.append(("Synthesis type (all peptides)", k, len(g)))
for label in TARGET_GROUPS.values():
    rows.append(("Function (all peptides)", label, int(df[label].sum())))
    rows.append(("Function (Natural peptides only)", label, int(nat[label].sum())))
sizes = pd.DataFrame(rows, columns=["grouping", "group", "n_peptides"])
sizes["tested"] = sizes["n_peptides"] >= MIN_N
sizes.to_csv(TABLES / "group_sizes.csv", index=False)
print(sizes.to_string(index=False))
print(f"\n{len(df):,} peptides; {int((df['N_FUNCTIONS'] == 0).sum())} without any function label; "
      f"median {df['N_FUNCTIONS'].median():.0f} functions per peptide")
