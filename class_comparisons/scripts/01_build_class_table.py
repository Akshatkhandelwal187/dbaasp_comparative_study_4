"""Class comparisons step 1 - one table with every peptide of the final set and its class labels, plus class sizes,
per-class FASTA files and the overlap of the function / target-object classes.

Peptides  deliverables/tables/peptide_metadata.csv (11,716 peptides; sequences as in the FASTA: upper case, X -> Ala)
Labels    the existing DBAASP columns of raw_data/dbaasp_datasets_corrected.zip (by_origin/*.csv): SYNTHESIS TYPE, KINGDOM,
          TARGET GROUP, TARGET OBJECT. Nothing new is calculated.
Writes    overview/tables/class_peptides.csv     ID, sequence, length, origin, synthesis type, kingdom, one True/False column per
                                                 target group (TG_...) and target object (TO_...)
          overview/tables/class_sizes.csv        every class of every comparison, n, and whether it is tested (n >= MIN_N)
          overview/tables/function_overlap_jaccard.csv, target_object_overlap_jaccard.csv
          <comparison>/fasta/<class>.fasta       the peptides of every tested class (same headers and sequences as deliverables/fasta)
"""
import io
import re
import zipfile

import pandas as pd

from common import (COMPARISONS, CC, DELIV, MIN_N, OVERVIEW, RAW_ZIP, TARGET_GROUPS, TARGET_OBJECTS, classes, slug)

(OVERVIEW / "tables").mkdir(parents=True, exist_ok=True)
meta = pd.read_csv(DELIV / "tables" / "peptide_metadata.csv")
zf = zipfile.ZipFile(RAW_ZIP)
raw = pd.concat([pd.read_csv(io.BytesIO(zf.read(f"by_origin/{o}.csv"))) for o in ("Natural", "Synthetic")])
raw["ID"] = "DBAASP_" + raw["ID"].astype(str)
raw = raw.drop_duplicates("ID").set_index("ID")

df = meta.rename(columns={"LENGTH_CLASS": "LENGTH_CATEGORY"}).copy()
for col in ["SYNTHESIS TYPE", "KINGDOM", "TARGET GROUP", "TARGET OBJECT"]:
    df[col] = df["ID"].map(raw[col])
assert df["SYNTHESIS TYPE"].notna().all()
assert (df["ORIGIN"] == df["SYNTHESIS TYPE"].map({"Ribosomal": "Natural", "Nonribosomal": "Natural", "Synthetic": "Synthetic"})).all()
df["KINGDOM_CLASS"] = df["KINGDOM"].where(~df["KINGDOM"].str.contains(",", na=False), "Multiple")

tokens = df["TARGET GROUP"].map(lambda v: [] if pd.isna(v) else [p.strip() for p in re.split(r"\s{2,}", str(v).strip()) if p.strip()])
for tok, lab in TARGET_GROUPS.items():
    df[f"TG_{lab}"] = tokens.map(lambda t, x=tok: x in t)
pattern = "|".join(re.escape(o) for o in sorted(TARGET_OBJECTS, key=len, reverse=True))
objs = df["TARGET OBJECT"].map(lambda v: [] if pd.isna(v) else re.findall(pattern, str(v)))
for o in TARGET_OBJECTS:
    df[f"TO_{o}"] = objs.map(lambda t, x=o: x in t)
df["N_TARGET_GROUPS"] = df[[f"TG_{l}" for l in TARGET_GROUPS.values()]].sum(axis=1)
df.drop(columns=["TARGET GROUP", "TARGET OBJECT"]).to_csv(OVERVIEW / "tables" / "class_peptides.csv", index=False)

rows = []
for key, spec in COMPARISONS.items():
    cl, pop = classes(df, key)
    fasta_dir = CC / key / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    for f in fasta_dir.glob("*.fasta"):
        f.unlink()
    for name, mask in cl:
        n = int(mask.sum())
        rows.append({"comparison": key, "title": spec["title"], "class": name, "n_peptides": n,
                     "share_of_population_%": round(100 * n / int(pop.sum()), 1), "tested": n >= MIN_N})
        if n >= MIN_N:
            with open(fasta_dir / f"{slug(name)}.fasta", "w") as fh:
                for i, s in zip(df.loc[mask, "ID"], df.loc[mask, "SEQUENCE"]):
                    fh.write(f">{i}\n{s}\n")
sizes = pd.DataFrame(rows)
sizes.to_csv(OVERVIEW / "tables" / "class_sizes.csv", index=False)


def jaccard(cols, names):
    j = pd.DataFrame(index=names, columns=names, dtype=float)
    for a, ca in zip(names, cols):
        for b, cb in zip(names, cols):
            j.loc[a, b] = round((df[ca] & df[cb]).sum() / (df[ca] | df[cb]).sum(), 3)
    return j


fn = [l for l in TARGET_GROUPS.values() if df[f"TG_{l}"].sum() >= MIN_N]
jaccard([f"TG_{l}" for l in fn], fn).to_csv(OVERVIEW / "tables" / "function_overlap_jaccard.csv")
to = [o for o in TARGET_OBJECTS if df[f"TO_{o}"].sum() >= MIN_N]
jaccard([f"TO_{o}" for o in to], to).to_csv(OVERVIEW / "tables" / "target_object_overlap_jaccard.csv")

print(sizes.to_string(index=False))
print(f"\n{len(df):,} peptides; {int((df['N_TARGET_GROUPS'] == 0).sum())} without a target group; "
      f"median {df['N_TARGET_GROUPS'].median():.0f} target groups per peptide; "
      f"{int((objs.map(len) == 0).sum())} without a target object")
