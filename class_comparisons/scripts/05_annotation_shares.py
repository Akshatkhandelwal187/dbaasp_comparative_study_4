"""Class comparisons step 5 - other DBAASP annotations per class (run 01 first).

For Natural vs Synthetic and for every tested class of the non-overlapping comparisons (source organism, synthesis type),
the share of peptides with a modified C-terminus, a modified N-terminus, at least one D-amino acid, a former X
(changed to Ala in the cleaning), and the share tested against every target group. Natural vs Synthetic differences are
tested with the chi-square test (2 x 2).
Columns used: C TERMINUS, N TERMINUS, SEQUENCE (lower case = D-amino acid, X = non-standard) of the raw subsets.
Writes overview/tables/annotation_shares.csv, overview/tables/annotation_tests_natural_vs_synthetic.csv and
       overview/tables/former_x_share_by_class.csv (share of residues that were X, i.e. how much Ala is inflated in each class)
"""
import io
import zipfile

import pandas as pd
from scipy import stats

from common import COMPARISONS, OVERVIEW, RAW_ZIP, TARGET_GROUPS, load_table, tested

df = load_table()
zf = zipfile.ZipFile(RAW_ZIP)
raw = pd.concat([pd.read_csv(io.BytesIO(zf.read(f"by_origin/{o}.csv"))) for o in ("Natural", "Synthetic")])
raw["ID"] = "DBAASP_" + raw["ID"].astype(str)
raw = raw.drop_duplicates("ID").set_index("ID")
seq = df["ID"].map(raw["SEQUENCE"]).astype(str)
props = {"C-terminus modified": df["ID"].map(raw["C TERMINUS"]).notna(),
         "N-terminus modified": df["ID"].map(raw["N TERMINUS"]).notna(),
         "contains a D-amino acid": seq.str.contains(r"[a-wyz]"),
         "contained X (now Ala)": seq.str.upper().str.contains("X")}
props.update({f"tested: {lab}": df[f"TG_{lab}"] for lab in TARGET_GROUPS.values()})
props["target object: Lipid Bilayer"] = df["TO_Lipid Bilayer"]

groups = [("Origin", "Natural", df["ORIGIN"] == "Natural"), ("Origin", "Synthetic", df["ORIGIN"] == "Synthetic")]
for key in ["synthesis_type", "source_organism"]:
    cl, _ = tested(df, key)
    groups += [(key, n, m) for n, m in cl]
rows = []
for comp, name, mask in groups:
    row = {"comparison": comp, "class": name, "n_peptides": int(mask.sum())}
    row.update({p: round(100 * v[mask].mean(), 1) for p, v in props.items()})
    rows.append(row)
shares = pd.DataFrame(rows)
shares.to_csv(OVERVIEW / "tables" / "annotation_shares.csv", index=False)

nat = df["ORIGIN"] == "Natural"
tests = []
for p, v in props.items():
    chi2, pval, _, _ = stats.chi2_contingency(pd.crosstab(nat, v))
    tests.append({"property": p, "Natural_%": round(100 * v[nat].mean(), 1), "Synthetic_%": round(100 * v[~nat].mean(), 1),
                  "chi2": round(chi2, 1), "p_value": pval})
pd.DataFrame(tests).to_csv(OVERVIEW / "tables" / "annotation_tests_natural_vs_synthetic.csv", index=False)
print(shares.T.to_string())
print(pd.DataFrame(tests).to_string(index=False))

# share of residues that were X (now Ala) in every tested class: Ala enrichment in classes with a high share is an artefact
nx = seq.str.upper().str.count("X")
xs = []
for key in COMPARISONS:
    cl, _ = tested(df, key)
    for name, mask in cl:
        xs.append({"comparison": key, "class": name, "n_peptides": int(mask.sum()),
                   "peptides_with_former_X_%": round(100 * (nx[mask] > 0).mean(), 1),
                   "former_X_%_of_residues": round(100 * nx[mask].sum() / df.loc[mask, "LENGTH"].sum(), 2)})
pd.DataFrame(xs).to_csv(OVERVIEW / "tables" / "former_x_share_by_class.csv", index=False)
print(pd.DataFrame(xs).sort_values("former_X_%_of_residues", ascending=False).head(6).to_string(index=False))
