"""Step 5 - compare the length classes with each other, using only columns that already exist.
Input : nonredundant_csv/by_length_class/*.csv  (after cleaning + CD-HIT 90 %)
Output: tables (csv) and plots (png) in analysis/length_class/
"""
import ast
import glob
import pandas as pd
import matplotlib.pyplot as plt

OUT = "analysis/length_class"
order = ["Ultrashort (2–9 aa)", "Short (10–24 aa)", "Medium (25–50 aa)", "Long (50–100 aa)",
         "Antimicrobial Protein (>100 aa)", "Out of range (<2 aa)"]

data = pd.concat([pd.read_csv(f) for f in glob.glob("nonredundant_csv/by_length_class/*.csv")], ignore_index=True)
by_class = data.groupby("Length_Class")

# ------------------------------------------- 1. how big is each class, what lengths
summary = by_class["Length_of_Sequence"].agg(["count", "min", "median", "mean", "max"]).round(1).reindex(order)
summary["Natural_%"] = by_class["ORIGIN"].apply(lambda s: (s == "Natural").mean() * 100).round(1)
summary["Synthetic_%"] = 100 - summary["Natural_%"]
summary.to_csv(f"{OUT}/class_summary.csv")
print(summary, "\n")

# all lengths in one histogram with the class borders
plt.figure(figsize=(8, 4))
plt.hist(data["Length_of_Sequence"], bins=range(0, 105, 1), color="grey")
for border in [9.5, 24.5, 50.5]:
    plt.axvline(border, color="red", linestyle="--")
plt.xlabel("Peptide length (aa, cut at 100)")
plt.ylabel("Number of peptides")
plt.title("All non-redundant peptides; red lines = class borders")
plt.tight_layout()
plt.savefig(f"{OUT}/length_histogram_all.png", dpi=150)
plt.close()

# ----------------------------------------------------- 2. natural vs synthetic per class
origin = pd.crosstab(data["Length_Class"], data["ORIGIN"], normalize="index").reindex(order) * 100
origin.plot.bar(stacked=True, figsize=(8, 4), rot=30)
plt.ylabel("% of peptides in the class")
plt.title("Natural / Synthetic share in every length class")
plt.tight_layout()
plt.savefig(f"{OUT}/origin_share_by_class.png", dpi=150)
plt.close()

# ------------------------------- 3. terminal modifications, D-amino acids, X residues
flags = pd.DataFrame({
    "N-terminus modified %": by_class["N TERMINUS"].apply(lambda s: s.notna().mean() * 100),
    "C-terminus modified %": by_class["C TERMINUS"].apply(lambda s: s.notna().mean() * 100),
    "D-amino acid %": by_class["SEQUENCE"].apply(lambda s: s.str.contains("[a-z]").mean() * 100),
    "contains X %": by_class["SEQUENCE"].apply(lambda s: s.str.upper().str.contains("X").mean() * 100),
}).round(1).reindex(order)
flags.to_csv(f"{OUT}/terminus_D_aa_X_by_class.csv")
print(flags, "\n")

# ----------------------------------------------- 4. target group per class (heat map)
rows = []
for name, d in by_class:
    labels = d["tokenised_TARGET GROUP"].apply(ast.literal_eval).explode().dropna()
    rows.append((labels.value_counts() / len(d) * 100).rename(name))
targets = pd.DataFrame(rows).fillna(0).reindex(order).round(1)
targets = targets[targets.mean().sort_values(ascending=False).index]      # most common targets first
targets.to_csv(f"{OUT}/target_group_by_class.csv")
print(targets, "\n")

plt.figure(figsize=(9, 4))
plt.imshow(targets.values, cmap="Blues", aspect="auto")
plt.xticks(range(targets.shape[1]), targets.columns, rotation=45, ha="right")
plt.yticks(range(targets.shape[0]), targets.index)
for i in range(targets.shape[0]):
    for j in range(targets.shape[1]):
        plt.text(j, i, f"{targets.values[i, j]:.0f}", ha="center", va="center", fontsize=7)
plt.colorbar(label="% of peptides in the class")
plt.title("Target group by length class")
plt.tight_layout()
plt.savefig(f"{OUT}/target_group_heatmap.png", dpi=150)
plt.close()
