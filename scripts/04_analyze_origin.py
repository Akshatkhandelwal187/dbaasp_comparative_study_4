"""Step 4 - Natural vs Synthetic, using only columns that already exist in the data.
Input : nonredundant_csv/by_origin/Natural.csv and Synthetic.csv  (after cleaning + CD-HIT 90 %)
Output: tables (csv) and plots (png) in analysis/origin/
"""
import ast
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

OUT = "analysis/origin"
nat = pd.read_csv("nonredundant_csv/by_origin/Natural.csv")
syn = pd.read_csv("nonredundant_csv/by_origin/Synthetic.csv")
data = {"Natural": nat, "Synthetic": syn}
tests = []          # every statistical test is collected here and saved at the end


def share_table(counts_nat, counts_syn):
    """% of peptides in each origin that have every category."""
    t = pd.DataFrame({"Natural_%": counts_nat / len(nat) * 100, "Synthetic_%": counts_syn / len(syn) * 100})
    t["Natural_n"] = counts_nat
    t["Synthetic_n"] = counts_syn
    return t.fillna(0).round(1)


# ---------------------------------------------------------------- 1. length
lengths = pd.DataFrame({name: d["Length_of_Sequence"].describe() for name, d in data.items()}).round(1)
lengths.to_csv(f"{OUT}/length_summary.csv")
print(lengths, "\n")

u = stats.mannwhitneyu(nat["Length_of_Sequence"], syn["Length_of_Sequence"])
tests.append({"test": "Length: Natural vs Synthetic (Mann-Whitney U)", "p_value": u.pvalue})

# histogram (share of peptides, so the two different sizes can be compared) and boxplot
plt.figure(figsize=(7, 4))
bins = range(0, 101, 2)
for name, d in data.items():
    plt.hist(d["Length_of_Sequence"], bins=bins, alpha=0.5, density=True, label=f"{name} (n={len(d)})")
plt.xlabel("Peptide length (aa)")
plt.ylabel("Density")
plt.title("Length distribution (x axis cut at 100 aa)")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/length_histogram.png", dpi=150)
plt.close()

plt.figure(figsize=(5, 4))
plt.boxplot([nat["Length_of_Sequence"], syn["Length_of_Sequence"]], tick_labels=["Natural", "Synthetic"])
plt.yscale("log")
plt.ylabel("Peptide length (aa, log scale)")
plt.title("Length boxplot")
plt.tight_layout()
plt.savefig(f"{OUT}/length_boxplot.png", dpi=150)
plt.close()

# ------------------------------------------------------------- 2. outliers
# outlier rule of a boxplot: longer than Q3 + 1.5 * IQR (or shorter than Q1 - 1.5 * IQR)
outlier_rows = []
for name, d in data.items():
    q1, q3 = d["Length_of_Sequence"].quantile([0.25, 0.75])
    low, high = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
    out = d[(d["Length_of_Sequence"] < low) | (d["Length_of_Sequence"] > high)]
    print(f"{name}: outlier fences {low:.1f} - {high:.1f} aa, {len(out)} outliers ({len(out) / len(d) * 100:.1f} %)")
    outlier_rows.append(out.assign(group=name, fence_low=low, fence_high=high)
                        [["group", "ID", "NAME", "Length_of_Sequence", "fence_low", "fence_high"]])
outliers = pd.concat(outlier_rows).sort_values(["group", "Length_of_Sequence"], ascending=[True, False])
outliers.to_csv(f"{OUT}/length_outliers.csv", index=False)
print(outliers.groupby("group").head(5).to_string(index=False), "\n")

# ------------------------------------------------------ 3. length classes
order = ["Ultrashort (2–9 aa)", "Short (10–24 aa)", "Medium (25–50 aa)", "Long (50–100 aa)",
         "Antimicrobial Protein (>100 aa)", "Out of range (<2 aa)"]
cls = share_table(nat["Length_Class"].value_counts(), syn["Length_Class"].value_counts()).reindex(order)
cls.to_csv(f"{OUT}/length_class_share.csv")
print(cls, "\n")
p = stats.chi2_contingency(cls[["Natural_n", "Synthetic_n"]].values)[1]
tests.append({"test": "Length class mix: Natural vs Synthetic (chi-square)", "p_value": p})
cls[["Natural_%", "Synthetic_%"]].plot.bar(figsize=(8, 4), rot=30)
plt.ylabel("% of peptides")
plt.title("Length class of Natural vs Synthetic peptides")
plt.tight_layout()
plt.savefig(f"{OUT}/length_class_share.png", dpi=150)
plt.close()

# ------------------------------------------------- 4. synthesis type (existing column)
synth = pd.concat([nat, syn]).groupby("SYNTHESIS TYPE")["Length_of_Sequence"].describe().round(1)
synth.to_csv(f"{OUT}/synthesis_type_length.csv")
print(synth, "\n")

# -------------------------------------- 5. terminal modifications, D-amino acids, X
# a terminus is "modified" if the N TERMINUS / C TERMINUS column is filled in.
# D-amino acids are the lowercase letters in SEQUENCE, X is a non-standard amino acid.
flags = pd.DataFrame(index=["N-terminus modified", "C-terminus modified",
                            "contains D-amino acid (lowercase)", "contains X (non-standard residue)"])
for name, d in data.items():
    flags[f"{name}_%"] = [d["N TERMINUS"].notna().mean() * 100, d["C TERMINUS"].notna().mean() * 100,
                          d["SEQUENCE"].str.contains("[a-z]").mean() * 100, d["SEQUENCE"].str.upper().str.contains("X").mean() * 100]
flags = flags.round(1)
flags.to_csv(f"{OUT}/terminus_D_aa_X_share.csv")
print(flags, "\n")
for label, column in [("N-terminus modified", "N TERMINUS"), ("C-terminus modified", "C TERMINUS")]:
    table = [[nat[column].notna().sum(), nat[column].isna().sum()], [syn[column].notna().sum(), syn[column].isna().sum()]]
    tests.append({"test": f"{label}: Natural vs Synthetic (chi-square)", "p_value": stats.chi2_contingency(table)[1]})
for label, pattern, upper in [("Contains D-amino acid", "[a-z]", False), ("Contains X", "X", True)]:
    a = nat["SEQUENCE"].str.upper() if upper else nat["SEQUENCE"]
    b = syn["SEQUENCE"].str.upper() if upper else syn["SEQUENCE"]
    table = [[a.str.contains(pattern).sum(), (~a.str.contains(pattern)).sum()],
             [b.str.contains(pattern).sum(), (~b.str.contains(pattern)).sum()]]
    tests.append({"test": f"{label}: Natural vs Synthetic (chi-square)", "p_value": stats.chi2_contingency(table)[1]})

# most common modifications
mods = {}
for column in ["N TERMINUS", "C TERMINUS"]:
    for name, d in data.items():
        mods[f"{column} - {name}"] = d[column].str.strip().value_counts().head(5)
pd.concat(mods, axis=1).to_csv(f"{OUT}/top_terminus_modifications.csv")

# --------------------------------------------- 6. target group and target object
# one peptide can have several targets, so the columns are split into single labels
def count_labels(d, column):
    labels = d[column].apply(ast.literal_eval).explode().dropna()
    return labels.value_counts()

for title, column, top in [("target_group", "tokenised_TARGET GROUP", 20), ("target_object", "tokenised_TARGET OBJECT", 8)]:
    t = share_table(count_labels(nat, column), count_labels(syn, column))
    t = t.sort_values("Synthetic_n", ascending=False).head(top)
    # chi-square for each label: has label / does not, Natural vs Synthetic
    t["p_value"] = [stats.chi2_contingency([[a, len(nat) - a], [b, len(syn) - b]])[1] for a, b in zip(t["Natural_n"], t["Synthetic_n"])]
    t.to_csv(f"{OUT}/{title}_share.csv")
    print(t, "\n")
    if title == "target_group":
        t[["Natural_%", "Synthetic_%"]].iloc[::-1].plot.barh(figsize=(7, 5))
        plt.xlabel("% of peptides with this target")
        plt.title("Target group of Natural vs Synthetic peptides")
        plt.tight_layout()
        plt.savefig(f"{OUT}/target_group_share.png", dpi=150)
        plt.close()

# ----------------------------------------- 7. kingdom (only informative for Natural)
kingdom = share_table(nat["KINGDOM"].value_counts(), syn["KINGDOM"].value_counts()).sort_values("Natural_n", ascending=False)
kingdom.to_csv(f"{OUT}/kingdom_share.csv")
print(kingdom, "\n")

pd.DataFrame(tests).to_csv(f"{OUT}/statistical_tests.csv", index=False)
print(pd.DataFrame(tests).to_string(index=False))
