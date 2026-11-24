"""Survey step 3 - length distribution and function overlap for every class (run 01 first).

Per class: n, median, quartiles, mean, share of peptides per length category (the categories of the main analysis:
Ultrashort 2-9, Short 10-24, Medium 25-50, Long 51-100 aa; peptides here are 3-100 aa).
Test: Mann-Whitney U of the class against all OTHER peptides of the same population (effect = rank-biserial r, in -1..1,
positive = class is longer). Kruskal-Wallis across the source-organism classes and across the synthesis types.
Also the overlap of the function classes (Jaccard index), because one peptide can have several functions.
Writes survey/tables/survey_length_summary.csv, survey_length_tests.csv, function_overlap_jaccard.csv
"""
import numpy as np
import pandas as pd
from scipy import stats
from common import TABLES, TARGET_GROUPS, MIN_N

df = pd.read_csv(TABLES / "survey_peptides.csv")
nat = df["ORIGIN"] == "Natural"
bins = [0, 9, 24, 50, 100]
labels = ["Ultrashort (3-9)", "Short (10-24)", "Medium (25-50)", "Long (51-100)"]
df["LENGTH_CATEGORY"] = pd.cut(df["LENGTH"], bins=bins, labels=labels)

groups = []
for k in ["Animalia", "Bacteria", "Plantae", "Fungi"]:
    groups.append(("Source organism (Natural)", k, nat & (df["KINGDOM"] == k), nat))
for k in ["Ribosomal", "Nonribosomal", "Synthetic"]:
    groups.append(("Synthesis type", k, df["SOURCE_CLASS"] == k, pd.Series(True, index=df.index)))
for label in TARGET_GROUPS.values():
    groups.append(("Function (all)", label, df[label], pd.Series(True, index=df.index)))
    groups.append(("Function (Natural)", label, nat & df[label], nat))

summ, tests = [], []
for grouping, name, mask, population in groups:
    n = int(mask.sum())
    if n < MIN_N:
        continue
    L = df.loc[mask, "LENGTH"]
    rest = df.loc[population & ~mask, "LENGTH"]
    row = {"grouping": grouping, "group": name, "n_peptides": n, "median": L.median(), "q1": L.quantile(.25),
           "q3": L.quantile(.75), "mean": round(L.mean(), 2)}
    shares = df.loc[mask, "LENGTH_CATEGORY"].value_counts(normalize=True).reindex(labels) * 100
    row.update({f"%{lab}": round(v, 1) for lab, v in shares.items()})
    summ.append(row)
    if len(rest) >= MIN_N:
        u = stats.mannwhitneyu(L, rest, alternative="two-sided")
        tests.append({"grouping": grouping, "group": name, "n_group": n, "n_rest": len(rest),
                      "median_group": L.median(), "median_rest": rest.median(),
                      "rank_biserial_r": round(2 * u.statistic / (n * len(rest)) - 1, 3), "pvalue": u.pvalue})
summ = pd.DataFrame(summ)
tests = pd.DataFrame(tests)
# Bonferroni over all tests in this table
tests["pvalue_bonferroni"] = np.minimum(tests["pvalue"] * len(tests), 1)
summ.to_csv(TABLES / "survey_length_summary.csv", index=False)
tests.to_csv(TABLES / "survey_length_tests.csv", index=False)

for title, parts in [("source organism (Natural)", [df.loc[nat & (df["KINGDOM"] == k), "LENGTH"] for k in ["Animalia", "Bacteria", "Plantae", "Fungi"]]),
                     ("synthesis type", [df.loc[df["SOURCE_CLASS"] == k, "LENGTH"] for k in ["Ribosomal", "Nonribosomal", "Synthetic"]])]:
    h = stats.kruskal(*parts)
    print(f"Kruskal-Wallis, length across {title}: H = {h.statistic:.1f}, p = {h.pvalue:.2e}")

fn = list(TARGET_GROUPS.values())
tested = [f for f in fn if df[f].sum() >= MIN_N]
jac = pd.DataFrame(index=tested, columns=tested, dtype=float)
for a in tested:
    for b in tested:
        inter, union = (df[a] & df[b]).sum(), (df[a] | df[b]).sum()
        jac.loc[a, b] = round(inter / union, 3)
jac.to_csv(TABLES / "function_overlap_jaccard.csv")
print(summ.to_string(index=False))
print(tests.assign(pvalue=tests["pvalue"].map("{:.1e}".format), pvalue_bonferroni=tests["pvalue_bonferroni"].map("{:.1e}".format)).to_string(index=False))
