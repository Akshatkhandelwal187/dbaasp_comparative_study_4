"""Survey step 4 - figures from the tables of steps 1-3.

enrichment_<grouping>.png   heatmaps: 20 amino acids (hydrophilic -> hydrophobic) x classes, one panel per background,
                            colour = % change vs background, '*' = significant (p < 0.05/20), grey cell = not testable
length_boxplots.png         length distribution per class
length_category_share.png   share of peptides per length category per class
function_overlap.png        Jaccard overlap of the function classes
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from common import TABLES, FIGS, AA_ORDER, TARGET_GROUPS, MIN_N

FIGS.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
enr = pd.read_csv(TABLES / "survey_aa_enrichment.csv")
summ = pd.read_csv(TABLES / "survey_length_summary.csv")
df = pd.read_csv(TABLES / "survey_peptides.csv")

PANELS = [
    ("source_and_synthesis", "Source organism (Natural peptides) and synthesis type",
     [("Source organism (Natural)", g) for g in ["Animalia", "Bacteria", "Plantae", "Fungi"]]
     + [("Synthesis type", g) for g in ["Ribosomal", "Nonribosomal", "Synthetic"]]),
    ("function_all", "Function: all peptides (a peptide can belong to several classes)",
     [("Function (all)", g) for g in TARGET_GROUPS.values()]),
    ("function_natural", "Function: Natural peptides only",
     [("Function (Natural)", g) for g in TARGET_GROUPS.values()]),
]
LIM = 400   # colour scale limit in % (values beyond are clipped in colour only; numbers stay in the csv)

for fname, title, cols in PANELS:
    cols = [c for c in cols if ((enr["grouping"] == c[0]) & (enr["group"] == c[1])).any()]
    fig, axes = plt.subplots(1, 2, figsize=(1.0 + 0.62 * len(cols) * 2, 7.2), sharey=True)
    for ax, bg in zip(axes, ["SwissProt", "Non-AMP"]):
        mat = np.full((20, len(cols)), np.nan); sig = np.zeros((20, len(cols)), bool)
        for j, (gp, g) in enumerate(cols):
            sub = enr[(enr["grouping"] == gp) & (enr["group"] == g) & (enr["background"] == bg)].set_index("amino_acid")
            for i, a in enumerate(AA_ORDER):
                mat[i, j] = sub.loc[a, "pct_change"]; sig[i, j] = sub.loc[a, "result"] != "Not significant"
        im = ax.imshow(np.clip(mat, -LIM, LIM), cmap="RdBu_r", norm=TwoSlopeNorm(0, -LIM, LIM), aspect="auto")
        for i in range(20):
            for j in range(len(cols)):
                if sig[i, j]:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=13, color="black")
        n = [int(enr[(enr["grouping"] == gp) & (enr["group"] == g)]["n_peptides"].iloc[0]) for gp, g in cols]
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([f"{g}\n(n={k:,})" for (_, g), k in zip(cols, n)], rotation=60, ha="right", fontsize=8.5)
        ax.set_yticks(range(20)); ax.set_yticklabels(AA_ORDER)
        ax.set_title(f"vs {'SwissProt proteins' if bg == 'SwissProt' else 'non-antimicrobial peptides'}", fontsize=10.5)
        for s in ax.spines.values(): s.set_visible(False)
    fig.colorbar(im, ax=axes, shrink=0.6, pad=0.02, label=f"% change vs background (colour clipped at ±{LIM}%)")
    fig.suptitle(title + "\n* = significant (permutation test, Bonferroni alpha = 0.05/20)", fontsize=11)
    fig.savefig(FIGS / f"enrichment_{fname}.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# length distributions
order = [("Source organism (Natural)", g) for g in ["Animalia", "Bacteria", "Plantae", "Fungi"]] + \
        [("Synthesis type", g) for g in ["Ribosomal", "Nonribosomal", "Synthetic"]] + \
        [("Function (all)", g) for g in TARGET_GROUPS.values() if (summ["group"] == g).any()]
masks = {("Source organism (Natural)", g): (df["ORIGIN"] == "Natural") & (df["KINGDOM"] == g) for g in ["Animalia", "Bacteria", "Plantae", "Fungi"]}
masks.update({("Synthesis type", g): df["SOURCE_CLASS"] == g for g in ["Ribosomal", "Nonribosomal", "Synthetic"]})
masks.update({("Function (all)", g): df[g] for g in TARGET_GROUPS.values()})
order = [o for o in order if masks[o].sum() >= MIN_N]
fig, ax = plt.subplots(figsize=(10, 4.8))
data = [df.loc[masks[o], "LENGTH"] for o in order]
bp = ax.boxplot(data, showfliers=False, patch_artist=True, medianprops={"color": "black"})
cmap = {"Source organism (Natural)": "#2a9d8f", "Synthesis type": "#e9c46a", "Function (all)": "#6c8ebf"}
for patch, o in zip(bp["boxes"], order):
    patch.set_facecolor(cmap[o[0]]); patch.set_alpha(0.85)
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels([f"{o[1]}\n(n={len(d):,})" for o, d in zip(order, data)], rotation=60, ha="right", fontsize=8.5)
ax.set_ylabel("peptide length (aa)"); ax.set_ylim(0, 75)
ax.set_title("Length by class (box = quartiles, line = median; teal = source organism, yellow = synthesis type, blue = function)", fontsize=10)
fig.savefig(FIGS / "length_boxplots.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# length categories
cats = ["%Ultrashort (3-9)", "%Short (10-24)", "%Medium (25-50)", "%Long (51-100)"]
colors = ["#f4a261", "#e9c46a", "#2a9d8f", "#264653"]
s = summ[summ["grouping"].isin(["Source organism (Natural)", "Synthesis type", "Function (all)"])].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(10, 4.8)); bottom = np.zeros(len(s))
for c, col in zip(cats, colors):
    ax.bar(range(len(s)), s[c], bottom=bottom, color=col, label=c[1:]); bottom += s[c].values
ax.set_xticks(range(len(s))); ax.set_xticklabels(s["group"], rotation=60, ha="right", fontsize=8.5)
ax.set_ylabel("% of peptides in the class"); ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.12), frameon=False)
fig.savefig(FIGS / "length_category_share.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# function overlap
jac = pd.read_csv(TABLES / "function_overlap_jaccard.csv", index_col=0)
fig, ax = plt.subplots(figsize=(6.5, 5.6))
im = ax.imshow(jac.values, cmap="Blues", vmin=0, vmax=1)
ax.set_xticks(range(len(jac))); ax.set_xticklabels(jac.columns, rotation=60, ha="right", fontsize=8.5)
ax.set_yticks(range(len(jac))); ax.set_yticklabels(jac.index, fontsize=8.5)
for i in range(len(jac)):
    for j in range(len(jac)):
        v = jac.values[i, j]; ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.5, color="white" if v > 0.55 else "black")
fig.colorbar(im, shrink=0.8, label="Jaccard index (shared peptides / all peptides in either class)")
ax.set_title("Overlap between function classes")
fig.savefig(FIGS / "function_overlap.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("figures written to", FIGS)
