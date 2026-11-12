"""Deliverable step 5 - amino acid enrichment figures (run 04_aa_enrichment.py first).

Reads   tables/aa_enrichment_natural_vs_synthetic.csv   (written by step 4: fold change, p-value and result per test)
Writes  figures/amino_acid_enrichment.png         the 20 amino acids: Natural (top) and Synthetic (bottom) in ONE image
        figures/amino_acid_group_enrichment.png   the 20 groups of the Composition Profiler: Natural (left), Synthetic (right)

Comparability
  * both panels of a figure have the same axes, the same scale and the same order of the items;
  * the amino acids follow the Composition Profiler's own hydrophobicity order (Kyte-Doolittle, `-X hydrophobicity_kyte`),
    from hydrophilic (R) to hydrophobic (I), and are colour-coded along that order;
  * significance is the one the Composition Profiler gives (permutation test, Bonferroni alpha = 0.05 / 40 = 0.00125):
    "*" = significant, "**" = p < 0.0001 (no random relabeling reached the observed effect; the smallest value
    10,000 iterations can give); not significant = gray bar and "ns".
Usage:  python 05_enrichment_figures.py [table.csv] [output_dir]     (defaults: the deliverable folders)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

DELIV = Path(__file__).resolve().parent.parent
TABLE = Path(sys.argv[1]) if len(sys.argv) > 1 else DELIV / "tables" / "aa_enrichment_natural_vs_synthetic.csv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else DELIV / "figures"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 200

# cprofiler.aminoacid.AminoAcid.get_order("hydrophobicity_kyte"): increasing hydrophobicity (hydrophilic -> hydrophobic)
AA_ORDER = "RKDENQHPYWSTGAMCFLVI"
BONFERRONI_ALPHA = 0.05 / 40
FIXED_YLIM = (-2.7, 3.3)         # log2 fold change; the same in every amino acid panel (0.16x ... 9.8x)

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
NS_GRAY = "#bdbcb4"
UP, DOWN = "#1b7f79", "#b8421e"                            # groups: enriched / depleted
HYDRO = LinearSegmentedColormap.from_list("hydro", ["#b8421e", "#e0a030", "#6fae7b", "#2a78d6", "#1d3f91"])
plt.rcParams.update({"font.family": "sans-serif", "figure.dpi": DPI, "savefig.dpi": DPI,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "font.size": 10})

table = pd.read_csv(TABLE)
ORIGINS = ["Natural", "Synthetic"]
N = {o: sum(1 for l in open(DELIV / "fasta" / f"amp_{o.lower()}_cdhit90.fasta") if l[0] == ">") for o in ORIGINS}


def significance(row, origin):
    """Mark from the Composition Profiler's own call: '' is never returned, 'ns' = not significant."""
    if row[f"result_{origin}"] == "Not significant":
        return "ns"
    return "**" if row[f"pvalue_{origin}"] < 1e-4 else "*"


def fold_ticks(ax, horizontal, lo, hi):
    ticks = list(range(lo, hi + 1))
    labels = [f"{2.0 ** t:g}×" for t in ticks]
    (ax.set_xticks if horizontal else ax.set_yticks)(ticks)
    (ax.set_xticklabels if horizontal else ax.set_yticklabels)(labels)


def style(ax, horizontal):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x" if horizontal else "y", color=GRID, linewidth=0.75)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=9)


def annotate(ax, pos, value, fold, mark, horizontal):
    """Fold change and significance mark at the tip of a bar."""
    color = MUTED if mark == "ns" else INK2
    text = f"{fold:.2f}\n{mark}" if not horizontal else f"{fold:.2f} {mark}"
    if horizontal:
        ax.annotate(text, (value, pos), xytext=(4 if value >= 0 else -4, 0), textcoords="offset points",
                    ha="left" if value >= 0 else "right", va="center", fontsize=7.5, color=color, annotation_clip=False)
    else:
        text = f"{mark}\n{fold:.2f}" if value >= 0 else f"{fold:.2f}\n{mark}"
        ax.annotate(text, (pos, value), xytext=(0, 3 if value >= 0 else -3), textcoords="offset points",
                    ha="center", va="bottom" if value >= 0 else "top", fontsize=7.5, color=color,
                    linespacing=1.0, annotation_clip=False)


# ------------------------------------------------------------------ 20 amino acids: one image, two panels
aa = table[table["type"] == "amino acid"].set_index("test_name").loc[list(AA_ORDER)].reset_index()
fig, axes = plt.subplots(2, 1, figsize=(11, 8.6), sharex=True, sharey=True)
fig.subplots_adjust(left=0.085, right=0.985, top=0.835, bottom=0.22, hspace=0.12)
for ax, origin in zip(axes, ORIGINS):
    for i, (_, row) in enumerate(aa.iterrows()):
        mark = significance(row, origin)
        value = np.log2(row[f"fold_{origin}"])
        color = NS_GRAY if mark == "ns" else HYDRO(i / (len(aa) - 1))
        ax.bar(i, value, width=0.7, color=color, zorder=3, linewidth=0)
        annotate(ax, i, value, row[f"fold_{origin}"], mark, horizontal=False)
    ax.axhline(0, color=AXIS, linewidth=0.75, zorder=2)
    ax.set_ylim(*FIXED_YLIM)
    ax.set_xlim(-0.6, len(aa) - 0.4)
    fold_ticks(ax, horizontal=False, lo=-2, hi=3)
    style(ax, horizontal=False)
    ax.set_ylabel("Frequency relative to SwissProt", labelpad=8)
    ax.text(0.005, 0.97, f"{origin}  (n = {N[origin]:,} peptides)", transform=ax.transAxes, ha="left", va="top",
            fontsize=11, fontweight="bold", color=INK)
axes[1].set_xticks(range(len(aa)))
axes[1].set_xticklabels(aa["test_name"], color=INK, fontsize=11)
fig.text(0.07, 0.955, "Amino acid enrichment in Natural and Synthetic peptides", fontsize=13, fontweight="bold", ha="left", va="top")
fig.text(0.07, 0.915, "Residue frequency in the peptides relative to SwissProt, log scale; same axes and same order in both panels.\n"
         "Amino acids ordered from hydrophilic to hydrophobic (Kyte-Doolittle order of the Composition Profiler) and coloured along it.",
         fontsize=9.5, color=INK2, ha="left", va="top", linespacing=1.5)
# colour key for the order
key = fig.add_axes([0.085, 0.16, 0.30, 0.016])
key.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=HYDRO)
key.set_yticks([]); key.set_xticks([])
for side in key.spines.values():
    side.set_visible(False)
fig.text(0.085, 0.145, "hydrophilic", fontsize=8.5, color=INK2, ha="left", va="top")
fig.text(0.385, 0.145, "hydrophobic", fontsize=8.5, color=INK2, ha="right", va="top")
fig.legend([Patch(facecolor=NS_GRAY)], ["gray: not significant (ns)"], frameon=False, fontsize=9, labelcolor=INK2,
           loc="lower left", bbox_to_anchor=(0.43, 0.145), handlelength=1.0, handleheight=1.0)
fig.text(0.07, 0.025, "Above/below each bar: fold change relative to SwissProt and significance. * = significant, ** = p < 0.0001, ns = not significant.\n"
         "Composition Profiler permutation test (10,000 relabelings of whole sequences), Bonferroni-corrected over 20 amino acids + 20 groups "
         "(alpha = 0.05 / 40 = 0.00125).\nBars above 1× = enriched, below 1× = depleted.",
         fontsize=8.5, color=MUTED, ha="left", va="bottom", linespacing=1.5)
fig.savefig(OUT / "amino_acid_enrichment.png")
plt.close(fig)

# ------------------------------------------------------------------ 20 groups: one image, two panels side by side
FIX = {"propencity": "propensity", "potenential": "potential"}              # typos in the tool's group names
grp = table[table["type"] == "group"].reset_index(drop=True)                  # the tool's own group order, same in both panels
labels = []
for t, r in zip(grp["test_name"], grp["residues"]):
    for wrong, right in FIX.items():
        t = t.replace(wrong, right)
    labels.append(f"{t}  ({r})")
fig, axes = plt.subplots(1, 2, figsize=(13, 8.2), sharex=True, sharey=True)
fig.subplots_adjust(left=0.30, right=0.985, top=0.82, bottom=0.17, wspace=0.08)
for ax, origin in zip(axes, ORIGINS):
    for i, (_, row) in enumerate(grp.iterrows()):
        mark = significance(row, origin)
        value = np.log2(row[f"fold_{origin}"])
        color = NS_GRAY if mark == "ns" else (UP if value > 0 else DOWN)
        ax.barh(i, value, height=0.7, color=color, zorder=3, linewidth=0)
        annotate(ax, i, value, row[f"fold_{origin}"], mark, horizontal=True)
    ax.axvline(0, color=AXIS, linewidth=0.75, zorder=2)
    ax.set_ylim(len(grp) - 0.4, -0.6)
    ax.set_xlim(-2.6, 2.6)
    fold_ticks(ax, horizontal=True, lo=-2, hi=2)
    style(ax, horizontal=True)
    ax.set_xlabel("Frequency relative to SwissProt", labelpad=8)
    ax.set_title(f"{origin}  (n = {N[origin]:,} peptides)", loc="left", fontsize=11, fontweight="bold", color=INK, pad=8)
axes[0].set_yticks(range(len(grp)))
axes[0].set_yticklabels(labels, color=INK2, fontsize=9)
fig.text(0.02, 0.955, "Amino acid group enrichment in Natural and Synthetic peptides", fontsize=13, fontweight="bold", ha="left", va="top")
fig.text(0.02, 0.915, "Groups of the Composition Profiler (residues in brackets), in the tool's order; same axes and same order in both panels.\n"
         "Fold change relative to SwissProt, log scale.", fontsize=9.5, color=INK2, ha="left", va="top", linespacing=1.5)
fig.legend([Patch(facecolor=UP), Patch(facecolor=DOWN), Patch(facecolor=NS_GRAY)],
           ["enriched", "depleted", "not significant (ns)"], frameon=False, fontsize=9.5, labelcolor=INK2, ncol=3,
           loc="lower left", bbox_to_anchor=(0.30, 0.04), handlelength=1.0, handleheight=1.0)
fig.text(0.02, 0.008, "Next to each bar: fold change and significance. * = significant, ** = p < 0.0001, ns = not significant "
         "(Composition Profiler permutation test, 10,000 relabelings, Bonferroni alpha = 0.05 / 40 = 0.00125).",
         fontsize=8.5, color=MUTED, ha="left", va="bottom")
fig.savefig(OUT / "amino_acid_group_enrichment.png")
plt.close(fig)
print(f"figures written to {OUT}")
