"""Class comparisons step 4 - enrichment figures and the overview figures (run 02 and 03 first).

For each comparison:
  figures/amino_acid_enrichment.png   same design as deliverables/figures/amino_acid_enrichment.png: one row per class,
                                      columns = background (SwissProt | non-AMP), identical axes in every panel, amino acids
                                      hydrophilic -> hydrophobic, bar = fractional difference, number = p-value, gray = not significant
  figures/enrichment_heatmap.png      compact view of the same table: classes x amino acids, colour = % change, * = significant
Overview (overview/figures/):
  length_boxplots.png                 length of every tested class of every comparison
  function_overlap.png, target_object_overlap.png   Jaccard overlap of the classes that can share peptides
Usage: python 04_enrichment_figures.py [comparison ...]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from common import AA_ORDER, CC, COMPARISONS, GRID, INK, INK2, MUTED, OVERVIEW, SERIES, load_table, tested

NS_GRAY = "#c4c3bb"
CAP = 5.0                                    # bars above this fractional difference are cut (value written at the bar)
HYDRO = LinearSegmentedColormap.from_list("hydro", ["#b8421e", "#e0a030", "#6fae7b", "#2a78d6", "#1d3f91"])
TITLE = {"SwissProt": "Background: SwissProt proteins", "Non-AMP": "Background: non-AMP peptides"}
plt.rcParams.update({"font.family": "sans-serif", "figure.dpi": 200, "savefig.dpi": 200, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white", "text.color": INK, "font.size": 10})


def p_text(p):
    return "<0.0001" if p == 0 else f"{p:.3g}"


def bars(key, t, out):
    names = list(dict.fromkeys(t["class"]))
    n_of = t.drop_duplicates("class").set_index("class")["n_peptides"]
    cap = CAP if t["effect"].max() > CAP else None                     # very large bars are cut at CAP and labelled
    top = (CAP + 1.9) if cap else np.ceil(t["effect"].max()) + 1.0
    bottom = min(np.floor(t["effect"].min()), -1) - 1.45
    row_h = 2.45
    H = row_h * len(names) + 1.6
    fig, axes = plt.subplots(len(names), 2, figsize=(15, H), sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(left=0.07, right=0.99, top=1 - 0.95 / H, bottom=0.85 / H, hspace=0.14, wspace=0.05)
    for r, name in enumerate(names):
        for c, bg in enumerate(("SwissProt", "Non-AMP")):
            ax = axes[r, c]
            d = t[(t["class"] == name) & (t["background"] == bg)].set_index("amino_acid").loc[list(AA_ORDER)]
            for i, (aa, row) in enumerate(d.iterrows()):
                v, sig = row["effect"], row["result"] != "Not significant"
                h = min(v, cap) if cap else v
                ax.bar(i, h, width=0.72, color=HYDRO(i / 19) if sig else NS_GRAY, zorder=3)
                label = p_text(row["pvalue"])
                if cap and v > cap:                                     # cut bar: zig-zag mark + the real value
                    ax.plot([i - 0.36, i - 0.12, i + 0.12, i + 0.36], [h - 0.35, h - 0.15, h - 0.45, h - 0.25],
                            color="white", linewidth=1.6, zorder=4)
                    ax.text(i, h - 0.75, f"{v:+.1f}", ha="center", va="top", rotation=90, fontsize=7, color="white",
                            fontweight="bold", zorder=5)
                ax.annotate(label, (i, h), xytext=(0, 2 if v >= 0 else -2), textcoords="offset points",
                            ha="center", va="bottom" if v >= 0 else "top", rotation=90, fontsize=6, color=INK2)
            ax.axhline(0, color=INK2, linewidth=1, zorder=4)
            ax.set_ylim(bottom, top)
            step = 1 if top - bottom <= 8 else 2
            ticks = np.arange(-1, top, step)
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{v:g}" for v in ticks])
            ax.set_xlim(-0.6, 19.6)
            ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.tick_params(length=0, labelsize=8.5)
            ax.text(0.01, 0.97, f"{name}  (n = {n_of[name]:,})", transform=ax.transAxes, ha="left", va="top",
                    fontsize=10.5, fontweight="bold")
            if r == 0:
                ax.set_title(TITLE[bg], loc="left", fontsize=11, color=INK2, pad=8)
            if c == 0:
                ax.set_ylabel("Fractional difference", color=INK2, fontsize=9)
    for ax in axes[-1]:
        ax.set_xticks(range(20))
        ax.set_xticklabels(list(AA_ORDER), fontsize=10.5)
    fig.text(0.07, 1 - 0.18 / H, f"Amino acid enrichment: {COMPARISONS[key]['title']}", fontsize=14, fontweight="bold",
             ha="left", va="top")
    key_ax = fig.add_axes([0.07, 0.30 / H, 0.20, 0.13 / H])
    key_ax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=HYDRO)
    key_ax.set_xticks([]); key_ax.set_yticks([])
    for sp in key_ax.spines.values():
        sp.set_visible(False)
    fig.text(0.07, 0.25 / H, "hydrophilic", fontsize=9, color=INK2, ha="left", va="top")
    fig.text(0.27, 0.25 / H, "hydrophobic", fontsize=9, color=INK2, ha="right", va="top")
    fig.text(0.31, 0.37 / H, "gray = not significant (Bonferroni, p >= 0.0025); number = p-value"
             + (f"; bars above +{CAP:g} are cut, their value is written inside the bar" if cap else ""), fontsize=9, color=MUTED,
             ha="left", va="center")
    fig.savefig(out)
    plt.close(fig)


def heatmap(key, t, out, lim=400):
    names = list(dict.fromkeys(t["class"]))
    n_of = t.drop_duplicates("class").set_index("class")["n_peptides"]
    fig, axes = plt.subplots(1, 2, figsize=(2.2 + 0.55 * len(names) * 2, 7.4), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.88, top=0.86, bottom=0.27, wspace=0.06)
    for ax, bg in zip(axes, ["SwissProt", "Non-AMP"]):
        sub = t[t["background"] == bg]
        mat = sub.pivot(index="amino_acid", columns="class", values="pct_change").loc[list(AA_ORDER), names].values
        sig = (sub.pivot(index="amino_acid", columns="class", values="result").loc[list(AA_ORDER), names] != "Not significant").values
        im = ax.imshow(np.clip(mat, -lim, lim), cmap="RdBu_r", norm=TwoSlopeNorm(0, -100, lim), aspect="auto")
        for i in range(20):
            for j in range(len(names)):
                if sig[i, j]:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=11, color=INK)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels([f"{n} ({n_of[n]:,})" for n in names], rotation=55, ha="right", fontsize=8.5)
        ax.set_yticks(range(20)); ax.set_yticklabels(list(AA_ORDER))
        ax.set_title(TITLE[bg], fontsize=10.5, color=INK2, loc="left")
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
    cax = fig.add_axes([0.90, 0.30, 0.015, 0.52])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(f"% change vs background (colour capped at +{lim} %)", fontsize=8.5)
    fig.text(0.06, 0.97, f"Amino acid enrichment overview: {COMPARISONS[key]['title']}", fontsize=12.5, fontweight="bold",
             ha="left", va="top")
    fig.text(0.06, 0.925, "* = significant (Composition Profiler permutation test, Bonferroni alpha = 0.05/20); "
             "exact values in tables/aa_enrichment.csv", fontsize=8.5, color=INK2, ha="left", va="top")
    fig.savefig(out)
    plt.close(fig)


keys = sys.argv[1:] or list(COMPARISONS)
for key in keys:
    t = pd.read_csv(CC / key / "tables" / "aa_enrichment.csv")
    bars(key, t, CC / key / "figures" / "amino_acid_enrichment.png")
    heatmap(key, t, CC / key / "figures" / "enrichment_heatmap.png")
    print("figures written for", key)

# ---------------------------------------------------------------- overview figures
df = load_table()
(OVERVIEW / "figures").mkdir(parents=True, exist_ok=True)
groups, labels, seps = [], [], []
for key in COMPARISONS:
    cl, _ = tested(df, key)
    if groups:
        seps.append(len(groups) + 0.5)
    for name, mask in cl:
        groups.append(df.loc[mask, "LENGTH"])
        labels.append(f"{name} ({int(mask.sum()):,})")
fig, ax = plt.subplots(figsize=(17, 5.6))
fig.subplots_adjust(left=0.05, right=0.99, top=0.86, bottom=0.37)
bp = ax.boxplot(groups, showfliers=False, patch_artist=True, widths=0.6, medianprops={"color": INK, "linewidth": 1.4},
                whiskerprops={"color": MUTED}, capprops={"color": MUTED}, boxprops={"linewidth": 0})
for patch in bp["boxes"]:
    patch.set_facecolor(SERIES); patch.set_alpha(0.45)
start = 0.5
for s, key in zip(seps + [len(groups) + 0.5], COMPARISONS):
    ax.text((start + s) / 2, 1.02, COMPARISONS[key]["title"].replace(": target group", "").replace(" (mechanism, all peptides)", ""),
            transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8.5, color=INK2)
    if s < len(groups):
        ax.axvline(s, color=GRID, linewidth=1.2)
    start = s
ax.set_xticks(range(1, len(groups) + 1))
ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
ax.set_ylabel("Peptide length (aa)", color=INK2)
ax.set_ylim(0, 75)
ax.grid(axis="y", color=GRID, linewidth=0.75); ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.text(0.05, 0.97, "Peptide length by class (box = quartiles, line = median, whiskers = 1.5 IQR; outliers not drawn)",
         fontsize=12.5, fontweight="bold", ha="left", va="top")
fig.savefig(OVERVIEW / "figures" / "length_boxplots.png")
plt.close(fig)

for name, csv, title in [("function_overlap", "function_overlap_jaccard.csv", "Overlap between target groups (all peptides)"),
                         ("target_object_overlap", "target_object_overlap_jaccard.csv", "Overlap between target objects (all peptides)")]:
    j = pd.read_csv(OVERVIEW / "tables" / csv, index_col=0)
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    fig.subplots_adjust(left=0.27, right=0.93, top=0.90, bottom=0.27)
    im = ax.imshow(j.values, cmap=LinearSegmentedColormap.from_list("b", ["#f4f8fd", "#86b6ef", "#2a78d6", "#0d366b"]), vmin=0, vmax=1)
    ax.set_xticks(range(len(j))); ax.set_xticklabels(j.columns, rotation=55, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(j))); ax.set_yticklabels(j.index, fontsize=8.5)
    for a in range(len(j)):
        for b in range(len(j)):
            v = j.values[a, b]
            ax.text(b, a, f"{v:.2f}", ha="center", va="center", fontsize=7.5, color="white" if v > 0.55 else INK)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    fig.colorbar(im, shrink=0.75, label="Jaccard index (shared / all peptides in either class)")
    fig.text(0.27, 0.97, title, fontsize=12, fontweight="bold", ha="left", va="top")
    fig.savefig(OVERVIEW / "figures" / f"{name}.png")
    plt.close(fig)
print("overview figures written")
