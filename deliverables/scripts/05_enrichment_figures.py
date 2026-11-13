"""Deliverable step 5 - the amino acid enrichment figure (run 04_aa_enrichment.py first).

Reads   tables/aa_enrichment.csv
Writes  figures/amino_acid_enrichment.png   ONE image, four panels with identical axes:
            columns = background (SwissProt proteins | non-AMP peptides), rows = Natural | Synthetic

Bar height = fractional difference (peptides - background) / background: 0 = same as the background, 1 = twice as often,
-0.5 = half as often. The number at each bar is the Composition Profiler's p-value. Amino acids are in one fixed order,
hydrophilic -> hydrophobic (Kyte-Doolittle order of the Composition Profiler), coloured along that order; gray = not significant.
Usage:  python 05_enrichment_figures.py [table.csv] [output_dir]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

DELIV = Path(__file__).resolve().parent.parent
TABLE = Path(sys.argv[1]) if len(sys.argv) > 1 else DELIV / "tables" / "aa_enrichment.csv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else DELIV / "figures"
OUT.mkdir(parents=True, exist_ok=True)

AA_ORDER = "RKDENQHPYWSTGAMCFLVI"      # cprofiler AminoAcid.get_order("hydrophobicity_kyte"): hydrophilic -> hydrophobic
NS_GRAY, GRID, INK, INK2, MUTED = "#c4c3bb", "#e1e0d9", "#0b0b0b", "#52514e", "#898781"
HYDRO = LinearSegmentedColormap.from_list("hydro", ["#b8421e", "#e0a030", "#6fae7b", "#2a78d6", "#1d3f91"])
plt.rcParams.update({"font.family": "sans-serif", "figure.dpi": 200, "savefig.dpi": 200, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white", "text.color": INK, "font.size": 10})

t = pd.read_csv(TABLE)
n_pep = {o: sum(1 for l in open(DELIV / "fasta" / f"amp_{o.lower()}_cdhit90.fasta") if l[0] == ">") for o in ("Natural", "Synthetic")}
n_bg = {"SwissProt": sum(1 for l in open(DELIV / "background" / "swissprot51_5k.fasta") if l[0] == ">"),
        "Non-AMP": sum(1 for l in open(DELIV / "background" / "pepbench_nonamp_background.fasta") if l[0] == ">")}
TITLE = {"SwissProt": "Background: SwissProt proteins", "Non-AMP": "Background: non-AMP peptides"}

def p_text(p):
    return "<0.0001" if p == 0 else f"{p:.3g}"


top = np.ceil(t["effect"].max()) + 1.0                                  # same y range in all four panels
bottom = min(np.floor(t["effect"].min()), -1) - 0.8
fig, axes = plt.subplots(2, 2, figsize=(15, 8.6), sharex=True, sharey=True)
fig.subplots_adjust(left=0.07, right=0.99, top=0.90, bottom=0.14, hspace=0.12, wspace=0.05)
for r, origin in enumerate(("Natural", "Synthetic")):
    for c, bg in enumerate(("SwissProt", "Non-AMP")):
        ax = axes[r, c]
        d = t[(t["origin"] == origin) & (t["background"] == bg)].set_index("amino_acid").loc[list(AA_ORDER)]
        for i, (aa, row) in enumerate(d.iterrows()):
            sig = row["result"] != "Not significant"
            v = row["effect"]
            ax.bar(i, v, width=0.72, color=HYDRO(i / 19) if sig else NS_GRAY, zorder=3)
            ax.annotate(p_text(row["pvalue"]), (i, v), xytext=(0, 3 if v >= 0 else -3), textcoords="offset points",
                        ha="center", va="bottom" if v >= 0 else "top", rotation=90, fontsize=7, color=INK2)
        ax.axhline(0, color=INK2, linewidth=1, zorder=4)
        ax.set_ylim(bottom, top)
        ticks = np.arange(-1, top, 1)
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{v:g}" for v in ticks])
        ax.set_xlim(-0.6, 19.6)
        ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
        ax.text(0.01, 0.97, f"{origin}  (n = {n_pep[origin]:,})", transform=ax.transAxes, ha="left", va="top",
                fontsize=11.5, fontweight="bold")
        if r == 0:
            ax.set_title(TITLE[bg], loc="left", fontsize=11, color=INK2, pad=8)
        if c == 0:
            ax.set_ylabel("Fractional difference", color=INK2)
for ax in axes[1]:
    ax.set_xticks(range(20))
    ax.set_xticklabels(list(AA_ORDER), fontsize=11)
fig.text(0.07, 0.965, "Amino acid enrichment", fontsize=14, fontweight="bold", ha="left", va="top")
key = fig.add_axes([0.07, 0.065, 0.20, 0.016])
key.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=HYDRO)
key.set_xticks([]); key.set_yticks([])
for sp in key.spines.values():
    sp.set_visible(False)
fig.text(0.07, 0.05, "hydrophilic", fontsize=9, color=INK2, ha="left", va="top")
fig.text(0.27, 0.05, "hydrophobic", fontsize=9, color=INK2, ha="right", va="top")
fig.text(0.31, 0.073, "gray = not significant", fontsize=9, color=MUTED, ha="left", va="center")
fig.savefig(OUT / "amino_acid_enrichment.png")
print(f"figure written to {OUT}")
