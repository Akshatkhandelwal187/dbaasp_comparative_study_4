"""Deliverable step 5 - the amino acid enrichment figure (run 04_aa_enrichment.py first).

Reads   tables/aa_enrichment.csv
Writes  figures/amino_acid_enrichment.png   ONE image, four panels with identical axes:
            columns = background (SwissProt proteins | non-AMP peptides), rows = Natural | Synthetic

How to read it: each bar says how much more or less often the amino acid occurs in the antimicrobial peptides than in the
background, in percent. 0 = same as the background, +100 % = twice as often, -50 % = half as often.
Amino acids are in one fixed order, hydrophilic -> hydrophobic (Kyte-Doolittle order of the Composition Profiler),
coloured along that order. Gray = not significant. * = significant (Composition Profiler, Bonferroni alpha = 0.05 / 20).
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
TITLE = {"SwissProt": f"Background: normal proteins (SwissProt, {n_bg['SwissProt']:,})",
         "Non-AMP": f"Background: non-antimicrobial peptides ({n_bg['Non-AMP']:,})"}

top = np.ceil(t["pct_change"].max() / 100) * 100          # same y range in all four panels
fig, axes = plt.subplots(2, 2, figsize=(15, 8.4), sharex=True, sharey=True)
fig.subplots_adjust(left=0.07, right=0.99, top=0.80, bottom=0.17, hspace=0.12, wspace=0.05)
for r, origin in enumerate(("Natural", "Synthetic")):
    for c, bg in enumerate(("SwissProt", "Non-AMP")):
        ax = axes[r, c]
        d = t[(t["origin"] == origin) & (t["background"] == bg)].set_index("amino_acid").loc[list(AA_ORDER)]
        for i, (aa, row) in enumerate(d.iterrows()):
            sig = row["result"] != "Not significant"
            ax.bar(i, row["pct_change"], width=0.72, color=HYDRO(i / 19) if sig else NS_GRAY, zorder=3)
            if sig:
                v = row["pct_change"]
                ax.annotate("*", (i, v), xytext=(0, -2 if v >= 0 else -9), textcoords="offset points", ha="center",
                            va="bottom", fontsize=12, color=INK2)
        ax.axhline(0, color=INK2, linewidth=1, zorder=4)
        ax.set_ylim(-100, top)
        ax.set_yticks(np.arange(-100, top + 1, 100))
        ax.set_yticklabels([f"{v:+.0f}%" if v else "0" for v in np.arange(-100, top + 1, 100)])
        ax.set_xlim(-0.6, 19.6)
        ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
        ax.text(0.01, 0.97, f"{origin}  (n = {n_pep[origin]:,})", transform=ax.transAxes, ha="left", va="top",
                fontsize=11.5, fontweight="bold")
        if r == 0:
            ax.set_title(TITLE[bg], loc="left", fontsize=11, color=INK2, pad=8)
        if c == 0:
            ax.set_ylabel("Difference from background", color=INK2)
for ax in axes[1]:
    ax.set_xticks(range(20))
    ax.set_xticklabels(list(AA_ORDER), fontsize=11)
fig.text(0.07, 0.955, "Amino acids in Natural and Synthetic antimicrobial peptides, compared with two backgrounds",
         fontsize=13.5, fontweight="bold", ha="left", va="top")
fig.text(0.07, 0.915, "Each bar: how much more or less often the amino acid occurs in the peptides than in the background.\n"
         "0 = same as the background  ·  +100% = twice as often  ·  −50% = half as often.",
         fontsize=10.5, color=INK2, ha="left", va="top", linespacing=1.6)
key = fig.add_axes([0.07, 0.085, 0.25, 0.018])
key.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=HYDRO)
key.set_xticks([]); key.set_yticks([])
for s in key.spines.values():
    s.set_visible(False)
fig.text(0.07, 0.068, "hydrophilic", fontsize=9, color=INK2, ha="left", va="top")
fig.text(0.32, 0.068, "hydrophobic", fontsize=9, color=INK2, ha="right", va="top")
fig.text(0.37, 0.094, "Amino acids are in the same order in every panel, coloured from hydrophilic to hydrophobic.\n"
         "Gray bar = not significant.   * = significant (Composition Profiler, 10,000 permutations, Bonferroni-corrected).",
         fontsize=9.5, color=MUTED, ha="left", va="center", linespacing=1.6)
fig.savefig(OUT / "amino_acid_enrichment.png")
print(f"figure written to {OUT}")
