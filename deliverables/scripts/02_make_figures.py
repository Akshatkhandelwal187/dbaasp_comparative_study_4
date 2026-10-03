"""Deliverable step 2 - length figures for the non-redundant AMP dataset (run 01_remove_redundancy.py first).

Reads  deliverables/tables/peptide_metadata.csv
Writes deliverables/figures/length_histogram_overlap.png         3 overlapping histograms, each as % of its own group
       deliverables/figures/length_histogram_overlap_counts.png  same, as peptide counts
       deliverables/figures/length_category_fraction_bar.png     share of Natural / Synthetic peptides per length category
       deliverables/tables/length_summary.csv                    numbers behind the histograms
       deliverables/tables/length_category_fractions.csv         numbers behind the bar chart

Colours: Natural = blue, Synthetic = orange (categorical slots 1 and 2, colour-blind validated); the entire dataset
is the grey outline.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path as MplPath

DELIV = Path(__file__).resolve().parent.parent
HIST_MAX = 100                      # histogram x-axis ends here; longer peptides are counted in the footnote
DPI = 200

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
NATURAL, SYNTHETIC, ENTIRE = "#2a78d6", "#eb6834", "#898781"

# length categories, same boundaries as the Length_Class column of the DBAASP files (ascending length)
CATEGORIES = [("1 aa", "(out of range)"), ("Ultrashort", "2–9 aa"), ("Short", "10–24 aa"), ("Medium", "25–50 aa"),
              ("Long", "51–100 aa"), (">100 aa", "")]
CLASS_NAMES = ["<2 aa (out of range)", "Ultrashort (2-9 aa)", "Short (10-24 aa)", "Medium (25-50 aa)",
               "Long (51-100 aa)", ">100 aa"]

plt.rcParams.update({"font.family": "sans-serif", "figure.dpi": DPI, "savefig.dpi": DPI,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "font.size": 10})

meta = pd.read_csv(DELIV / "tables" / "peptide_metadata.csv")
groups = {"Entire dataset": meta["LENGTH"], "Natural": meta.loc[meta["ORIGIN"] == "Natural", "LENGTH"],
          "Synthetic": meta.loc[meta["ORIGIN"] == "Synthetic", "LENGTH"]}


def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color=GRID, linewidth=0.75)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=9)


def heading(fig, title, subtitle):
    fig.text(0.07, 0.955, title, fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(0.07, 0.905, subtitle, fontsize=9.5, color=INK2, ha="left", va="top")


# ------------------------------------------------------------------ tables
rows = []
for name, L in groups.items():
    rows.append({"group": name, "n": len(L), "mean": round(L.mean(), 1), "median": L.median(), "Q1": L.quantile(.25),
                 "Q3": L.quantile(.75), "min": L.min(), "max": L.max(), "n_longer_than_100": int((L > HIST_MAX).sum())})
pd.DataFrame(rows).to_csv(DELIV / "tables" / "length_summary.csv", index=False)

frac = pd.DataFrame({"length_category": CLASS_NAMES})
for name, L in groups.items():
    key = name.split()[0]
    counts = L.map(lambda n: 0 if n < 2 else 1 if n < 10 else 2 if n < 25 else 3 if n < 51 else 4 if n < 101 else 5)
    frac[f"{key}_n"] = counts.value_counts().reindex(range(6), fill_value=0).values
    frac[f"{key}_fraction"] = (frac[f"{key}_n"] / len(L)).round(4)
frac.to_csv(DELIV / "tables" / "length_category_fractions.csv", index=False)
# the categories computed here have to agree with the LENGTH_CLASS column written in step 1
check = meta.groupby("LENGTH_CLASS").size()
assert all(check[c] == n for c, n in zip(CLASS_NAMES, frac["Entire_n"]) if n > 0)


# ------------------------------------------------------------------ overlapping histogram
def histogram(filename, as_percent):
    fig, ax = plt.subplots(figsize=(9, 5.2))
    fig.subplots_adjust(left=0.085, right=0.97, top=0.82, bottom=0.15)
    bins = np.arange(0.5, HIST_MAX + 1.5, 1)                    # one bin per amino-acid length
    handles = []
    for z, (name, color) in enumerate([("Natural", NATURAL), ("Synthetic", SYNTHETIC), ("Entire dataset", INK2)], start=1):
        L = groups[name]
        w = np.full(len(L), 100 / len(L)) if as_percent else None
        if name != "Entire dataset":                                       # filled, translucent, with a 2 px outline
            ax.hist(L, bins=bins, weights=w, histtype="stepfilled", color=color, alpha=0.40, linewidth=0, zorder=z)
            ax.hist(L, bins=bins, weights=w, histtype="step", color=color, linewidth=1.5, zorder=z)
            swatch = Patch(facecolor=color, alpha=0.55, edgecolor=color, linewidth=1.5)
        else:                                                              # outline only, drawn on top
            ax.hist(L, bins=bins, weights=w, histtype="step", color=color, linewidth=1.5, zorder=z + 1)
            swatch = Line2D([], [], color=color, linewidth=1.5)
        handles.append((swatch, f"{name}  (n = {len(L):,}; median {L.median():g} aa)"))
    ax.set_xlim(0, HIST_MAX)
    ax.set_xticks(range(0, HIST_MAX + 1, 10))
    ax.set_xlabel("Peptide length (amino acids)", labelpad=8)
    ax.set_ylabel("% of the group's peptides per 1-aa bin" if as_percent else "Number of peptides per 1-aa bin",
                  labelpad=8)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    style_axes(ax)
    ax.legend(*zip(*handles[::-1]), loc="upper right", frameon=False, fontsize=9.5, labelcolor=INK2, handlelength=1.6)
    n_long = {k: int((v > HIST_MAX).sum()) for k, v in groups.items()}
    fig.text(0.07, 0.035, f"Not shown: {n_long['Entire dataset']} peptides longer than {HIST_MAX} aa "
             f"(Natural {n_long['Natural']}, Synthetic {n_long['Synthetic']}; longest {groups['Entire dataset'].max()} aa).",
             fontsize=8.5, color=MUTED, ha="left", va="bottom")
    heading(fig, "Peptide length distribution: entire dataset vs Natural vs Synthetic",
            "Non-redundant DBAASP peptides (CD-HIT, 90 % identity). " +
            ("Each histogram is scaled to 100 % of its own group." if as_percent else "Counts of peptides; the entire dataset is Natural + Synthetic."))
    fig.savefig(DELIV / "figures" / filename)
    plt.close(fig)


histogram("length_histogram_overlap.png", as_percent=True)
histogram("length_histogram_overlap_counts.png", as_percent=False)


# ------------------------------------------------------------------ grouped bar chart
def rounded_top_bar(ax, x0, x1, height, color, radius_pt=3):
    """Bar with a ~4 px rounded data end and a square base (radius is measured in points, not data units)."""
    if height <= 0:
        return
    inv = ax.transData.inverted()
    px = ax.figure.dpi / 72 * radius_pt
    origin = inv.transform((0, 0))
    rx = abs(inv.transform((px, 0))[0] - origin[0])
    ry = min(abs(inv.transform((0, px))[1] - origin[1]), height)
    t = np.linspace(0, np.pi / 2, 10)
    left = np.c_[x0 + rx - rx * np.cos(t), height - ry + ry * np.sin(t)]      # 180 -> 90 degrees
    right = np.c_[x1 - rx + rx * np.sin(t), height - ry + ry * np.cos(t)]     # 90 -> 0 degrees
    verts = np.vstack([[x0, 0], left, right, [x1, 0], [x0, 0]])      # last vertex = CLOSEPOLY marker
    ax.add_patch(PathPatch(MplPath(verts, closed=True), facecolor=color, linewidth=0, zorder=3))


def bar_chart(filename):
    fig, ax = plt.subplots(figsize=(9, 5.2))
    fig.subplots_adjust(left=0.085, right=0.97, top=0.82, bottom=0.17)
    n_cat = len(CATEGORIES)
    top = np.ceil(max(frac["Natural_fraction"].max(), frac["Synthetic_fraction"].max()) * 100 / 10) * 10
    ax.set_xlim(-0.6, n_cat - 0.4)
    ax.set_ylim(0, top)
    fig.canvas.draw()                                              # fixes the layout so points can be turned into data units
    inv = ax.transData.inverted()
    one_pt = abs(inv.transform((ax.figure.dpi / 72, 0))[0] - inv.transform((0, 0))[0])
    bar_w, gap = 18 * one_pt, 1.5 * one_pt                         # bars <= 24 px thick, 2 px apart
    for i in range(n_cat):
        for side, (key, color) in enumerate([("Natural", NATURAL), ("Synthetic", SYNTHETIC)]):
            x0 = i - gap / 2 - bar_w if side == 0 else i + gap / 2
            pct = frac[f"{key}_fraction"][i] * 100
            rounded_top_bar(ax, x0, x0 + bar_w, pct, color)
            label = "0" if pct == 0 else "<0.1" if pct < 0.05 else f"{pct:.1f}"
            ax.annotate(label, (x0 + bar_w / 2, pct), xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8, color=INK2)
    ax.set_xticks(range(n_cat))
    ax.set_xticklabels([f"{a}\n{b}" if b else a for a, b in CATEGORIES], color=INK2, fontsize=9.5, linespacing=1.5)
    ax.set_xlabel("Length category", labelpad=8)
    ax.set_ylabel("% of the group's peptides", labelpad=8)
    style_axes(ax)
    ax.tick_params(axis="x", pad=6)
    n_nat, n_syn = len(groups["Natural"]), len(groups["Synthetic"])
    ax.legend([Line2D([], [], marker="s", linestyle="", markersize=8, color=c) for c in (NATURAL, SYNTHETIC)],
              [f"Natural  (n = {n_nat:,})", f"Synthetic  (n = {n_syn:,})"], loc="upper right", frameon=False,
              fontsize=9.5, labelcolor=INK2, handletextpad=0.3)
    heading(fig, "Natural vs Synthetic peptides across length categories",
            "Share of each group's non-redundant peptides (CD-HIT, 90 % identity) in each category; each group sums to 100 %.")
    fig.savefig(DELIV / "figures" / filename)
    plt.close(fig)


bar_chart("length_category_fraction_bar.png")
print(pd.DataFrame(rows).to_string(index=False))
print(frac.to_string(index=False))
