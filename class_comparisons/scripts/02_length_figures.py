"""Class comparisons step 2 - length figures and tables for every comparison (run 01 first).

For each comparison folder (source_organism, synthesis_type, function_all, function_natural, target_object):
  figures/length_histogram_overlap.png         one panel per class: the class as bars + smooth line (% of the class per 1-aa bin),
                                               the whole population as a gray outline + gray line, medians as dashed lines.
                                               Same bins, smoothing and colours as deliverables/figures/length_histogram_overlap.png
  figures/length_histogram_overlap_counts.png  the same as peptide counts
  figures/length_category_fraction_bar.png     share of every class in the four length categories (100 % stacked bars)
  tables/length_summary.csv                    n, mean, median, quartiles, min, max per class (+ the population)
  tables/length_category_fractions.csv         numbers behind the bar chart
  tables/length_tests.csv                      each class vs all OTHER peptides of its population: Mann-Whitney U,
                                               rank-biserial r (+ = class longer), p, Bonferroni p (over the classes of the comparison);
                                               Kruskal-Wallis across the classes where classes do not overlap (organism, synthesis type)
Only classes with n >= MIN_N are drawn and tested.
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import (AXIS, CATEGORIES, CC, COMPARISONS, GRID, INK, INK2, LENGTH_RAMP, MUTED, SERIES, SURFACE, load_table, tested)

HIST_MAX, SMOOTH_BW, DPI = 100, 3.0, 200
plt.rcParams.update({"font.family": "sans-serif", "figure.dpi": DPI, "savefig.dpi": DPI, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "font.size": 10})
NON_OVERLAPPING = {"source_organism", "synthesis_type"}
df = load_table()


def smooth_curve(lengths, x, bandwidth=SMOOTH_BW):
    """Gaussian kernel density as % per 1-aa bin, mirrored at 0.5 aa (same as deliverables/scripts/02_make_figures.py)."""
    values, n = np.unique(lengths, return_counts=True)
    kernel = lambda c: np.exp(-0.5 * ((x[:, None] - c[None, :]) / bandwidth) ** 2)
    return 100 * (kernel(values) + kernel(1 - values)) @ n / (n.sum() * bandwidth * np.sqrt(2 * np.pi))


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color=GRID, linewidth=0.75)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=8.5)


def category(lengths):
    return pd.cut(lengths, [0, 9, 24, 50, 100], labels=[c[0] for c in CATEGORIES])


def histograms(key, cl, pop_name, popL, out, as_percent):
    ncol = 2 if len(cl) <= 4 else 3
    nrow = math.ceil(len(cl) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 2.9 * nrow + 0.9), sharex=True, squeeze=False)
    fig.subplots_adjust(left=0.07 if ncol == 3 else 0.09, right=0.98, top=1 - 1.15 / (2.9 * nrow + 0.9),
                        bottom=0.55 / (2.9 * nrow + 0.9), hspace=0.38, wspace=0.16)
    bins = np.arange(0.5, HIST_MAX + 1.5, 1)
    x = np.linspace(0.5, HIST_MAX, 400)
    ymax = 0
    for ax, (name, mask) in zip(axes.flat, cl):
        L = df.loc[mask, "LENGTH"]
        scale = (lambda n: 100 / n) if as_percent else (lambda n: 1)
        ax.hist(popL, bins=bins, weights=np.full(len(popL), scale(len(popL)) * (1 if as_percent else len(L) / len(popL))),
                histtype="step", color=MUTED, linewidth=0.7, zorder=1)
        ax.plot(x, smooth_curve(popL, x) * (1 if as_percent else len(L) / 100), color=INK2, linewidth=1.3, zorder=4)
        ax.hist(L, bins=bins, weights=np.full(len(L), scale(len(L))), histtype="stepfilled", color=SERIES, alpha=0.35,
                linewidth=0, zorder=2)
        ax.plot(x, smooth_curve(L, x) * (1 if as_percent else len(L) / 100), color=SERIES, linewidth=2, zorder=5)
        cm, pm = L.median(), popL.median()
        marks = [(cm, SERIES, -3, "right"), (pm, INK2, 3, "left")] if cm <= pm else [(cm, SERIES, 3, "left"), (pm, INK2, -3, "right")]
        for med, col, dx, ha in marks:
            if col == INK2 and pm == cm:
                continue
            ax.axvline(med, color=col, linewidth=1.2, linestyle="--", zorder=6)
            ax.annotate(f"{med:g}", (med, 1.0), xycoords=ax.get_xaxis_transform(), xytext=(dx, -2), textcoords="offset points",
                        ha=ha, va="top", fontsize=8.5, color=col, fontweight="bold")
        ax.set_title(f"{name}  (n = {len(L):,}; median {L.median():g} aa)", loc="left", fontsize=10, color=INK, pad=14)
        style(ax)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.08)
    for ax in axes.flat[len(cl):]:
        ax.axis("off")
    for ax in axes.flat[:len(cl)]:
        ax.set_xlim(0, HIST_MAX)
        ax.set_xticks(range(0, HIST_MAX + 1, 20))
    for r in range(nrow):
        axes[r, 0].set_ylabel("% of class per 1-aa bin" if as_percent else "peptides per 1-aa bin", fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("Peptide length (amino acids)", fontsize=9)
    for c in range(ncol):                                  # x labels on the last filled row of each column
        filled = [r for r in range(nrow) if r * ncol + c < len(cl)]
        if filled and filled[-1] != nrow - 1:
            axes[filled[-1], c].set_xlabel("Peptide length (amino acids)", fontsize=9)
            axes[filled[-1], c].tick_params(labelbottom=True)
    h = 2.9 * nrow + 0.9
    fig.text(0.07 if ncol == 3 else 0.09, 1 - 0.18 / h, f"Peptide length distribution: {COMPARISONS[key]['title']}",
             fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(0.07 if ncol == 3 else 0.09, 1 - 0.50 / h,
             f"blue = class (bars + smooth line), gray = {pop_name} (n = {len(popL):,}; median {popL.median():g} aa)"
             + ("" if as_percent else ", scaled to the class size") + "; dashed lines = medians",
             fontsize=9, color=INK2, ha="left", va="top")
    fig.savefig(out)
    plt.close(fig)


def category_bars(key, frac, pop_name, out):
    labels = list(frac["class"])
    fig, ax = plt.subplots(figsize=(max(6.5, 1.0 + 0.95 * len(labels)), 5.4))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.84, bottom=0.25)
    bottom = np.zeros(len(labels))
    for (cname, rng, _, _), col in zip(CATEGORIES, LENGTH_RAMP):
        v = frac[f"%{cname}"].values
        ax.bar(range(len(labels)), v, bottom=bottom, width=0.62, color=col, edgecolor=SURFACE, linewidth=1.5, zorder=3,
               label=f"{cname}  {rng}")
        for i, (b, h) in enumerate(zip(bottom, v)):
            if h >= 6:
                ax.text(i, b + h / 2, f"{h:.0f}", ha="center", va="center", fontsize=8,
                        color="white" if col in LENGTH_RAMP[2:] else INK)
        bottom += v
    ax.axvline(0.5, color=AXIS, linewidth=1, linestyle=":", zorder=1)       # first bar = the whole population (reference)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([f"{c}\n(n = {n:,})" for c, n in zip(labels, frac["n_peptides"])], rotation=35, ha="right",
                       fontsize=8.5, color=INK2)
    ax.get_xticklabels()[0].set_fontweight("bold")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of the class's peptides", labelpad=6)
    style(ax)
    ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, fontsize=9, labelcolor=INK2, handlelength=1.2)
    fig.text(0.09, 0.975, f"Peptide length categories: {COMPARISONS[key]['title']}", fontsize=13, fontweight="bold",
             ha="left", va="top")
    fig.savefig(out)
    plt.close(fig)


overview = []
for key, spec in COMPARISONS.items():
    (CC / key / "figures").mkdir(parents=True, exist_ok=True)
    (CC / key / "tables").mkdir(parents=True, exist_ok=True)
    cl, pop = tested(df, key)
    pop_name = "all Natural peptides" if spec["population"] == "Natural" else "all peptides"
    popL = df.loc[pop, "LENGTH"]

    summ, fr, tests = [], [], []
    for name, mask in [(pop_name, pop)] + cl:
        L = df.loc[mask, "LENGTH"]
        summ.append({"class": name, "n_peptides": len(L), "mean": round(L.mean(), 1), "median": L.median(), "Q1": L.quantile(.25),
                     "Q3": L.quantile(.75), "min": L.min(), "max": L.max()})
        shares = category(L).value_counts(normalize=True).reindex([c[0] for c in CATEGORIES], fill_value=0) * 100
        fr.append({"class": name, "n_peptides": len(L), **{f"%{c}": round(v, 1) for c, v in shares.items()}})
        if name != pop_name:
            rest = df.loc[pop & ~mask, "LENGTH"]
            u = stats.mannwhitneyu(L, rest, alternative="two-sided")
            tests.append({"test": "Mann-Whitney U, class vs rest of population", "class": name, "n_class": len(L), "n_rest": len(rest),
                          "median_class": L.median(), "median_rest": rest.median(),
                          "rank_biserial_r": round(2 * u.statistic / (len(L) * len(rest)) - 1, 3), "p_value": u.pvalue})
    tests = pd.DataFrame(tests)
    tests["p_bonferroni"] = np.minimum(tests["p_value"] * len(tests), 1)
    if key in NON_OVERLAPPING:
        h = stats.kruskal(*[df.loc[m, "LENGTH"] for _, m in cl])
        tests = pd.concat([tests, pd.DataFrame([{"test": "Kruskal-Wallis across the classes", "class": "all tested classes",
                                                 "n_class": sum(int(m.sum()) for _, m in cl), "p_value": h.pvalue,
                                                 "H": round(h.statistic, 1)}])], ignore_index=True)
    pd.DataFrame(summ).to_csv(CC / key / "tables" / "length_summary.csv", index=False)
    frac = pd.DataFrame(fr)
    frac.to_csv(CC / key / "tables" / "length_category_fractions.csv", index=False)
    tests.to_csv(CC / key / "tables" / "length_tests.csv", index=False)

    histograms(key, cl, pop_name, popL, CC / key / "figures" / "length_histogram_overlap.png", as_percent=True)
    histograms(key, cl, pop_name, popL, CC / key / "figures" / "length_histogram_overlap_counts.png", as_percent=False)
    category_bars(key, frac, pop_name, CC / key / "figures" / "length_category_fraction_bar.png")
    print(f"\n== {key}")
    print(pd.DataFrame(summ).to_string(index=False))
    print(tests.assign(p_value=tests["p_value"].map("{:.1e}".format)).to_string(index=False))
