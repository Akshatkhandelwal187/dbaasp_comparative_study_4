"""Deliverable step 4 - amino acid enrichment of the Natural and Synthetic peptides against SwissProt
(run 01_remove_redundancy.py first).

Uses the official Composition Profiler (github.com/vvacic/cprofiler, v2.0) with the same functions and settings as
scripts/06_run_composition_profiler.py, i.e. what
    cprof discover -Q <query.fasta> -B <background.fasta> -I 10000 -A 0.05 -b
does, but keeps the table:
  * effect = (frequency in the peptides - frequency in SwissProt) / frequency in SwissProt (residues pooled over all sequences)
  * p-value = share of 10,000 random relabelings of whole sequences (query <-> background) with an effect at least as large
  * Bonferroni: alpha = 0.05 / 40 tests (20 amino acids + 20 amino acid groups)
Setup (once):  git clone https://github.com/vvacic/cprofiler  &&  pip install -e cprofiler

Query       fasta/amp_natural_cdhit90.fasta, fasta/amp_synthetic_cdhit90.fasta
Background  background/swissprot51_5k.fasta (5,002 proteins; byte-identical to the tool's own `-D sprot` file)
Writes      tables/aa_enrichment_natural_vs_synthetic.csv   all 40 tests, both datasets
            tables/aa_enrichment_background_check.csv       same effects with AMP-named proteins removed from the background
            figures/amino_acid_enrichment.png               the 20 amino acids
            figures/amino_acid_group_enrichment.png         the 20 amino acid groups of the tool
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path as MplPath

from cprofiler.aminoacid import AminoAcid
from cprofiler.fasta import Fasta
from cprofiler.profile import CompositionProfiler

DELIV = Path(__file__).resolve().parent.parent
QUERIES = {"Natural": DELIV / "fasta" / "amp_natural_cdhit90.fasta",
           "Synthetic": DELIV / "fasta" / "amp_synthetic_cdhit90.fasta"}
BACKGROUND = DELIV / "background" / "swissprot51_5k.fasta"
ITERATIONS = 10000                 # default of the tool
ALPHA = 0.05                       # default of the tool
DPI = 200

# SwissProt entries whose protein name is an antimicrobial peptide / protein (name screen, used for the background check only)
AMP_NAME = re.compile(r"defensin|cecropin|magainin|antimicrobial|bacteriocin|cathelicidin|lantibiotic|\bnisin|microcin|attacin|"
                      r"protegrin|dermaseptin|brevinin|temporin|indolicidin|histatin|hepcidin|bombinin|gaegurin|esculentin|"
                      r"melittin|mastoparan|azurocidin|\bthionin|diptericin|drosocin|apidaecin|tachyplesin", re.I)

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
NATURAL, SYNTHETIC = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "figure.dpi": DPI, "savefig.dpi": DPI,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "text.color": INK, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.edgecolor": AXIS, "axes.linewidth": 0.75, "font.size": 10})

alphabet = AminoAcid.get_order("alpha")                  # A C D E F G H I K L M N P Q R S T V W Y
groups = AminoAcid.get_groups()
group_names = AminoAcid.get_group_names()
alpha = ALPHA / (len(alphabet) + len(groups))            # Bonferroni correction (the -b option)
residues_of = {**{a: a for a in alphabet}, **{group_names[k]: groups[k] for k in groups}}


# ------------------------------------------------------------------ enrichment
def read_counts(sequences):
    return Fasta.count_chars(sequences, alphabet)         # only the 20 letters are counted, X etc. are skipped


def pooled_frequencies(counts):
    """% of all counted residues for each of the 40 tests (the tool's own pooling, written out again as a check)."""
    total = counts.sum()
    by_group = [counts[:, [alphabet.index(x) for x in groups[k]]].sum() for k in groups]
    return 100 * np.r_[counts.sum(axis=0), by_group] / total


def enrichment(query_counts, background_counts):
    np.random.seed(1)                                     # the p-values come from random sampling, so fix the seed
    table = CompositionProfiler.discover(query_counts, background_counts, alphabet=alphabet, groups=groups,
                                         group_names=group_names, iterations=ITERATIONS, alpha_value=alpha)
    fq, fb = pooled_frequencies(query_counts), pooled_frequencies(background_counts)
    assert np.allclose(table["effect"], (fq - fb) / fb), "tool and own frequency calculation disagree"
    table["freq_peptides_%"], table["freq_swissprot_%"] = fq, fb
    return table


background_seqs = Fasta.read(str(BACKGROUND))
amp_named = [s.header for s in background_seqs if AMP_NAME.search(s.header)]
background_counts = read_counts(background_seqs)
background_counts_no_amp = read_counts([s for s in background_seqs if not AMP_NAME.search(s.header)])

main, check = {}, {}
for name, path in QUERIES.items():
    query_counts = read_counts(Fasta.read(str(path)))
    main[name] = enrichment(query_counts, background_counts)
    check[name] = enrichment(query_counts, background_counts_no_amp)
    print(f"{name}: {len(query_counts):,} peptides, {int(query_counts.sum()):,} residues counted")
print(f"SwissProt background: {len(background_counts):,} proteins, {int(background_counts.sum()):,} residues; "
      f"{len(amp_named)} with an AMP-type name")

# ------------------------------------------------------------------ tables
wide = pd.DataFrame({"test_name": main["Natural"]["test_name"],
                     "type": ["amino acid"] * len(alphabet) + ["group"] * len(groups),
                     "residues": [residues_of[t] for t in main["Natural"]["test_name"]],
                     "freq_swissprot_%": main["Natural"]["freq_swissprot_%"]})
for name, table in main.items():
    wide[f"freq_{name}_%"] = table["freq_peptides_%"]
    wide[f"effect_{name}"] = table["effect"]
    wide[f"fold_{name}"] = 1 + table["effect"]
    wide[f"pvalue_{name}"] = table["pvalue"]
    wide[f"result_{name}"] = table["test_result"]
wide.round({c: 4 for c in wide.columns if wide[c].dtype == float}).to_csv(
    DELIV / "tables" / "aa_enrichment_natural_vs_synthetic.csv", index=False)

bg = pd.DataFrame({"test_name": wide["test_name"]})
for name in main:
    bg[f"effect_{name}"] = main[name]["effect"]
    bg[f"effect_{name}_without_AMP_named"] = check[name]["effect"]
    bg[f"result_{name}"] = main[name]["test_result"]
    bg[f"result_{name}_without_AMP_named"] = check[name]["test_result"]
bg.round({c: 4 for c in bg.columns if bg[c].dtype == float}).to_csv(
    DELIV / "tables" / "aa_enrichment_background_check.csv", index=False)
for name in main:
    delta = (check[name]["effect"] - main[name]["effect"]).abs()
    changed = bg.loc[bg[f"result_{name}"] != bg[f"result_{name}_without_AMP_named"], "test_name"].tolist()
    print(f"background check, {name}: largest change of an effect {delta.max():.3f}; results that change: {changed or 'none'}")


# ------------------------------------------------------------------ figures
def rounded_bar(ax, u0, u1, value, color, alpha=1.0, horizontal=False, radius_pt=3):
    """Bar from 0 to value (either sign) with a ~4 px rounded data end and a square base.
    u0..u1 is the bar's extent across its length; the radius is in points, not data units."""
    if value == 0:
        return
    inv, px = ax.transData.inverted(), ax.figure.dpi / 72
    origin = inv.transform((0, 0))
    per_pt_x, per_pt_y = abs(inv.transform((px, 0))[0] - origin[0]), abs(inv.transform((0, px))[1] - origin[1])
    ru, rv = (per_pt_y, per_pt_x) if horizontal else (per_pt_x, per_pt_y)
    ru, rv = radius_pt * ru, min(radius_pt * rv, abs(value))
    h, t = abs(value), np.linspace(0, np.pi / 2, 10)
    left = np.c_[u0 + ru - ru * np.cos(t), h - rv + rv * np.sin(t)]
    right = np.c_[u1 - ru + ru * np.sin(t), h - rv + rv * np.cos(t)]
    uv = np.vstack([[u0, 0], left, right, [u1, 0], [u0, 0]])    # last vertex = CLOSEPOLY marker
    uv[:, 1] *= np.sign(value)
    ax.add_patch(PathPatch(MplPath(uv[:, ::-1] if horizontal else uv, closed=True), facecolor=color, alpha=alpha,
                           linewidth=0, zorder=3))


def fold_ticks(ax, horizontal, lo=-2, hi=2):
    """Axis in log2 of the fold change, labelled as 0.25x ... 4x."""
    ticks = list(range(lo, hi + 1))
    labels = [f"{2.0 ** t:g}×" for t in ticks]
    (ax.set_xticks if horizontal else ax.set_yticks)(ticks)
    (ax.set_xticklabels if horizontal else ax.set_yticklabels)(labels)


def style_axes(ax, horizontal=False):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x" if horizontal else "y", color=GRID, linewidth=0.75)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=9)


def draw_bars(ax, rows, horizontal):
    """Two bars (Natural, Synthetic) per row of `rows`; faded when not significant after Bonferroni."""
    ax.figure.canvas.draw()                                   # layout is fixed, so points can be turned into data units
    inv = ax.transData.inverted()
    px = ax.figure.dpi / 72
    one_pt = abs(inv.transform((0, px))[1] - inv.transform((0, 0))[1]) if horizontal else \
        abs(inv.transform((px, 0))[0] - inv.transform((0, 0))[0])
    width, gap = min(16 * one_pt, 0.42), 1.5 * one_pt           # bars <= 24 px thick, 2 px apart
    for i, (_, row) in enumerate(rows.iterrows()):
        for side, (name, color) in enumerate([("Natural", NATURAL), ("Synthetic", SYNTHETIC)]):
            u0 = i - gap / 2 - width if side == 0 else i + gap / 2
            significant = row[f"result_{name}"] != "Not significant"
            rounded_bar(ax, u0, u0 + width, np.log2(row[f"fold_{name}"]), color, alpha=1.0 if significant else 0.35,
                        horizontal=horizontal)


def series_key(ax, x, y, name, color, transform, ha):
    """A small square in the series colour followed by the series name (the text itself stays in ink)."""
    ax.plot([x], [y], marker="s", markersize=6, color=color, transform=transform, clip_on=False, linestyle="")
    ax.annotate(name, (x, y), xycoords=transform, xytext=(8 if ha == "left" else -8, 0), textcoords="offset points",
                ha=ha, va="center", fontsize=8.5, color=INK2, annotation_clip=False)


def legend(ax, faded, **kwargs):
    handles = [Patch(facecolor=NATURAL), Patch(facecolor=SYNTHETIC)] + ([Patch(facecolor=MUTED, alpha=0.35)] if faded else [])
    n = {k: len(QUERIES_SEQ[k]) for k in QUERIES}
    labels = [f"Natural  (n = {n['Natural']:,})", f"Synthetic  (n = {n['Synthetic']:,})"] + \
        (["faded: not significant"] if faded else [])
    ax.legend(handles, labels, frameon=False, fontsize=9.5, labelcolor=INK2, handlelength=1.0, handleheight=1.0,
              handletextpad=0.6, columnspacing=2.0, **kwargs)


def heading(fig, title, subtitle, x=0.07):
    fig.text(x, 0.955, title, fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(x, 0.905, subtitle, fontsize=9.5, color=INK2, ha="left", va="top", linespacing=1.5)


def fold_text(ax, x, y, row, name, transform):
    """The fold change as a number; muted when the result is not significant."""
    significant = row[f"result_{name}"] != "Not significant"
    ax.text(x, y, f"{row[f'fold_{name}']:.2f}", transform=transform, ha="center", va="center", fontsize=8.5,
            color=INK2 if significant else MUTED, clip_on=False)


QUERIES_SEQ = {name: Fasta.read(str(path)) for name, path in QUERIES.items()}
NOTE = (f"Residue frequency in the peptides relative to SwissProt ({len(background_seqs):,} proteins), log scale. "
        "Peptides after CD-HIT (90 % identity).")
FOOT = (f"Numbers under the axis: frequency in the peptides / frequency in SwissProt (grey = not significant).\n"
        f"Significance: permutation test, {ITERATIONS:,} iterations, Bonferroni (alpha = {ALPHA} / {len(alphabet) + len(groups)}).")

# ---- 20 amino acids, alphabetical, with the numbers in a small table under the axis
aa = wide[wide["type"] == "amino acid"].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11, 5.6))
fig.subplots_adjust(left=0.085, right=0.985, top=0.8, bottom=0.25)
ax.set_xlim(-0.6, len(aa) - 0.4)
ax.set_ylim(-2.4, 2.7)
fold_ticks(ax, horizontal=False)
style_axes(ax)
draw_bars(ax, aa, horizontal=False)
ax.axhline(0, color=AXIS, linewidth=0.75, zorder=2)
ax.set_xticks(range(len(aa)))
ax.set_xticklabels(aa["test_name"], color=INK2, fontsize=10.5)
ax.set_ylabel("Frequency relative to SwissProt", labelpad=8)
below = ax.get_xaxis_transform()                                   # x in data units, y as a fraction of the axes
for name, color, y in (("Natural", NATURAL, -0.125), ("Synthetic", SYNTHETIC, -0.19)):
    series_key(ax, -0.012, y, name, color, ax.transAxes, "right")
    for i, (_, row) in enumerate(aa.iterrows()):
        fold_text(ax, i, y, row, name, below)
legend(ax, faded=True, loc="lower left", bbox_to_anchor=(0, 1.01), ncol=3)
fig.text(0.07, 0.025, FOOT, fontsize=8.5, color=MUTED, ha="left", va="bottom", linespacing=1.5)
heading(fig, "Amino acid enrichment in Natural and Synthetic peptides", NOTE)
fig.savefig(DELIV / "figures" / "amino_acid_enrichment.png")
plt.close(fig)

# ---- 20 amino acid groups of the tool, most enriched at the top, numbers in two columns on the right
FIX = {"propencity": "propensity", "potenential": "potential"}          # typos in the tool's group names
grp = wide[wide["type"] == "group"].copy()
grp["mean_log2"] = (np.log2(grp["fold_Natural"]) + np.log2(grp["fold_Synthetic"])) / 2
grp = grp.sort_values("mean_log2", ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(11, 8.4))
fig.subplots_adjust(left=0.40, right=0.86, top=0.85, bottom=0.075)
ax.set_ylim(len(grp) - 0.4, -0.6)                                       # first row on top
ax.set_xlim(-2.1, 1.6)
fold_ticks(ax, horizontal=True, hi=1)
style_axes(ax, horizontal=True)
draw_bars(ax, grp, horizontal=True)
ax.axvline(0, color=AXIS, linewidth=0.75, zorder=2)
ax.set_yticks(range(len(grp)))
labels = [re.sub("|".join(FIX), lambda m: FIX[m.group()], t) + f"  ({r})" for t, r in zip(grp["test_name"], grp["residues"])]
ax.set_yticklabels(labels, color=INK2, fontsize=9)
ax.set_xlabel("Frequency relative to SwissProt", labelpad=8)
side = ax.get_yaxis_transform()                                          # x as a fraction of the axes, y in data units
for name, color, x in (("Natural", NATURAL, 1.06), ("Synthetic", SYNTHETIC, 1.19)):
    series_key(ax, x - 0.045, -1.15, name, color, side, "left")
    for i, (_, row) in enumerate(grp.iterrows()):
        fold_text(ax, x, i, row, name, side)
legend(ax, faded=bool((grp[["result_Natural", "result_Synthetic"]] == "Not significant").any().any()),
       loc="lower right", bbox_to_anchor=(1.0, 0.0))
heading(fig, "Amino acid group enrichment in Natural and Synthetic peptides",
        "Groups of the Composition Profiler (residues in brackets), sorted by the mean of both datasets.\n"
        "Numbers on the right: frequency in the peptides / frequency in SwissProt (grey = not significant).", x=0.02)
fig.savefig(DELIV / "figures" / "amino_acid_group_enrichment.png")
plt.close(fig)

print(wide[wide["type"] == "amino acid"][["test_name", "freq_swissprot_%", "freq_Natural_%", "fold_Natural", "result_Natural",
                                          "freq_Synthetic_%", "fold_Synthetic", "result_Synthetic"]].round(2).to_string(index=False))
