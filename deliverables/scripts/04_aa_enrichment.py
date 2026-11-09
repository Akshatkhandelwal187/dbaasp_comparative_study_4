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
The figures are drawn from the table by 05_enrichment_figures.py.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

from cprofiler.aminoacid import AminoAcid
from cprofiler.fasta import Fasta
from cprofiler.profile import CompositionProfiler

DELIV = Path(__file__).resolve().parent.parent
QUERIES = {"Natural": DELIV / "fasta" / "amp_natural_cdhit90.fasta",
           "Synthetic": DELIV / "fasta" / "amp_synthetic_cdhit90.fasta"}
BACKGROUND = DELIV / "background" / "swissprot51_5k.fasta"
ITERATIONS = 10000                 # default of the tool
ALPHA = 0.05                       # default of the tool

# SwissProt entries whose protein name is an antimicrobial peptide / protein (name screen, used for the background check only)
AMP_NAME = re.compile(r"defensin|cecropin|magainin|antimicrobial|bacteriocin|cathelicidin|lantibiotic|\bnisin|microcin|attacin|"
                      r"protegrin|dermaseptin|brevinin|temporin|indolicidin|histatin|hepcidin|bombinin|gaegurin|esculentin|"
                      r"melittin|mastoparan|azurocidin|\bthionin|diptericin|drosocin|apidaecin|tachyplesin", re.I)

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


print(wide[wide["type"] == "amino acid"][["test_name", "freq_swissprot_%", "freq_Natural_%", "fold_Natural", "result_Natural",
                                          "freq_Synthetic_%", "fold_Synthetic", "result_Synthetic"]].round(2).to_string(index=False))
