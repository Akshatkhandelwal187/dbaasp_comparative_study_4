"""Deliverable step 4 - amino acid enrichment of the Natural and Synthetic peptides against TWO backgrounds
(run 01_remove_redundancy.py first).

Uses the official Composition Profiler (github.com/vvacic/cprofiler, v2.0): the same functions as
    cprof discover -Q <query.fasta> -B <background.fasta> -I 10000 -A 0.05 -b
but only the 20 single amino acids are tested (no amino acid groups) and the table is kept.
  * effect = (frequency in the peptides - frequency in the background) / frequency in the background (residues pooled)
  * p-value = share of 10,000 random relabelings of whole sequences (query <-> background) with an effect at least as large
  * Bonferroni: alpha = 0.05 / 20 tests = 0.0025
Setup (once):  git clone https://github.com/vvacic/cprofiler  &&  pip install -e cprofiler

Queries      fasta/amp_natural_cdhit90.fasta, fasta/amp_synthetic_cdhit90.fasta
Backgrounds  SwissProt     background/swissprot51_5k.fasta            5,002 normal proteins (the tool's own `-D sprot` file)
             Non-AMP       background/pepbench_nonamp_background.fasta  26,635 non-antimicrobial peptides (PepBench)
Writes       tables/aa_enrichment.csv   20 amino acids x 2 origins x 2 backgrounds (long format)
The figure is drawn from that table by 05_enrichment_figures.py.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from cprofiler.aminoacid import AminoAcid
from cprofiler.fasta import Fasta
from cprofiler.profile import CompositionProfiler

DELIV = Path(__file__).resolve().parent.parent
QUERIES = {"Natural": DELIV / "fasta" / "amp_natural_cdhit90.fasta",
           "Synthetic": DELIV / "fasta" / "amp_synthetic_cdhit90.fasta"}
BACKGROUNDS = {"SwissProt": DELIV / "background" / "swissprot51_5k.fasta",
               "Non-AMP": DELIV / "background" / "pepbench_nonamp_background.fasta"}
ITERATIONS = 10000                 # default of the tool
ALPHA = 0.05                       # default of the tool

alphabet = AminoAcid.get_order("alpha")                  # A C D E F G H I K L M N P Q R S T V W Y
alpha = ALPHA / len(alphabet)                            # Bonferroni correction over the 20 tests


def read_counts(path):
    seqs = Fasta.read(str(path))
    return Fasta.count_chars(seqs, alphabet)             # one row per sequence, one column per amino acid


rows = []
for bname, bpath in BACKGROUNDS.items():
    bcounts = read_counts(bpath)
    print(f"background {bname}: {len(bcounts):,} sequences, {int(bcounts.sum()):,} residues")
    for qname, qpath in QUERIES.items():
        qcounts = read_counts(qpath)
        np.random.seed(1)                                # the p-values come from random sampling, so fix the seed
        table = CompositionProfiler.discover(qcounts, bcounts, alphabet=alphabet, groups={}, group_names={},
                                             iterations=ITERATIONS, alpha_value=alpha)
        fq, fb = 100 * qcounts.sum(axis=0) / qcounts.sum(), 100 * bcounts.sum(axis=0) / bcounts.sum()
        assert np.allclose(table["effect"], (fq - fb) / fb), "tool and own frequency calculation disagree"
        assert list(table["test_name"]) == list(alphabet)
        rows.append(pd.DataFrame({"background": bname, "origin": qname, "amino_acid": list(alphabet),
                                  "freq_peptides_%": fq, "freq_background_%": fb,
                                  "effect": table["effect"].values, "fold": 1 + table["effect"].values,
                                  "pct_change": 100 * table["effect"].values, "pvalue": table["pvalue"].values,
                                  "result": table["test_result"].values}))
        print(f"  {qname}: {len(qcounts):,} peptides, {int(qcounts.sum()):,} residues")

out = pd.concat(rows, ignore_index=True)
out.round({c: 4 for c in out.columns if out[c].dtype == float}).to_csv(DELIV / "tables" / "aa_enrichment.csv", index=False)
wide = out.pivot_table(index="amino_acid", columns=["background", "origin"], values="pct_change").round(0)
print((wide.astype(int)).to_string())
