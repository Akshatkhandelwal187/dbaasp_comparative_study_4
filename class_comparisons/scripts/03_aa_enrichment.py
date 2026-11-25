"""Class comparisons step 3 - amino acid enrichment of every tested class against two backgrounds (run 01 first).

Same test and settings as deliverables/scripts/04_aa_enrichment.py, with the official Composition Profiler
(github.com/vvacic/cprofiler, v2.0; the functions behind `cprof discover`):
  * effect  = (frequency in the class - frequency in the background) / frequency in the background (residues pooled)
  * p-value = share of 10,000 random relabelings of whole sequences (class <-> background) at least as extreme (two-tailed)
  * Bonferroni over the 20 amino acids: alpha = 0.05 / 20 = 0.0025; np.random.seed(1) before every test
Backgrounds: SwissProt (5,002 proteins) and non-AMP peptides (26,635), the files in deliverables/background/.
Writes <comparison>/tables/aa_enrichment.csv (long format) for every comparison.
Usage: python 03_aa_enrichment.py [comparison ...]      (default: all; about 15-20 minutes in total)
"""
import sys

import numpy as np
import pandas as pd
from cprofiler.aminoacid import AminoAcid
from cprofiler.fasta import Fasta
from cprofiler.profile import CompositionProfiler

from common import ALPHA, BACKGROUNDS, CC, COMPARISONS, DELIV, ITERATIONS, load_table, slug, tested

alphabet = AminoAcid.get_order("alpha")
alpha = ALPHA / len(alphabet)


def counts(path):
    return Fasta.count_chars(Fasta.read(str(path)), alphabet)


df = load_table()
bg = {name: counts(path) for name, path in BACKGROUNDS.items()}
keys = sys.argv[1:] or list(COMPARISONS)
for key in keys:
    rows = []
    cl, _ = tested(df, key)
    for name, mask in cl:
        q = counts(CC / key / "fasta" / f"{slug(name)}.fasta")
        assert len(q) == int(mask.sum())
        for bname, b in bg.items():
            np.random.seed(1)
            table = CompositionProfiler.discover(q, b, alphabet=alphabet, groups={}, group_names={},
                                                 iterations=ITERATIONS, alpha_value=alpha)
            fq, fb = 100 * q.sum(axis=0) / q.sum(), 100 * b.sum(axis=0) / b.sum()
            assert np.allclose(table["effect"], (fq - fb) / fb) and list(table["test_name"]) == list(alphabet)
            eff = table["effect"].values
            sig = table["pvalue"].values < alpha
            rows.append(pd.DataFrame({"comparison": key, "class": name, "n_peptides": len(q), "n_residues": int(q.sum()),
                                      "background": bname, "amino_acid": list(alphabet), "freq_class_%": fq,
                                      "freq_background_%": fb, "effect": eff, "pct_change": 100 * eff,
                                      "pvalue": table["pvalue"].values,
                                      "result": np.where(sig, np.where(eff > 0, "Enriched", "Depleted"), "Not significant")}))
        print(f"{key} / {name}: n = {len(q):,}", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.round({c: 4 for c in out.columns if out[c].dtype == float}).to_csv(CC / key / "tables" / "aa_enrichment.csv", index=False)

# consistency check: the Synthetic class must reproduce deliverables/tables/aa_enrichment.csv exactly
if "synthesis_type" in keys:
    mine = pd.read_csv(CC / "synthesis_type" / "tables" / "aa_enrichment.csv").query("`class` == 'Synthetic'")
    ref = pd.read_csv(DELIV / "tables" / "aa_enrichment.csv").query("origin == 'Synthetic'")
    m = mine.merge(ref, on=["background", "amino_acid"])
    assert len(m) == 40 and np.allclose(m["effect_x"], m["effect_y"]) and np.allclose(m["pvalue_x"], m["pvalue_y"])
    assert (m["result_x"].str[:3] == m["result_y"].str[:3]).all()
    print("check passed: Synthetic reproduces deliverables/tables/aa_enrichment.csv (effects, p-values, results)")
