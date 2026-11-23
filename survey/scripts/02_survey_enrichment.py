"""Survey step 2 - amino acid enrichment for every class with at least MIN_N peptides (run 01 first).

Same test as deliverables/scripts/04_aa_enrichment.py (the Composition Profiler 'discover' test, v2.0):
  effect  = (frequency in the group - frequency in the background) / frequency in the background (residues pooled)
  p-value = share of random relabelings of WHOLE sequences (group <-> background) whose |effect| is at least as large
            (two-sided; p = 0 means < 1/ITERATIONS)
  Bonferroni over the 20 amino acids: alpha = 0.05 / 20
Only difference to the tool call: the permutations are vectorised in numpy (the tool is not installed here), and the
group is drawn without the tool's off-by-one on the background slice. Backgrounds: SwissProt (5,002 proteins) and non-AMP peptides.
Writes survey/tables/survey_aa_enrichment.csv (long format).
"""
import numpy as np
import pandas as pd
from common import DELIV, TABLES, TARGET_GROUPS, MIN_N, ITERATIONS, ALPHA, SEED, AA_ORDER

AAS = AA_ORDER
BACKGROUNDS = {"SwissProt": DELIV / "background" / "swissprot51_5k.fasta",
               "Non-AMP": DELIV / "background" / "pepbench_nonamp_background.fasta"}


def counts_matrix(seqs):
    m = np.zeros((len(seqs), 20), dtype=np.float64)
    idx = {a: i for i, a in enumerate(AAS)}
    for r, s in enumerate(seqs):
        for ch in s.upper():
            j = idx.get(ch)
            if j is not None:                       # B, X, U, ... are skipped, like the tool does
                m[r, j] += 1
    return m


def read_fasta(path):
    seqs, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def effect(q_sum, b_sum):
    return (q_sum / q_sum.sum() - b_sum / b_sum.sum()) / (b_sum / b_sum.sum())


def permutation_test(q, b, rng):
    obs = effect(q.sum(0), b.sum(0))
    pool, nq, total = np.vstack([q, b]), len(q), q.sum(0) + b.sum(0)
    hits = np.zeros(20)
    for _ in range(ITERATIONS):
        qs = pool[rng.permutation(len(pool))[:nq]].sum(0)
        hits += np.abs(effect(qs, total - qs)) >= np.abs(obs)
    return obs, hits / ITERATIONS


df = pd.read_csv(TABLES / "survey_peptides.csv")
nat = df["ORIGIN"] == "Natural"
groups = []   # (grouping, group, boolean mask)
for k in ["Animalia", "Bacteria", "Plantae", "Fungi", "Archaea", "Protista", "Multiple", "Unclassified"]:
    groups.append(("Source organism (Natural)", k, nat & (df["KINGDOM"] == k)))
for k in ["Ribosomal", "Nonribosomal", "Synthetic"]:
    groups.append(("Synthesis type", k, df["SOURCE_CLASS"] == k))
for label in TARGET_GROUPS.values():
    groups.append(("Function (all)", label, df[label]))
    groups.append(("Function (Natural)", label, nat & df[label]))
groups.append(("Reference", "All peptides", pd.Series(True, index=df.index)))

bg_counts = {n: counts_matrix(read_fasta(p)) for n, p in BACKGROUNDS.items()}
rows = []
for grouping, name, mask in groups:
    n = int(mask.sum())
    if n < MIN_N:
        print(f"skip {grouping} / {name}: n = {n} < {MIN_N}")
        continue
    q = counts_matrix(df.loc[mask, "SEQUENCE"].tolist())
    for bname, b in bg_counts.items():
        rng = np.random.default_rng(SEED)
        obs, p = permutation_test(q, b, rng)
        fq, fb = 100 * q.sum(0) / q.sum(), 100 * b.sum(0) / b.sum()
        for i, a in enumerate(AAS):
            sig = p[i] < ALPHA
            rows.append({"grouping": grouping, "group": name, "n_peptides": n, "background": bname, "amino_acid": a,
                         "freq_group_%": fq[i], "freq_background_%": fb[i], "effect": obs[i], "pct_change": 100 * obs[i],
                         "pvalue": p[i], "result": ("Enriched" if obs[i] > 0 else "Depleted") if sig else "Not significant"})
    print(f"done {grouping} / {name} (n = {n})", flush=True)

out = pd.DataFrame(rows)
out = out.round({c: 4 for c in out.columns if out[c].dtype == float})
out.to_csv(TABLES / "survey_aa_enrichment.csv", index=False)
print(f"\n{len(out):,} rows written")
