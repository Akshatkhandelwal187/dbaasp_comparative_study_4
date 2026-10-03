"""Deliverable step 3 - checks on the CD-HIT result and two what-if runs (run 01_remove_redundancy.py first).

Everything is rebuilt from deliverables/tables, nothing is read from outside deliverables/.
Writes deliverables/tables/cdhit_checks.csv  (check, result)

Checks
  * every peptide CD-HIT removed has >= 90 % identity to a kept representative that is at least as long
  * clustering the final >= 10 aa set again removes nothing (the result is stable)
What-ifs (not used for the deliverables, only to show how much the choices matter)
  * Natural and Synthetic clustered separately (the old pipeline's way) instead of together
  * peptides of 5-9 aa clustered too (the deliverables leave peptides < 10 aa out of CD-HIT)
"""
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

DELIV = Path(__file__).resolve().parent.parent
meta = pd.read_csv(DELIV / "tables" / "peptide_metadata.csv")
rem = pd.read_csv(DELIV / "tables" / "removed_records.csv")
summary = pd.read_csv(DELIV / "tables" / "redundancy_removal_summary.csv", dtype={"step": str})
rem["SEQUENCE"] = rem["SEQUENCE"].str.upper()
cd_lost = rem[rem["STEP"].str.startswith("6")]
results = []


def cdhit(table, min_len):
    """cd-hit at 90 % on table (columns ID, SEQUENCE); returns the IDs it keeps. min_len = shortest length it may read."""
    with tempfile.TemporaryDirectory() as tmp:
        with open(f"{tmp}/in.fasta", "w") as f:
            for i, s in zip(table["ID"], table["SEQUENCE"]):
                f.write(f">{i}\n{s}\n")
        res = subprocess.run(["cd-hit", "-i", f"{tmp}/in.fasta", "-o", f"{tmp}/out.fasta", "-c", "0.9", "-n", "5",
                              "-d", "0", "-l", str(min_len - 1), "-M", "0", "-T", "1"],
                             capture_output=True, text=True, check=True)
        assert int(re.search(r"total seq:\s+(\d+)", res.stdout).group(1)) == len(table), "cd-hit dropped sequences"
        return {l[1:].strip() for l in open(f"{tmp}/out.fasta") if l[0] == ">"}


# ---- checks
length = dict(zip(meta["ID"], meta["LENGTH"]))
results.append(("CD-HIT removed peptides", len(cd_lost)))
results.append(("lowest identity to the kept representative (%)", cd_lost["IDENTITY_%"].min()))
results.append(("every representative is in the final FASTA", bool(cd_lost["COVERED_BY"].isin(meta["ID"]).all())))
results.append(("no representative is shorter than the peptide it covers",
                bool((cd_lost["COVERED_BY"].map(length) >= cd_lost["LENGTH"]).all())))
final_ge10 = meta[meta["LENGTH"] >= 10]
results.append(("re-clustering the final >= 10 aa set removes", len(final_ge10) - len(cdhit(final_ge10, 10))))

# ---- peptides as they were just before CD-HIT (final set + what CD-HIT removed)
before = pd.concat([meta[["ID", "ORIGIN", "SEQUENCE"]],
                    cd_lost.rename(columns={"DBAASP_ID": "ID"})[["ID", "ORIGIN", "SEQUENCE"]]], ignore_index=True)
before["LENGTH"] = before["SEQUENCE"].str.len()

# ---- what-if 1: cluster each origin on its own (needs the sets before the Natural/Synthetic duplicate step)
same_seq_dups = rem[rem["STEP"].str.startswith("5")].rename(columns={"DBAASP_ID": "ID"})[["ID", "ORIGIN", "SEQUENCE"]]
clean = pd.concat([before, same_seq_dups.assign(LENGTH=same_seq_dups["SEQUENCE"].str.len())], ignore_index=True)
assert len(clean) == int(summary.loc[summary["step"] == "4", "remaining_Total"].iloc[0])
per_origin = {}
for origin in ["Natural", "Synthetic"]:
    sub = clean[clean["ORIGIN"] == origin]
    kept = cdhit(sub[sub["LENGTH"] >= 10], 10) | set(sub.loc[sub["LENGTH"] < 10, "ID"])
    per_origin[origin] = len(kept)
    median = sub.loc[sub["ID"].isin(kept), "LENGTH"].median()
    results.append((f"what-if, {origin} clustered on its own: kept (of {len(sub):,}) / median length",
                    f"{len(kept)} / {median:g}"))
results.append(("what-if, separate clustering: total (still contains Natural/Synthetic near-duplicates)",
                sum(per_origin.values())))
medians = meta.groupby("ORIGIN")["LENGTH"].median()
results.append(("combined clustering (deliverable): Natural / Synthetic / total; median length Natural / Synthetic",
                f"{(meta['ORIGIN'] == 'Natural').sum()} / {(meta['ORIGIN'] == 'Synthetic').sum()} / {len(meta)}; "
                f"{medians['Natural']:g} / {medians['Synthetic']:g}"))
cross_origin = pd.crosstab(cd_lost["ORIGIN"], cd_lost["COVERED_BY"].map(dict(zip(before["ID"], before["ORIGIN"]))))
results.append(("CD-HIT removals: Synthetic peptide covered by a Natural one", int(cross_origin.loc["Synthetic", "Natural"])))
results.append(("CD-HIT removals: Natural peptide covered by a Synthetic one", int(cross_origin.loc["Natural", "Synthetic"])))

# ---- what-if 2: also cluster the 5-9 aa peptides (cd-hit cannot read sequences shorter than its word size, 5)
ge5 = before[before["LENGTH"] >= 5]
kept = cdhit(ge5, 5)
extra = ge5[(ge5["LENGTH"] < 10) & ~ge5["ID"].isin(kept)]
results.append(("what-if, 5-9 aa peptides clustered too: additionally removed", len(extra)))
results.append(("   of which Natural / Synthetic",
                f"{(extra['ORIGIN'] == 'Natural').sum()} / {(extra['ORIGIN'] == 'Synthetic').sum()}"))
results.append(("peptides of 5-9 aa / 2-4 aa / 1 aa in the final FASTA",
                f"{meta['LENGTH'].between(5, 9).sum()} / {meta['LENGTH'].between(2, 4).sum()} / {(meta['LENGTH'] == 1).sum()}"))

out = pd.DataFrame(results, columns=["check", "result"])
out.to_csv(DELIV / "tables" / "cdhit_checks.csv", index=False)
print(out.to_string(index=False))
