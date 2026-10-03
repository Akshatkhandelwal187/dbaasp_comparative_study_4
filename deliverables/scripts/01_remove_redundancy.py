"""Deliverable step 1 - one non-redundant AMP dataset (Natural + Synthetic) with an exact loss count.

Starts from the untouched raw data (raw_data/dbaasp_datasets_corrected.zip, by_origin/*.csv) and nothing in the
repository outside deliverables/ is modified.

  1. cleaning, same four rules as scripts/01_clean_datasets.py (per origin)
  2. identical sequence in Natural AND Synthetic -> keep the Natural entry
  3. CD-HIT, >= 90 % identity, on the combined set, so the final FASTA is non-redundant as a whole
       cd-hit -c 0.9 -n 5 -d 0 -l 9      (peptides < 10 aa skip CD-HIT, as in scripts/02_run_cdhit.py)
     NOTE -l 9: cd-hit's default (-l 10) silently throws away every sequence of 10 aa or fewer.
     scripts/02_run_cdhit.py sent the 10-aa peptides to cd-hit with the default, so they were dropped and
     counted as "redundant". Here the number of sequences cd-hit actually read is checked against the number sent.

Writes (all under deliverables/):
  fasta/amp_nonredundant.fasta          >DBAASP_<ID>  one-line upper-case sequence
  tables/peptide_metadata.csv           ID -> origin, length, length class, sequence
  tables/removed_records.csv            every record lost, why, and which kept entry covers it
  tables/redundancy_removal_summary.csv the step-by-step counts (Natural / Synthetic / Total)
  cdhit/                                cd-hit input, output, cluster file and log
"""
import io
import re
import subprocess
import zipfile
from pathlib import Path

import pandas as pd

DELIV = Path(__file__).resolve().parent.parent
RAW_ZIP = DELIV.parent / "raw_data" / "dbaasp_datasets_corrected.zip"
CDHIT_DIR = DELIV / "cdhit"
IDENTITY, MIN_LEN = 0.9, 10          # peptides shorter than MIN_LEN are not clustered
ORIGINS = ["Natural", "Synthetic"]   # Natural first: it wins when the same sequence is in both

# length classes exactly as the Length_Class column of the data defines them
CLASSES = [(1, 1, "<2 aa (out of range)"), (2, 9, "Ultrashort (2-9 aa)"), (10, 24, "Short (10-24 aa)"),
           (25, 50, "Medium (25-50 aa)"), (51, 100, "Long (51-100 aa)"), (101, 10**6, ">100 aa")]


def length_class(n):
    return next(name for lo, hi, name in CLASSES if lo <= n <= hi)


removed = []   # one dict per lost record


def log_removed(df, step, origin, kept_id=None):
    for i, row in df.iterrows():
        removed.append({"DBAASP_ID": f"DBAASP_{row['ID']}", "ORIGIN": origin, "SEQUENCE": row["SEQUENCE"],
                        "LENGTH": len(row["SEQUENCE"]), "STEP": step,
                        "COVERED_BY": None if kept_id is None else kept_id(row), "IDENTITY_%": None})


# ---------------------------------------------------------------- 1. cleaning
zf = zipfile.ZipFile(RAW_ZIP)
clean, counts = {}, {}
for origin in ORIGINS:
    df = pd.read_csv(io.BytesIO(zf.read(f"by_origin/{origin}.csv")))
    df["SEQUENCE"] = df["SEQUENCE"].astype(str).str.strip()
    c = {"raw rows": len(df)}

    dup_row = df.duplicated()                                   # same record entered twice
    log_removed(df[dup_row], "1 exact duplicate row", origin, lambda r: f"DBAASP_{r['ID']}")
    df = df[~dup_row]
    c["after duplicate rows"] = len(df)

    multi = df["SEQUENCE"].str.contains(r"\s")                  # several chains in one cell
    log_removed(df[multi], "2 multi-chain entry", origin)
    df = df[~multi]
    c["after multi-chain"] = len(df)

    all_x = df["SEQUENCE"].str.upper().str.strip("X") == ""     # no real amino acid
    log_removed(df[all_x], "3 sequence of only X", origin)
    df = df[~all_x]
    c["after all-X"] = len(df)

    key = df["SEQUENCE"].str.upper()                            # D-amino acids are lowercase
    dup_seq = key.duplicated()
    first = dict(zip(key[~dup_seq], df.loc[~dup_seq, "ID"]))
    log_removed(df[dup_seq], "4 duplicate sequence (same origin)", origin,
                lambda r: f"DBAASP_{first[r['SEQUENCE'].upper()]}")
    df = df[~dup_seq].copy()
    c["after duplicate sequences"] = len(df)

    df["ORIGIN"] = origin
    assert df["ID"].is_unique and (df["SEQUENCE"].str.len() == df["Length_of_Sequence"]).all()
    clean[origin], counts[origin] = df, c

data = pd.concat([clean[o] for o in ORIGINS], ignore_index=True)
data["KEY"] = data["SEQUENCE"].str.upper()
assert data["ID"].is_unique, "IDs overlap between Natural and Synthetic"

# ---------------------------------------------------------------- 2. same sequence in Natural and Synthetic
dup_cross = data["KEY"].duplicated()          # Natural rows come first, so the Synthetic copy is the one dropped
first = dict(zip(data.loc[~dup_cross, "KEY"], data.loc[~dup_cross, "ID"]))
for _, r in data[dup_cross].iterrows():
    removed.append({"DBAASP_ID": f"DBAASP_{r['ID']}", "ORIGIN": r["ORIGIN"], "SEQUENCE": r["SEQUENCE"],
                    "LENGTH": len(r["SEQUENCE"]), "STEP": "5 identical sequence in Natural and Synthetic",
                    "COVERED_BY": f"DBAASP_{first[r['KEY']]}", "IDENTITY_%": None})
data = data[~dup_cross].copy()
data["HEADER"] = "DBAASP_" + data["ID"].astype(str)

# ---------------------------------------------------------------- 3. CD-HIT on the combined set
CDHIT_DIR.mkdir(exist_ok=True)
to_cluster = data[data["KEY"].str.len() >= MIN_LEN]
too_short = data[data["KEY"].str.len() < MIN_LEN]
with open(CDHIT_DIR / "input_ge10aa.fasta", "w") as f:
    for h, s in zip(to_cluster["HEADER"], to_cluster["KEY"]):
        f.write(f">{h}\n{s}\n")

cmd = ["cd-hit", "-i", str(CDHIT_DIR / "input_ge10aa.fasta"), "-o", str(CDHIT_DIR / "cdhit_ge10aa.fasta"),
       "-c", str(IDENTITY), "-n", "5", "-d", "0", "-l", str(MIN_LEN - 1), "-M", "0", "-T", "1"]
res = subprocess.run(cmd, capture_output=True, text=True, check=True)
(CDHIT_DIR / "cdhit.log").write_text(res.stdout)
n_read = int(re.search(r"total seq:\s+(\d+)", res.stdout).group(1))
assert n_read == len(to_cluster), f"cd-hit read {n_read} of {len(to_cluster)} sequences - some were discarded"

# cluster file: every non-representative line ends in "at <identity>%"
rep_of, identity_of, cluster = {}, {}, []
for line in open(CDHIT_DIR / "cdhit_ge10aa.fasta.clstr"):
    if line.startswith(">Cluster"):
        if cluster:
            rep = next(m for m in cluster if m[1] is None)[0]
            for name, pid in cluster:
                if pid is not None:
                    rep_of[name], identity_of[name] = rep, pid
        cluster = []
    else:
        m = re.search(r">(\S+?)\.\.\. (\*|at [^%]*%)", line)
        name = m.group(1)
        cluster.append((name, None if m.group(2) == "*" else float(re.search(r"([\d.]+)%", m.group(2)).group(1))))
rep = next(m for m in cluster if m[1] is None)[0]
for name, pid in cluster:
    if pid is not None:
        rep_of[name], identity_of[name] = rep, pid

cd_kept = {l[1:].strip() for l in open(CDHIT_DIR / "cdhit_ge10aa.fasta") if l.startswith(">")}
cd_lost = to_cluster[~to_cluster["HEADER"].isin(cd_kept)]
assert len(cd_lost) == len(rep_of) and len(cd_kept) + len(cd_lost) == len(to_cluster)
origin_of = dict(zip(data["HEADER"], data["ORIGIN"]))
for _, r in cd_lost.iterrows():
    removed.append({"DBAASP_ID": r["HEADER"], "ORIGIN": r["ORIGIN"], "SEQUENCE": r["SEQUENCE"],
                    "LENGTH": len(r["SEQUENCE"]), "STEP": "6 CD-HIT >= 90 % identity",
                    "COVERED_BY": rep_of[r["HEADER"]], "IDENTITY_%": identity_of[r["HEADER"]]})

final = pd.concat([too_short, to_cluster[to_cluster["HEADER"].isin(cd_kept)]]).copy()
final["ID"] = final["ID"].astype(int)
final = final.sort_values("ID").reset_index(drop=True)
assert final["HEADER"].is_unique and final["KEY"].is_unique

# ---------------------------------------------------------------- outputs
with open(DELIV / "fasta" / "amp_nonredundant.fasta", "w") as f:
    for h, s in zip(final["HEADER"], final["KEY"]):
        f.write(f">{h}\n{s}\n")

meta = pd.DataFrame({"ID": final["HEADER"], "ORIGIN": final["ORIGIN"], "LENGTH": final["KEY"].str.len(),
                     "LENGTH_CLASS": final["KEY"].str.len().map(length_class), "SEQUENCE": final["KEY"]})
meta.to_csv(DELIV / "tables" / "peptide_metadata.csv", index=False)

rem = pd.DataFrame(removed).sort_values(["STEP", "DBAASP_ID"]).reset_index(drop=True)
rem.to_csv(DELIV / "tables" / "removed_records.csv", index=False)

# step-by-step summary (steps are numbered in the STEP text, so sorting by it gives the processing order)
left = {o: counts[o]["raw rows"] for o in ORIGINS}
rows = [["0", "raw rows (DBAASP export)", None, None, None, *left.values(), sum(left.values())]]
for step, group in rem.groupby("STEP"):
    lost = group.groupby("ORIGIN").size().reindex(ORIGINS, fill_value=0)
    left = {o: left[o] - int(lost[o]) for o in ORIGINS}
    rows.append([step[0], step[2:], *lost.tolist(), int(lost.sum()), *left.values(), sum(left.values())])
summary = pd.DataFrame(rows, columns=["step", "what is removed", "lost_Natural", "lost_Synthetic", "lost_Total",
                                      "remaining_Natural", "remaining_Synthetic", "remaining_Total"])
summary[summary.columns[2:]] = summary[summary.columns[2:]].astype("Int64")
summary.to_csv(DELIV / "tables" / "redundancy_removal_summary.csv", index=False)

# the final counts have to agree with what was written
assert left == meta["ORIGIN"].value_counts().reindex(ORIGINS).to_dict()
assert sum(left.values()) == len(final) == sum(1 for l in open(DELIV / "fasta" / "amp_nonredundant.fasta") if l[0] == ">")
assert len(rem) + len(final) == sum(counts[o]["raw rows"] for o in ORIGINS)

print(summary.to_string(index=False))
print(f"\nCD-HIT: sent {len(to_cluster)} (>= {MIN_LEN} aa), read {n_read}, kept {len(cd_kept)}, lost {len(cd_lost)}; "
      f"{len(too_short)} peptides < {MIN_LEN} aa bypassed CD-HIT")
cd = rem[rem["STEP"].str.startswith("6")]
cross = pd.crosstab(cd["ORIGIN"], cd["COVERED_BY"].map(origin_of).rename("kept representative is"))
print("\nCD-HIT losses by origin of the lost peptide (rows) and of the representative that covers it (columns):")
print(cross.to_string())
print(f"final: {len(final)}  ({meta['ORIGIN'].value_counts().to_dict()})")
