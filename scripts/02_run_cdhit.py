"""Step 2 - remove redundancy with CD-HIT (>= 90 % identity) for EVERY cleaned dataset.

For each CSV in by_*/ this makes a folder cdhit/<dimension>/<group>/ containing
    input.fasta              all cleaned sequences (upper case, header = ID)
    input_ge10aa.fasta       the sequences >= 10 aa (these go into cd-hit)
    cdhit_ge10aa.fasta       cd-hit output (one representative per cluster)
    cdhit_ge10aa.fasta.clstr cd-hit cluster file
    nr90.fasta               final non-redundant set (cd-hit output + peptides < 10 aa)
and the matching non-redundant table nonredundant_csv/<dimension>/<group>.csv

Why peptides < 10 aa skip cd-hit:
  * cd-hit does not handle very short sequences reliably (its own minimum is 10 aa)
  * two different peptides of <= 9 aa can never be >= 90 % identical
    (one mismatch in 9 residues = 88.9 %), and identical ones were already removed in step 1.

cd-hit settings: -c 0.9 (90 % identity)  -n 5 (word size required for 0.9)  -d 0 (keep full ID in output)
"""
import glob
import os
import subprocess
import pandas as pd

MIN_LEN = 10
summary = []

for path in sorted(glob.glob("by_*/*.csv")):
    dimension = path.split("/")[0]
    group = os.path.basename(path)[:-4]
    folder = f"cdhit/{dimension}/{group}"
    os.makedirs(folder, exist_ok=True)
    os.makedirs(f"nonredundant_csv/{dimension}", exist_ok=True)

    df = pd.read_csv(path)
    df["seq"] = df["SEQUENCE"].str.upper()      # D-amino acids (lowercase) -> upper case

    def write_fasta(table, filename):
        with open(filename, "w") as f:
            for id_, seq in zip(table["ID"], table["seq"]):
                f.write(f">{id_}\n{seq}\n")

    long_seqs = df[df["seq"].str.len() >= MIN_LEN]
    short_seqs = df[df["seq"].str.len() < MIN_LEN]
    write_fasta(df, f"{folder}/input.fasta")
    write_fasta(long_seqs, f"{folder}/input_ge10aa.fasta")

    # run cd-hit on the >= 10 aa sequences
    kept_ids = set(short_seqs["ID"])
    if len(long_seqs) > 0:
        result = subprocess.run(
            ["cd-hit", "-i", f"{folder}/input_ge10aa.fasta", "-o", f"{folder}/cdhit_ge10aa.fasta",
             "-c", "0.9", "-n", "5", "-d", "0", "-M", "0", "-T", "0"],
            capture_output=True, text=True, check=True)
        open(f"{folder}/cdhit.log", "w").write(result.stdout)
        with open(f"{folder}/cdhit_ge10aa.fasta") as f:
            kept_ids |= {line[1:].strip() for line in f if line.startswith(">")}

    # final non-redundant table + fasta (same order as the cleaned csv)
    nr = df[df["ID"].astype(str).isin({str(i) for i in kept_ids})]
    nr.drop(columns="seq").to_csv(f"nonredundant_csv/{dimension}/{group}.csv", index=False)
    write_fasta(nr, f"{folder}/nr90.fasta")

    summary.append({"dimension": dimension, "group": group, "rows_clean": len(df),
                    "shorter_than_10aa_not_clustered": len(short_seqs),
                    "sent_to_cdhit": len(long_seqs),
                    "kept_from_cdhit": len(nr) - len(short_seqs),
                    "rows_nr90": len(nr),
                    "removed_by_cdhit": len(df) - len(nr)})

summary = pd.DataFrame(summary)
summary.to_csv("reports/cdhit_summary.csv", index=False)
print(summary.to_string(index=False))

# refresh manifest.csv (same file name) with raw / clean / non-redundant row counts
old = pd.read_csv("manifest.csv")
label = dict(zip(old["file"].apply(os.path.basename), old["value"]))
raw = pd.read_csv("reports/cleaning_report.csv")
raw["dimension"] = raw["file"].str.split("/").str[0]
raw["group"] = raw["file"].apply(lambda p: os.path.basename(p)[:-4])
m = summary.merge(raw[["dimension", "group", "rows_raw"]], on=["dimension", "group"])
m["value"] = [label.get(g + ".csv", g) for g in m["group"]]
m["file"] = m["dimension"] + "/" + m["group"] + ".csv"
m = m.rename(columns={"rows_raw": "n_rows_raw", "rows_clean": "n_rows_clean", "rows_nr90": "n_rows_nr90"})
m[["dimension", "value", "n_rows_raw", "n_rows_clean", "n_rows_nr90", "file"]].to_csv("manifest.csv", index=False)
