"""Deliverable step 6 - the loss table as a spreadsheet with a plain-language note for every step (run 01 first).

Reads   tables/redundancy_removal_summary.csv, tables/x_to_alanine_counts.csv, tables/peptide_metadata.csv, the raw zip
Writes  sequences_lost.xlsx
"""
import io
import zipfile
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

DELIV = Path(__file__).resolve().parent.parent
summary = pd.read_csv(DELIV / "tables" / "redundancy_removal_summary.csv", dtype={"step": str})
xala = pd.read_csv(DELIV / "tables" / "x_to_alanine_counts.csv", index_col=0)
meta = pd.read_csv(DELIV / "tables" / "peptide_metadata.csv")

# peptides of the delivered files that carried an X in the raw data
zf = zipfile.ZipFile(DELIV.parent / "raw_data" / "dbaasp_datasets_corrected.zip")
raw = pd.concat([pd.read_csv(io.BytesIO(zf.read(f"by_origin/{o}.csv"))) for o in ("Natural", "Synthetic")])
raw_seq = dict(zip("DBAASP_" + raw["ID"].astype(str), raw["SEQUENCE"].astype(str).str.upper()))
had_x = meta[meta["ID"].map(lambda i: "X" in raw_seq[i])]
had_x_by = had_x["ORIGIN"].value_counts().reindex(["Natural", "Synthetic"], fill_value=0)

s = summary.set_index("step")
fin = s.loc["9"]
n_short = int((meta["LENGTH"] < 10).sum())
entering = int(s.loc["8", "remaining_Total"])
lost9 = int(fin["lost_Total"])

NOTES = {
    "1": ["Step 1 - exact duplicate row: the same entry appears twice in the export. The entry itself is kept,",
          "so the same ID can appear here and in the final set. With the complete export of November 2026 there are none",
          "(the earlier export had one page downloaded twice, which produced 2,000 such rows)."],
    "2": ["Step 2 - multi-chain entry: the sequence cell holds more than one chain. Its length and",
          "amino acid composition are those of several chains together, not of one peptide, so it",
          "cannot be compared with the single-chain peptides or clustered against them."],
    "3": ["Step 3 - sequence of only X: X stands for an amino acid that was not identified. A",
          "sequence made only of X records no amino acid at all, so there is nothing to measure."],
    "4": ["Step 4 - fewer than 3 real amino acids: a peptide needs at least 3 amino acids that are",
          "not X. The count is made before X is changed to Alanine, so X does not count."],
    "5": ["Step 5 - longer than 100 aa: sequences of more than 100 amino acids are proteins rather",
          "than peptides and are removed. The delivered files hold peptides of 3 to 100 aa."],
    "6": ["Step 6 - X is 30 % or more of the peptide: such a peptide is removed, because too much of",
          "it would be invented once X is changed to Alanine. Where X is below 30 %, every X is",
          f"changed to Alanine (this is not a loss): {xala.loc['peptides with X -> Ala', 'Natural']:,} Natural and "
          f"{xala.loc['peptides with X -> Ala', 'Synthetic']:,} Synthetic peptides,",
          f"{int(xala.loc['X residues -> Ala'].sum()):,} residues in all at this step; {len(had_x):,} of these peptides "
          f"({had_x_by['Natural']:,} Natural, {had_x_by['Synthetic']:,} Synthetic) are in the delivered files."],
    "7": ["Step 7 - duplicate sequence (same origin): searched after X -> Alanine, because the change",
          "can make two peptides identical. Upper and lower case are treated as the same; the first",
          "entry is kept."],
    "8": ["Step 8 - identical sequence in Natural and in Synthetic: the Natural entry is kept."],
    "9": ["Step 9 - CD-HIT at 90 % identity removes near-duplicate sequences, keeping one of each",
          "group (the longest). Natural and Synthetic are clustered together. Peptides shorter than",
          f"10 aa are not clustered: of the {entering:,} entering this step, {n_short:,} skip CD-HIT, {entering - n_short:,} go in,",
          f"{lost9:,} are removed and {entering - n_short - lost9:,} are kept; {entering - n_short - lost9:,} + {n_short:,} = {len(meta):,}."],
}
LABELS = {"0": "raw rows (DBAASP export)", "1": "exact duplicate row", "2": "multi-chain entry", "3": "sequence of only X",
          "4": "fewer than 3 real amino acids (not X)", "5": "longer than 100 aa", "6": "X is 30 % or more of the peptide",
          "7": "duplicate sequence (same origin, after X -> Ala)", "8": "identical sequence in Natural and Synthetic",
          "9": "CD-HIT >= 90 % identity"}

wb = Workbook()
ws = wb.active
ws.title = "Sequences lost"
ws.append(["Step", "What is removed", "Natural lost", "Synthetic lost", "Total lost", "Natural left", "Synthetic left", "Total left"])
for step, r in s.iterrows():
    vals = [None if pd.isna(r[c]) else int(r[c]) for c in ["lost_Natural", "lost_Synthetic", "lost_Total",
                                                           "remaining_Natural", "remaining_Synthetic", "remaining_Total"]]
    ws.append([step, LABELS[step], *vals])
ws.append([None, "FINAL (delivered FASTA files)", None, None, int(s["lost_Total"].sum()),
           int(fin["remaining_Natural"]), int(fin["remaining_Synthetic"]), int(fin["remaining_Total"])])
for c in ws[1]:
    c.font = Font(bold=True)
for c in ws[ws.max_row]:
    c.font = Font(bold=True)
for step in [str(i) for i in range(1, 10)]:
    ws.append([])
    for line in NOTES[step]:
        ws.append([line])
ws.column_dimensions["A"].width = 7
ws.column_dimensions["B"].width = 52
for col in "CDEFGH":
    ws.column_dimensions[col].width = 14
for row in ws.iter_rows(min_row=2, max_row=s.shape[0] + 2):
    for c in row[2:]:
        c.number_format = "#,##0"
        c.alignment = Alignment(horizontal="right")
ws.freeze_panes = "A2"
wb.save(DELIV / "sequences_lost.xlsx")
print("written", DELIV / "sequences_lost.xlsx", f"({len(had_x):,} delivered peptides carry a former X)")
