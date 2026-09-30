"""Step 1 - clean every CSV in the by_* folders.

The CSV files are overwritten in place (same names).
The untouched originals are kept in raw_data/dbaasp_datasets_corrected.zip.
A count of what was removed at every step goes to reports/cleaning_report.csv.
"""
import glob
import pandas as pd

report = []

for path in sorted(glob.glob("by_*/*.csv")):
    df = pd.read_csv(path)
    row = {"file": path, "rows_raw": len(df)}

    # tidy up spaces around the sequence
    df["SEQUENCE"] = df["SEQUENCE"].astype(str).str.strip()

    # basic checks (only reported, nothing removed here)
    row["missing_sequence"] = int((df["SEQUENCE"].isin(["", "nan"])).sum())
    row["length_mismatch_monomers"] = int(
        ((df["SEQUENCE"].str.len() != df["Length_of_Sequence"]) & ~df["SEQUENCE"].str.contains(" ")).sum()
    )

    # 1. exact duplicate rows
    df = df.drop_duplicates()
    row["after_exact_duplicates"] = len(df)

    # 2. multi-chain entries: several chains separated by spaces.
    #    Their Length_of_Sequence also counts the spaces, so it is wrong.
    df = df[~df["SEQUENCE"].str.contains(r"\s")]
    row["after_multichain_removed"] = len(df)

    # 3. sequences made only of X (no real amino acid at all)
    df = df[df["SEQUENCE"].str.upper().str.strip("X") != ""]
    row["after_all_X_removed"] = len(df)

    # 4. duplicate sequences (D-amino acids are lowercase, so compare in upper case)
    df = df[~df["SEQUENCE"].str.upper().duplicated()]
    row["rows_clean"] = len(df)

    row["ids_unique"] = bool(df["ID"].is_unique)
    df.to_csv(path, index=False)
    report.append(row)

report = pd.DataFrame(report)
report.to_csv("reports/cleaning_report.csv", index=False)
print(report[["file", "rows_raw", "rows_clean"]].to_string(index=False))
print("\nTotals:", report[["rows_raw", "after_exact_duplicates", "after_multichain_removed",
                           "after_all_X_removed", "rows_clean"]].sum().to_dict())
