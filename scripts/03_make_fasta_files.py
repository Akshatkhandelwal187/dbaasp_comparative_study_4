"""Step 3 - collect the non-redundant FASTA files of the two datasets we study
(by_origin and by_length_class) into fasta_files/<dimension>/<subgroup>/<subgroup>.fasta.
These are the input files for Composition Profiler.
"""
import os
import shutil

for dimension in ["by_origin", "by_length_class"]:
    for group in sorted(os.listdir(f"cdhit/{dimension}")):
        target = f"fasta_files/{dimension}/{group}"
        os.makedirs(target, exist_ok=True)
        shutil.copy(f"cdhit/{dimension}/{group}/nr90.fasta", f"{target}/{group}.fasta")
        print(f"{target}/{group}.fasta")
