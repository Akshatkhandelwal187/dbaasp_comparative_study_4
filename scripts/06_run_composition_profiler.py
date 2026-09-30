"""Step 6 - amino-acid enrichment with the official Composition Profiler (github.com/vvacic/cprofiler)
for the Natural and the Synthetic dataset, each against the same SwissProt background.

Setup (once):
    git clone https://github.com/vvacic/cprofiler
    pip install -e cprofiler        # -e is needed: the normal pip install does not ship the data folder

This script does exactly what
    cprof discover -Q <query.fasta> -B <background.fasta> -I 10000 -A 0.05 -b
does (same functions from the official package), but saves the table as a CSV
(the command line tool only prints it, cut off, on the screen).
The picture is made with the official command
    cprof plot -Q <query.fasta> -B <background.fasta> -X alpha -O <png>
"""
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cprofiler.aminoacid import AminoAcid
from cprofiler.fasta import Fasta
from cprofiler.profile import CompositionProfiler

BACKGROUND = "composition_profiler/background/swissprot51_5k.fasta"
OUT = "composition_profiler/results"
ITERATIONS = 10000                 # default of the tool
ALPHA = 0.05                       # default of the tool

np.random.seed(1)                  # so that the (random) p-values can be repeated

alphabet = AminoAcid.get_order("alpha")            # A C D E F G H I K L M N P Q R S T V W Y
groups = AminoAcid.get_groups()
group_names = AminoAcid.get_group_names()
alpha = ALPHA / (len(alphabet) + len(groups))      # Bonferroni correction (the -b option)

background = Fasta.count_chars(Fasta.read(BACKGROUND), alphabet)

results = {}
for name in ["Natural", "Synthetic"]:
    query_file = f"fasta_files/by_origin/{name}/{name}.fasta"
    query = Fasta.count_chars(Fasta.read(query_file), alphabet)

    table = CompositionProfiler.discover(query, background, alphabet=alphabet, groups=groups,
                                         group_names=group_names, iterations=ITERATIONS,
                                         alpha_value=alpha)
    table.to_csv(f"{OUT}/{name}_vs_swissprot.csv", index=False)
    results[name] = table

    # bar plot from the official tool
    subprocess.run(["cprof", "plot", "-Q", query_file, "-B", BACKGROUND, "-X", "alpha",
                    "-Y", f"({name} - SwissProt) / SwissProt", "-O", f"{OUT}/{name}_vs_swissprot.png"],
                   check=True)

# Natural and Synthetic side by side
side = results["Natural"].merge(results["Synthetic"], on="test_name", suffixes=("_Natural", "_Synthetic"))
side.to_csv(f"{OUT}/natural_vs_synthetic_side_by_side.csv", index=False)

# one simple figure: the 20 amino acids, Natural vs Synthetic
aa = side.iloc[:20]
x = np.arange(20)
plt.figure(figsize=(10, 4))
plt.bar(x - 0.2, aa["effect_Natural"], width=0.4, label="Natural")
plt.bar(x + 0.2, aa["effect_Synthetic"], width=0.4, label="Synthetic")
plt.axhline(0, color="black", linewidth=0.8)
plt.xticks(x, aa["test_name"])
plt.ylabel("(peptides - SwissProt) / SwissProt")
plt.title("Amino acid enrichment (+) / depletion (-) vs SwissProt")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT}/natural_vs_synthetic_amino_acids.png", dpi=150)

print(side.round(3).to_string(index=False))
