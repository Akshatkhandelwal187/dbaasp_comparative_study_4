# DBAASP comparative study of antimicrobial peptide classes

This repository builds a cleaned, redundancy-reduced set of antimicrobial peptides (AMPs) from DBAASP and compares classes of peptides
by **length** and **amino acid composition** (Composition Profiler enrichment against SwissProt proteins and non-antimicrobial peptides):

* **Natural vs Synthetic**: [`deliverables/`](deliverables/README.md)
* **Source organism, synthesis type, function (target group) and target object**: [`class_comparisons/`](class_comparisons/README.md)

All labels are existing DBAASP columns; no new peptide properties are calculated. Scripts are short, numbered and read top to bottom.
A full explanation of the study (methods, decisions, results, limitations) is in [`docs/`](docs/README.md): an in-depth research guide and a draft article.

## Folder map

```
data collection/
    combined_df/peptides *.csv      DBAASP peptide list, 13 export pages (25,542 peptides)
    kingdom/peptides *.csv          DBAASP peptide lists per kingdom
    dbaasp_data_final.ipynb         the original notebook (annotation and subsets)
    build_subsets.py                script version of the notebook: rebuilds raw_data/ from the export pages
raw_data/dbaasp_datasets_corrected.zip   subsets by origin, length class, domain, kingdom, target group, target object (+ manifest)
manifest.csv                        rows per subset
deliverables/                       reference dataset: cleaning, CD-HIT 90 %, loss accounting, Natural vs Synthetic figures
class_comparisons/                  the same figures and tables for the other classifications
docs/                               research guide and draft article (LaTeX + PDF)
```

## Pipeline

```bash
pip install -r requirements.txt
sudo apt-get install cd-hit                          # tested with 4.8.1
git clone https://github.com/vvacic/cprofiler && pip install -e cprofiler      # -e is required

python "data collection/build_subsets.py"            # export pages -> raw_data/dbaasp_datasets_corrected.zip
python deliverables/scripts/01_remove_redundancy.py  # cleaning + CD-HIT -> 11,716 peptides (2,552 Natural, 9,164 Synthetic)
python deliverables/scripts/02_make_figures.py
python deliverables/scripts/03_check_cdhit_result.py
python deliverables/scripts/04_aa_enrichment.py
python deliverables/scripts/05_enrichment_figures.py
python deliverables/scripts/06_sequences_lost_xlsx.py
python class_comparisons/scripts/01_build_class_table.py
python class_comparisons/scripts/02_length_figures.py
python class_comparisons/scripts/03_aa_enrichment.py # about 30 minutes
python class_comparisons/scripts/04_enrichment_figures.py
```

## History

* **September 2026**: first pass (per-subset cleaning and CD-HIT, Natural vs Synthetic and length-class statistics). It had a CD-HIT setting
  that silently dropped all 10-aa peptides and was based on an incomplete export; it was **removed in November 2026** (see git history).
* **October-November 2026**: `deliverables/` (stricter cleaning, combined CD-HIT with the fix, two backgrounds).
* **November 2026**: the export page `peptides 16000-18000.csv` turned out to be a second copy of `peptides 18000-20000.csv`
  (2,000 DBAASP entries, IDs 16253-18252, were missing). It was downloaded again and everything was re-run; the class survey became `class_comparisons/`.
