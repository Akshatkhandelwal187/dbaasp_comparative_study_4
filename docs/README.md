# Documents

| File | What it is |
|---|---|
| [`research_guide/research_guide.pdf`](research_guide/research_guide.pdf) | In-depth guide to the whole study (26 pages): data collection, cleaning, CD-HIT, length and enrichment methods, survey of classes, decision log, limitations, questions to prepare for, file guide. |
| [`article/article.pdf`](article/article.pdf) | Draft research article (11 pages) with the figures and findings. Red text marks placeholders (authors, affiliation, PepBench citation, repository URL) and the re-run needed before submission. |

Both are built from the `.tex` files next to them; the figures are read from `deliverables/figures/`, `survey/figures/` and `analysis/`:

```bash
cd docs/research_guide && pdflatex research_guide.tex && pdflatex research_guide.tex
cd docs/article && pdflatex article.tex && pdflatex article.tex
```

**Data issue found while writing these:** `data collection/combined_df/peptides 16000-18000.csv` and `peptides 18000-20000.csv` are identical
(both hold IDs 18253-20252). The export page for DBAASP IDs of about 16253-18252 (roughly 2,000 entries, at least 87 of them natural peptides with a kingdom)
is missing from the study. Re-export that page and re-run the pipeline to fix it (guide, section 3.3).
