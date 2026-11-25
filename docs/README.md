# Documents

| File | What it is |
|---|---|
| [`research_guide/research_guide.pdf`](research_guide/research_guide.pdf) | In-depth guide to the whole study (24 pages): data collection, cleaning, CD-HIT, length and enrichment methods, Natural vs Synthetic, the class comparisons (source organism, synthesis type, function, target object), decision log, limitations, questions to prepare for, file guide. |
| [`article/article.pdf`](article/article.pdf) | Draft research article (11 pages) with the figures and findings. Red text marks placeholders (authors, affiliation, PepBench citation, repository URL). |

Both are built from the `.tex` files next to them; the figures are read from `deliverables/figures/` and `class_comparisons/`:

```bash
cd docs/research_guide && pdflatex research_guide.tex && pdflatex research_guide.tex
cd docs/article && pdflatex article.tex && pdflatex article.tex
```

Version 2 (November 2026) uses the complete DBAASP export (the page `peptides 16000-18000.csv` had been a copy of `peptides 18000-20000.csv`
in version 1) and adds the class comparisons. All conclusions of version 1 hold; numbers were updated (11,716 instead of 10,824 peptides).
