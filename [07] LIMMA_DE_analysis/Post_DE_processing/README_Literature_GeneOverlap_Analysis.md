# README — Literature_GeneOverlap_Analysis.R

This document accompanies `Literature_GeneOverlap_Analysis.R` and provides an **academic-style** rationale and instructions for **validating DEG-derived feature sets** against **prior biological knowledge** from the literature. The script quantifies overlaps between your **All DEGs** and **DEG∩HighVariance** gene sets and a **user-supplied literature gene list**, summarizing results with **Jaccard index**, optional **Fisher’s exact enrichment**, and an **overlap bar figure**.

---

## 1) Scientific context (short)

Comparing data-driven DEGs to **established gene signatures** (e.g., interferon-stimulated genes in SLE) is a standard practice to **corroborate biological validity**. A sizable and statistically enriched overlap indicates that the feature sets are **consistent with known mechanisms**, strengthening the interpretability and credibility of subsequent analyses.

---

## 2) Inputs & configuration

Edit the **CONFIG** block at the top of the script:

```r
RESULTS_DIR <- "C:/.../Results"
ALLDEG_TSV  <- file.path(RESULTS_DIR, "ML", "ML_Features_AllDEGs.tsv")
DEG_HV_TSV  <- file.path(RESULTS_DIR, "ML", "ML_Features_DEG_HighVariance.tsv")

# Provide a literature list:
LITERATURE_TSV <- NULL            # (optional) path to TSV/CSV; must have a 'Gene' column
LITERATURE_GENES <- c("IFI44L","IFIT3","MX1","OAS1","ISG15","STAT1","IRF7","RSAD2","IFIT1","OASL")  # editable
```

**Requirements (R packages):** `tidyverse`, `data.table`, `ggplot2`.

**Assumptions:**  
- The gene list TSVs contain a column named `Gene` (or `gene`, `SYMBOL`, `symbol`).  
- If you know the **array universe size** (≈20,000 for many human arrays), set `UNIVERSE_SIZE` to enable **Fisher enrichment testing**.

---

## 3) How to run

From R:

```r
source("Literature_GeneOverlap_Analysis.R")
```

Outputs will be written under: `Results/Literature_Overlap/`

---

## 4) Outputs (with descriptions)

| File | Description |
|---|---|
| `Overlap_Summary.tsv` | Table summarizing set sizes, intersection, **Union**, **Jaccard index**, and optional **Fisher’s exact p-value** (if `UNIVERSE_SIZE` is set). |
| `Overlap_AllDEGs__Literature.txt` | Gene names in the intersection between **AllDEGs** and **Literature**. |
| `Overlap_DEG_HV__Literature.txt` | Gene names in the intersection between **DEG∩HighVariance** and **Literature**. |
| `Overlap_Bar.png` | Bar chart showing the counts of overlapping genes for each comparison. |

---

## 5) Interpretation guidance

- **Jaccard index** close to 0.1–0.3 can already be meaningful for large gene spaces; interpret relative to list sizes.  
- **Fisher’s exact p-value** (if computed) tests whether the overlap is larger than expected by chance.  
- Larger overlaps for **DEG∩HighVariance** versus **AllDEGs** suggest that adding a **variance filter** concentrated biologically salient signals.

---

## 6) Customization tips

- Replace `LITERATURE_GENES` with a curated SLE signature from your review (e.g., interferon score genes).  
- Point `LITERATURE_TSV` to a file with a `Gene` column to avoid editing the script.  
- If you have **multiple literature gene sets**, run the script repeatedly or extend it to loop across sets and produce a multi-bar figure.

---

## 7) Reporting & citation

When reporting, state the **size of each set**, the **intersection**, **Jaccard**, and (if used) **Fisher’s p-value**, and cite the **source** of the literature signature (paper or database). This strengthens the argument that your features are **mechanistically grounded** and not artifacts of statistical selection.
