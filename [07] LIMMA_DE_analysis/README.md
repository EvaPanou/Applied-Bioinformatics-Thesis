# Differential Gene Expression and Feature Selection from Temporal and Conditional Contrasts in SLE Pregnancy

## 1. Repository Organization
```
├── DE_Analysis/
├── Post_DE_processing/ -> THIS FOLDER IS STIL UNDER PROCESSING, AND VIRTUALLY WITH NO IMPORTANT INFORMATION AT THE MOMENT
└── README.md
```
---

## 2. Overview
This repository contains the processed results, intermediate outputs, and statistical summaries from the **Differential Expression (DEG) analysis** performed on the **GSE108497** dataset. The study explores **Systemic Lupus Erythematosus (SLE) in pregnancy**, using both **pairwise (conditional)** and **temporal (longitudinal)** comparisons to identify genes that are differentially expressed between SLE and healthy pregnancies and over time.

This pipeline provides a complete framework for analyzing **time-resolved and condition-dependent transcriptomic data**. The approach identifies genes that are:
- Differentially expressed across both **disease** and **time** dimensions.
- Statistically validated through **limma’s empirical Bayes moderation**.
- Biologically interpretable via clustering and enrichment.
- Ready for use in **machine learning-based classification**.

By combining classical DEG analysis with multivariate variance filtering, this workflow ensures both biological interpretability and predictive relevance.

All analyses were performed in **R** using the **limma** package, integrating both condition- and time-dependent models to extract statistically and biologically meaningful DEGs. Subsequent steps included gene clustering, enrichment analysis, and feature selection for machine learning (ML) applications.

---

## 3. Methodological Summary

### 3.1 Experimental Design
- **Dataset:** GSE108497 (Illumina HumanHT-12 v4 microarray, log2-normalized expression).
- **Factors:**
  - *Condition*: SLE vs. Healthy.
  - *Time*: Multiple pregnancy stages (longitudinal sampling).
  - *Interaction*: Condition × Time.

### 3.2 Statistical Model
The **limma** linear modeling framework was used to fit an expression model that included both main and interaction effects:
- **Main effects:** Condition and Time.
- **Interaction effect:** Condition × Time.

Two types of analyses were performed:
1. **Pairwise (conditional)** contrasts: SLE vs. Healthy at each timepoint.
2. **Temporal (longitudinal)** contrasts: differences between timepoints within each condition.

**Significance thresholds:**
- Adjusted p-value (FDR) < 0.05
- |log2 Fold-Change| > 1

### 3.3 Unification of DEGs
DEG lists from all contrasts were combined and deduplicated to produce a unified gene list:
- FDR thresholds tested: 0.05 and 0.01.
- Only unique genes were kept.
- Downstream analyses (heatmaps, enrichment, ML features) used **FDR < 0.05** and **|logFC| > 1**.
- ML feature gene lists/matrices were extracted

### Further Analysis: 
- Enrichement analysis was performed, and various graphs were produced, such as heatmaps and upset plots. 


