# README — Visualization_PCA_Volcano_MA.R

### important
*"In the PCA based on the LIMMA-derived DEG subset, SLE and Healthy samples separate primarily along PC2 (22.5% of total variance), indicating that disease status is captured in the second principal component. PC1 (36.4%) mainly reflects other biological or technical sources of variation, such as interindividual heterogeneity or gestational stage. In contrast, when using the larger consensus DEG set (500 genes), disease status becomes the dominant source of variance (PC1 = 37.9%), resulting in stronger visual separation between groups."*

This document accompanies `Visualization_PCA_Volcano_MA.R` and provides **academic-style** context and practical instructions for generating **diagnostic visualizations** from the SLE pregnancy transcriptomic dataset. The script produces (i) a **volcano plot** and **MA plot** from a LIMMA `FULL` table, and (ii) two **PCAs**: one using the **entire expression matrix** and one restricted to **DEG genes**.

---

## 1) Scientific context (short)

The figures produced by this script serve as **global quality and interpretability checks** that complement differential-expression results. The **volcano** and **MA** plots visualize the **effect-size vs. significance** and **expression-intensity vs. fold-change** landscapes for a selected contrast, respectively. The **PCA analyses** (all genes vs. DEG subset) assess whether the **dominant axes of variance** segregate samples by **Condition** and/or **Gestational Stage**, and whether restricting to DEGs enhances biological separation.

---

## 2) Inputs & configuration

Open the script and adjust the **CONFIG** block at the top:

```r
EXPR_FILE <- "C:/.../08_df_no_low_genes.tsv"          # rows = genes, cols = samples (log2)
META_FILE <- "C:/.../GSE108497_updated_metadata_II.csv" # must contain Sample, Condition, time_point
OUTDIR    <- "C:/.../Results"

# (optional) Pick a specific LIMMA FULL table; otherwise first FULL TSV is auto-detected under Results/
CONTRAST_TSV <- NULL  # e.g. file.path(OUTDIR, "A_SLE_vs_Healthy", "A_FULL_SLE.vs.Healthy.at.<16.weeks.tsv")

# (optional) DEG list used for the DEG-restricted PCA
DEG_LIST_TSV <- file.path(OUTDIR, "ML", "ML_Features_AllDEGs.tsv")
```

**Requirements (R packages):** `tidyverse`, `data.table`, `ggplot2`, `ggrepel`.

**Assumptions:**  
- Expression matrix is **log2**-scaled and samples align to metadata via `Sample`.  
- LIMMA `FULL` table contains `logFC`, `adj.P.Val`, and (for MA) `AveExpr` columns.  
- DEG list TSV has a column named `Gene` (or `gene`, `SYMBOL`, `symbol`).

---

## 3) How to run

From R:

```r
source("Visualization_PCA_Volcano_MA.R")
```

The script will:  
1. Load inputs and **auto-select** a `FULL` table if `CONTRAST_TSV` is not specified.  
2. Generate **Volcano.png** and **MA.png** under `Results/QA_Plots/`.  
3. Run **PCA** on **All Genes** and on the **DEG subset**, writing `PCA_AllGenes.png` and `PCA_DEG.png` to `Results/QA_Plots/`.

---

## 4) Outputs (with figure legends)

| File | Description | Notes |
|---|---|---|---|
| `QA_Plots/Volcano.png` | Volcano plot (log2FC vs. −log10(FDR)) | Dashed lines indicate |log2FC| cut and FDR=0.05; labeled top genes by FDR. |
| `QA_Plots/MA.png` | MA plot (AveExpr vs. log2FC) | Dashed lines show |log2FC| cut; helps detect intensity-dependent trends. |
| `QA_Plots/PCA_AllGenes.png` | PCA on all genes | PC1/PC2 colored by **Condition**, shaped by **time_point**; evaluates dominant variance structure. |
| `QA_Plots/PCA_DEG.png` | PCA on DEG genes | Tests whether DEGs enhance SLE vs. Healthy separation and stage structure. |

**Interpretation tips:**  
- **Volcano:** look for distinct clouds of significant genes; symmetry suggests balanced up/down calls.  
- **MA:** check for bias at low intensity; strong trends may indicate technical effects.  
- **PCA (All vs. DEG):** clearer grouping in the **DEG PCA** supports that your DEG selection captures disease/stage signal.

---

## 5) Parameters & customization

- Change FDR and |log2FC| thresholds within `plot_volcano()` and `plot_MA()` calls if needed.  
- Set `CONTRAST_TSV` to a specific `A_FULL_...tsv` to target a particular timepoint/contrast.  
- Increase `top_n_labels` in `plot_volcano()` to annotate more genes.

---

## 6) Troubleshooting

- **“No FULL topTable TSVs found”**: ensure you have run the LIMMA pipeline and `Results/` contains `A_FULL_*.tsv` or `B_FULL_*.tsv`.  
- **MA requires `AveExpr`**: confirm your `FULL` table includes the `AveExpr` column (produced by `topTable(..., number=Inf)`).  
- **PCA errors (too few genes)**: verify the DEG list exists or allow fallback to significant genes derived from the selected FULL table.

---

## 7) Citation note

When reporting these diagnostics, name the figures explicitly (e.g., *“Volcano plot for SLE vs. Healthy at <timepoint>”*) and reference the limma framework appropriately in your Methods chapter.
