# Quality Assessment Plots

This directory contains visualization outputs generated during the **differential expression analysis workflow** of the thesis project.

The plots were generated using the R script:

`PCA_and_plots_script.R`

This script loads the normalized expression matrix, metadata, and differential expression results from the **LIMMA DEG analysis**, and produces several diagnostic and exploratory plots including:

* Principal Component Analysis (PCA)
* Volcano plots
* MA plots

These visualizations are used to assess **global transcriptomic structure**, **differential expression behavior**, and **data quality**.

---

# Directory Structure

```
[08]_QA_Plots
│
├── Input/
│
├── PCA_Plots/
│   ├── PCA_AllGenes.png
│   ├── PCA_Old_DEG_Subset.png
│   ├── PCA_LIMMA_DEG_Subset.png
│
├── Volcano_MA_Plots/
│   ├── MA_*.png
│   ├── Volcano_*.png
│
└── PCA_and_plots_script.R
```

---

# Folder Description

## Input

This folder contains the **input files required for generating the plots**.

These files originate from the **LIMMA differential expression analysis** and include:

* normalized expression matrices
* sample metadata
* differential expression result tables
* DEG lists used for PCA subset analyses

These files are loaded by the visualization script and used to generate the PCA, volcano, and MA plots.

---

## PCA_Plots

This folder contains **Principal Component Analysis (PCA)** visualizations generated from different gene subsets.

The PCA plots are used to examine the **global structure of the transcriptomic dataset** and evaluate whether samples cluster according to:

* disease condition (SLE vs Healthy)
* pregnancy time point

The PCA analyses are the **main focus of this README**.

---

## Volcano_MA_Plots

This folder contains additional plots generated for the differential expression contrasts:

### Volcano Plots

These plots visualize the relationship between:

* **log₂ fold change**
* **statistical significance (-log₁₀ FDR)**

and highlight genes that pass differential expression thresholds.

### MA Plots

MA plots display:

* **average gene expression (A)**
* **log₂ fold change (M)**

and are commonly used to inspect the distribution of expression changes across the transcriptome.

Because multiple contrasts were generated during the LIMMA analysis, this folder contains **multiple volcano and MA plots**.

A dedicated README inside this folder provides a detailed explanation of these plots.

---

# Principal Component Analysis (PCA)

Principal Component Analysis (PCA) was performed to explore global transcriptomic variation across samples and evaluate whether gene expression patterns distinguish **Systemic Lupus Erythematosus (SLE)** patients from **healthy pregnant controls**.

The analysis was conducted using normalized **log₂-transformed expression values** after preprocessing and filtering of low-expression genes.

Samples are colored by **Condition (Healthy vs SLE)** and shaped by **pregnancy time point**:

* <16 weeks
* 16–23 weeks
* 24–31 weeks
* 32–40 weeks
* Postpartum (PP)

---

# Software and Packages

All analyses and visualizations were performed in **R** using the following packages:

* `tidyverse`
* `data.table`
* `ggplot2`
* `ggrepel`

Principal Component Analysis was computed using the base R function **`prcomp()`**, with centering and scaling applied to the expression matrix.

---

# 1. PCA Using All Genes

<p align="center">
<img src="PCA_Plots/PCA_AllGenes.png" width="650">
</p>

**Figure 1. PCA using all genes after preprocessing.**

This PCA represents the global transcriptional structure of the dataset using all genes retained after preprocessing.

### Observations

* Substantial overlap exists between **SLE** and **healthy** samples.
* Pregnancy timepoints are broadly distributed across the PCA space.
* This indicates that most genes are **not strongly condition-specific**.

This analysis primarily serves as a **quality control step** to assess overall sample structure and detect potential outliers.

---

# 2. PCA Using Condition-Based DEG Subset

<p align="center">
<img src="PCA_Plots/PCA_Old_DEG_Subset.png" width="650">
</p>

**Figure 2. PCA based on a previously defined DEG subset.**

This PCA uses a gene subset identified through **between-condition differential expression analysis (SLE vs Healthy)** without explicitly modeling the longitudinal structure of the dataset.

### Observations

* Separation between SLE and healthy samples becomes more apparent.
* The clustering primarily reflects **disease-associated transcriptional variation**.
* Temporal effects related to pregnancy stage are not explicitly accounted for.

---

# 3. PCA Using LIMMA Longitudinal DEG Subset

<p align="center">
<img src="PCA_Plots/PCA_LIMMA_DEG_Subset.png" width="650">
</p>

**Figure 3. PCA based on the final LIMMA-derived DEG subset.**

This PCA uses genes identified through a **LIMMA differential expression model** that incorporates both:

* **between-condition effects** (SLE vs Healthy)
* **longitudinal effects** across pregnancy timepoints.

### Observations

* Samples display clearer separation between disease conditions.
* The structure of the PCA reflects both **disease-related** and **time-dependent transcriptional variation**.

This gene set represents the **final DEG subset used for downstream analyses**.

---

# Summary

| PCA Analysis     | Gene Set                    | Model Type                     | Purpose                              |
| ---------------- | --------------------------- | ------------------------------ | ------------------------------------ |
| All Genes        | Full filtered transcriptome | None                           | Global dataset structure             |
| Old DEG subset   | Condition-based DEGs        | Between-condition analysis     | Highlight disease signal             |
| LIMMA DEG subset | Final DEG set               | Longitudinal + condition model | Final biologically relevant features |

Together, these analyses demonstrate how **feature selection and modeling strategy influence the ability of PCA to separate biological groups in transcriptomic data**.
