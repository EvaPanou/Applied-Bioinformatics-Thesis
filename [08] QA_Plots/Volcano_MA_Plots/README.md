# Volcano and MA Plots

This directory contains visualizations used to inspect **differential gene expression results** derived from transcriptomic analysis. Two complementary plot types are included:

- **Volcano plots**
- **MA plots**

Both plot types visualize the **same differential expression results**, but emphasize different statistical properties of the data.

---

# Plot Examples

## Volcano Plot Examples

<div align="center">

<img src="Volcano_Plots/Volcano_A_FULL_meta..CONDITION..SLE.meta..TIME..24.31.weeks...meta..CONDITION..Healthy.meta..TIME..24.31.weeks.png" width="700">

**Figure 1.** Example volcano plot showing the relationship between log₂ fold change and statistical significance.

</div>

<br>

<div align="center">

<img src="Volcano_Plots/Volcano_B_FULL_meta..CONDITION..SLE.meta..TIME..32.40.weeks...meta..CONDITION..SLE.meta..TIME..24.31.weeks.png" width="700">

**Figure 2.** Volcano plot illustrating differential expression patterns across another comparison.

</div>

Volcano plots display the relationship between **effect size** and **statistical significance**.

Each point corresponds to a gene and is positioned according to:

| Axis | Description |
|-----|-------------|
| X-axis | log₂ Fold Change |
| Y-axis | −log₁₀(FDR or adjusted p-value) |

Dashed guide lines indicate the thresholds used to highlight potentially relevant genes:

- **|log₂FC| ≥ 1**
- **FDR < 0.05**

Genes located toward the **upper left or upper right regions** of the plot typically represent transcripts with both:

- large expression changes  
- strong statistical significance

These genes often represent candidate **differentially expressed genes**.

---

## MA Plot Examples

<div align="center">

<img src="MA_Plots/MA_A_FULL_meta..CONDITION..SLE.meta..TIME..16.23.weeks...meta..CONDITION..Healthy.meta..TIME..16.23.weeks.png" width="700">

**Figure 3.** MA plot showing the relationship between average expression and log₂ fold change.

</div>

<br>

<div align="center">

<img src="MA_Plots/MA_B_FULL_meta..CONDITION..SLE.meta..TIME..32.40.weeks...meta..CONDITION..SLE.meta..TIME..24.31.weeks.png" width="700">

**Figure 4.** MA plot illustrating the distribution of fold changes across expression levels.

</div>

MA plots visualize the relationship between **average gene expression** and **log₂ fold change**.

| Axis | Description |
|-----|-------------|
| X-axis | Average expression (log₂ scale) |
| Y-axis | log₂ Fold Change |

The name **MA plot** derives from:
 - M = log ratio (log₂ fold change)
 - A = average expression


Horizontal dashed lines at **log₂FC = ±1** indicate approximate **two-fold expression change boundaries**.

Most genes form a dense band around **log₂FC ≈ 0**, indicating that the majority of genes show minimal expression differences between the compared groups.

---

# Relationship Between the Two Plot Types

Although MA plots and volcano plots display different axes, they represent the **same underlying differential expression results**.

| Plot | X-axis | Y-axis | Purpose |
|----|----|----|----|
| MA Plot | Average expression | log₂ Fold Change | Inspect fold change across expression levels |
| Volcano Plot | log₂ Fold Change | −log₁₀(FDR) | Highlight statistically significant genes |

Together, these visualizations provide **complementary perspectives on gene expression patterns**.

---

# Summary

The plots in this directory provide a **visual overview of gene expression differences across the analyzed comparisons**. MA plots emphasize the distribution of fold changes across expression intensities, while volcano plots highlight genes that combine strong effect sizes with statistical significance.

These visualizations serve as **exploratory tools for interpreting differential expression results** and identifying genes that may warrant further investigation.

---

# Notes

The dashed boundaries shown in the plots correspond to the thresholds used in the plotting script:

- **log₂ fold change cutoff:** |log₂FC| ≥ 1  
- **significance cutoff:** FDR < 0.05

These thresholds are used to visually highlight genes that exhibit both substantial expression change and statistical support in the differential expression results.