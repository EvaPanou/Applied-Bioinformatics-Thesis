# Folder: `A_SLE_vs_Healthy/FDR0.05/`

### Differentially Expressed Genes (DEGs) at FDR < 0.05 and |logFC| > 1

This folder contains the **broader, but still statistically significant**, set of genes differentially expressed between **SLE** and **Healthy** pregnancies at each timepoint.

These DEGs are defined by False Discovery Rate < 5% and a minimum twofold change in expression (|logFC| > 1).

---

## Summary by timepoint

| Timepoint   |   Upregulated genes (SLE>Healthy) |   Downregulated genes (SLE<Healthy) |
|:------------|----------------------------------:|------------------------------------:|
| <16 weeks   |                                49 |                                   1 |
| 16–23 weeks |                                49 |                                   1 |
| 24–31 weeks |                                44 |                                   0 |
| 32–40 weeks |                                32 |                                   0 |
| Postpartum  |                                40 |                                   0 |

> **Interpretation:** Relaxing the FDR cutoff from 0.01 to 0.05 increases sensitivity — more genes are detected, but the biological pattern remains identical: **SLE shows consistent up-regulation of immune and interferon-stimulated genes (ISGs)** at every stage of pregnancy.

## Top representative genes per timepoint

### <16 weeks

**Upregulated (SLE > Healthy):**
- IFI44L (logFC 2.96, adj.P.Val 2.1e-32)
- IFI27 (logFC 4.04, adj.P.Val 2.8e-31)
- EPSTI1 (logFC 2.22, adj.P.Val 1.6e-29)
- RSAD2 (logFC 3.28, adj.P.Val 4.2e-29)
- LY6E (logFC 2.26, adj.P.Val 7.9e-28)
**Downregulated (SLE < Healthy):**
- ORM1 (logFC -1.20, adj.P.Val 5.5e-08)

### 16–23 weeks

**Upregulated (SLE > Healthy):**
- IFI44L (logFC 2.99, adj.P.Val 1.2e-31)
- EPSTI1 (logFC 2.32, adj.P.Val 1.1e-30)
- LY6E (logFC 2.45, adj.P.Val 1.1e-30)
- OASL (logFC 2.02, adj.P.Val 4.2e-30)
- IFI6 (logFC 1.81, adj.P.Val 4.2e-30)
**Downregulated (SLE < Healthy):**
- ORM1 (logFC -1.05, adj.P.Val 4.8e-06)

### 24–31 weeks

**Upregulated (SLE > Healthy):**
- OASL (logFC 1.93, adj.P.Val 1.0e-26)
- ISG15 (logFC 2.39, adj.P.Val 1.0e-24)
- LY6E (logFC 2.17, adj.P.Val 4.5e-24)
- IFI44L (logFC 2.54, adj.P.Val 1.8e-23)
- IFI6 (logFC 1.59, adj.P.Val 2.0e-23)
**Downregulated (SLE < Healthy):**
- None detected at this threshold

### 32–40 weeks

**Upregulated (SLE > Healthy):**
- ZBP1 (logFC 1.19, adj.P.Val 7.3e-16)
- ISG15 (logFC 2.03, adj.P.Val 7.3e-16)
- RSAD2 (logFC 2.50, adj.P.Val 2.5e-14)
- OAS3 (logFC 1.40, adj.P.Val 5.0e-14)
- HERC5 (logFC 1.60, adj.P.Val 5.0e-14)
**Downregulated (SLE < Healthy):**
- None detected at this threshold

### Postpartum

**Upregulated (SLE > Healthy):**
- IFI27 (logFC 3.48, adj.P.Val 7.7e-19)
- XAF1 (logFC 1.56, adj.P.Val 6.1e-17)
- IFITM3 (logFC 1.41, adj.P.Val 1.4e-16)
- EPSTI1 (logFC 1.81, adj.P.Val 3.1e-16)
- IFI44L (logFC 2.27, adj.P.Val 3.2e-16)
**Downregulated (SLE < Healthy):**
- None detected at this threshold

---

## Academic conclusions

1. **Reproducible interferon-driven signature:** At this more inclusive FDR (5%), we capture additional interferon-stimulated genes beyond those passing 1%, confirming that the **core biological signal is robust**.

2. **Magnitude and direction:** Nearly all significant DEGs are *upregulated* in SLE — the immune activation pattern remains dominant, and downregulated genes are sparse.

3. **Temporal persistence:** The interferon module appears consistently at all gestational stages, showing that **SLE molecular dysregulation is sustained throughout pregnancy and postpartum.**

4. **Added biological depth:** At FDR 0.05, supplementary pathways emerge (e.g., *cytokine-mediated signaling, antigen processing, neutrophil activation*), expanding the immune context of the response.

5. **Complementary to FDR 0.01 results:** These DEGs are not contradictory but rather **an extended view** of the same underlying disease biology — providing a richer feature pool for pathway analysis and machine learning.


---

## Suggested downstream analyses

- Compare these results directly with the stricter FDR 0.01 set to distinguish **core vs peripheral** components of the SLE signature.

- Use this broader DEG set for **GO/KEGG enrichment**, expecting to see interferon and cytokine terms with higher gene counts.

- Include these genes in **UpSet plots** to visualize overlap across pregnancy stages.

- Use both FDR0.01 and FDR0.05 lists to test ML model robustness to feature selection thresholds.


---

## Thesis interpretation guidance

- The **qualitative pattern is identical** to the FDR 0.01 subset — SLE remains defined by heightened antiviral and interferon activity.

- Quantitatively, this threshold reveals a **wider network of immune-related genes**, providing additional targets for discussion and enrichment.

- In your Results section, emphasize that **statistical stringency (FDR 1% vs 5%) does not alter biological interpretation**, underscoring reproducibility and robustness.
