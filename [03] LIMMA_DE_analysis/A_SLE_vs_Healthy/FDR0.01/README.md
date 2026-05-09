# Folder: `A_SLE_vs_Healthy/FDR0.01/`

### Differentially Expressed Genes (DEGs) at FDR < 0.01 and |logFC| > 1

This folder contains the **high-confidence subset** of the conditional DEG analysis results, comparing **SLE vs Healthy** at each timepoint of pregnancy.

These are genes passing very strict thresholds (False Discovery Rate < 1% and a minimum twofold change in expression).

---

## Summary by timepoint

| Timepoint   |   Upregulated genes (SLE>Healthy) |   Downregulated genes (SLE<Healthy) |
|:------------|----------------------------------:|------------------------------------:|
| <16 weeks   |                                49 |                                   1 |
| 16–23 weeks |                                49 |                                   1 |
| 24–31 weeks |                                44 |                                   0 |
| 32–40 weeks |                                32 |                                   0 |
| Postpartum  |                                40 |                                   0 |

> **Interpretation:** At this stringent cutoff, only the **most reproducible, biologically robust** genes remain. Across all timepoints, SLE shows a strong *up-regulation* signature, while the number of *down-regulated* genes remains minimal.

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

## Conclusions

1. **Persistent interferon activation:** The upregulated genes are overwhelmingly **interferon-stimulated genes (ISGs)** — e.g., *IFI44L, IFI27, IFIT3, RSAD2, ISG15, MX1*. Their presence at every timepoint confirms a **chronic interferon-driven immune activation** in SLE pregnancies.

2. **Temporal consistency:** The signature is observed from <16 weeks through postpartum, indicating that **SLE immune dysregulation persists throughout gestation**.

3. **Directionality:** Almost no strong down-regulated genes appear, highlighting that SLE’s transcriptomic profile reflects **activation, not suppression**, of immune pathways.

4. **Statistical strength:** Because these genes remain significant at FDR < 0.01, they represent **extremely confident findings**, suitable for mechanistic interpretation and biomarker development.

5. **Biological meaning:** These DEGs implicate *Type I interferon signaling, antiviral response, cytokine-mediated signaling,* and *immune effector activation* — pathways central to SLE pathogenesis.
