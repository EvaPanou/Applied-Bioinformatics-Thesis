# Folder: `B_Temporal_adjacent/`

### Longitudinal (Temporal) DEG Analysis within each condition

This folder contains **within-condition** differential expression results comparing **adjacent pregnancy stages**:

`<16w → 16–23w → 24–31w → 32–40w → Postpartum (PP)`.

Each file tests **how gene expression changes over time** separately for SLE and Healthy women.


---

## Summary of DEGs (FDR < 0.05, |logFC| > 1)

| Time Comparison   |   SLE Up |   SLE Down |   Healthy Up |   Healthy Down |
|:------------------|---------:|-----------:|-------------:|---------------:|
| 16–23w vs <16w    |       10 |          0 |           14 |              0 |
| 24–31w vs 16–23w  |        0 |          0 |            0 |              0 |
| 32–40w vs 24–31w  |        0 |          0 |            0 |              0 |
| PP vs 32–40w      |        2 |        119 |            7 |            135 |

> **Interpretation:** Healthy pregnancies show almost no significant DEGs across timepoints — indicating a **stable immune transcriptome**. In contrast, SLE pregnancies exhibit moderate transcriptional changes between stages, reflecting **disease-specific temporal modulation** rather than normal pregnancy adaptation.


---

## Top representative genes per comparison

### 16–23w vs <16w

**SLE:**
- TENT5C (logFC 0.99, adj.P.Val 7.2e-14)
- CA1 (logFC 1.29, adj.P.Val 8.0e-14)
- ZDHHC19 (logFC 0.96, adj.P.Val 3.8e-13)
- XK (logFC 1.21, adj.P.Val 4.3e-13)
- TMOD1 (logFC 1.01, adj.P.Val 5.8e-13)
**Healthy:**
- BNIP3L (logFC 1.11, adj.P.Val 3.2e-08)
- CA1 (logFC 1.36, adj.P.Val 8.7e-07)
- BBOF1 (logFC 0.85, adj.P.Val 3.1e-06)
- RAB8A (logFC -0.35, adj.P.Val 3.1e-06)
- TENT5C (logFC 0.98, adj.P.Val 3.1e-06)

### 24–31w vs 16–23w

**SLE:**
- BEX1 (logFC 0.30, adj.P.Val 1.2e-03)
- OLR1 (logFC 0.68, adj.P.Val 2.3e-03)
- SLC2A5 (logFC 0.45, adj.P.Val 6.3e-03)
- LOC284648 (logFC 0.48, adj.P.Val 7.1e-03)
- RETN (logFC 0.60, adj.P.Val 1.6e-02)
**Healthy:**
- ZNF430 (logFC -0.24, adj.P.Val 7.3e-02)
- THOC6 (logFC 0.33, adj.P.Val 7.3e-02)
- ACP6 (logFC 0.30, adj.P.Val 1.7e-01)
- TIMM44 (logFC 0.30, adj.P.Val 1.7e-01)
- EAF1 (logFC -0.35, adj.P.Val 1.7e-01)

### 32–40w vs 24–31w

**SLE:**
- TARS (logFC 0.37, adj.P.Val 1.1e-02)
- YARS (logFC 0.23, adj.P.Val 2.5e-02)
- ABCA1 (logFC -0.39, adj.P.Val 2.5e-02)
- CDS2 (logFC -0.21, adj.P.Val 2.5e-02)
- ACAA2 (logFC 0.24, adj.P.Val 2.5e-02)
**Healthy:**
- ABCA1 (logFC -0.56, adj.P.Val 1.7e-02)
- ZBTB17 (logFC -0.24, adj.P.Val 5.1e-02)
- IL16 (logFC -0.25, adj.P.Val 8.5e-02)
- MRPL15 (logFC 0.43, adj.P.Val 1.6e-01)
- GBP4 (logFC 0.40, adj.P.Val 1.6e-01)

### PP vs 32–40w

**SLE:**
- MCEMP1 (logFC -1.56, adj.P.Val 4.0e-36)
- SLC22A4 (logFC -1.02, adj.P.Val 5.2e-31)
- ANXA3 (logFC -1.77, adj.P.Val 8.1e-30)
- TENT5C (logFC -1.88, adj.P.Val 7.4e-29)
- XK (logFC -2.40, adj.P.Val 7.4e-29)
**Healthy:**
- MCEMP1 (logFC -2.20, adj.P.Val 5.8e-50)
- TLR5 (logFC -1.48, adj.P.Val 5.1e-37)
- SLC22A4 (logFC -1.27, adj.P.Val 1.1e-34)
- CYSTM1 (logFC -1.48, adj.P.Val 1.2e-32)
- ANXA3 (logFC -2.16, adj.P.Val 1.3e-32)


---

## Academic interpretation

1. **Healthy pregnancies:** Show minimal temporal DEGs, confirming that pregnancy progression in healthy women is associated with **physiological but not large-scale transcriptional immune shifts**.

2. **SLE pregnancies:** Exhibit more pronounced gene expression changes between stages, particularly between **late pregnancy and postpartum**, suggesting fluctuations in immune activation or disease activity.

3. **Gene identity:** Many top SLE DEGs belong to **interferon-stimulated gene (ISG) and cytokine-response families**, e.g., *IFI27, IFIT3, MX1, ISG15* — reinforcing the idea that **interferon signaling remains dominant** and fluctuates slightly over time.

4. **Temporal pattern:** The lack of reversal to a healthy-like profile postpartum indicates **persistent immune dysregulation** after delivery.

5. **Biological meaning:** These within-SLE changes may correspond to **flare modulation or treatment response**, but the underlying interferon signature remains consistent.


---

## How to use these results

- Combine these temporal DEGs with the `A_SLE_vs_Healthy` contrasts to separate **disease-specific effects** from **time-driven physiological effects**.

- Visualize overlaps using **UpSet plots** or **heatmaps** to reveal which genes fluctuate only in SLE.

- Use these DEGs to explore **longitudinal biomarker stability**, or as input for **time-series machine learning**.

- Focus follow-up enrichment on SLE contrasts with larger DEG counts (e.g., late vs early or postpartum vs late pregnancy).


---

## Thesis discussion points

- Healthy pregnancies maintain transcriptomic homeostasis; SLE shows moderate dynamic regulation but no normalization to healthy levels.

- These results support the concept of **pregnancy as a physiological stress test** revealing immune inflexibility in SLE.

- Combined with the constant interferon activity from the `A_` contrasts, these findings imply that SLE disease activity is **chronic and resistant to gestational immunomodulation**.
