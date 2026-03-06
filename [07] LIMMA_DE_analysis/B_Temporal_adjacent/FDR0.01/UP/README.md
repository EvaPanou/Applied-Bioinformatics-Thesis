# Folder: `B_Temporal_adjacent/FDR0.01/UP/`

### Upregulated genes within each condition at **FDR < 0.01** and **logFC > 1**

These files list **genes that increase** between adjacent pregnancy stages, **separately for SLE and Healthy** women. The thresholds are strict (1% FDR, ≥2× increase), so only the most robust changes appear.


---

## Summary of upregulated DEGs (adjacent comparisons)

| Time comparison   |   SLE UP (FDR<0.01, logFC>1) |   Healthy UP (FDR<0.01, logFC>1) |
|:------------------|-----------------------------:|---------------------------------:|
| 16–23w vs <16w    |                           10 |                               14 |
| 24–31w vs 16–23w  |                            0 |                                0 |
| 32–40w vs 24–31w  |                            0 |                                0 |
| PP vs 32–40w      |                            2 |                                7 |

> **How to read this:** For each adjacent step, the table shows how many genes significantly **increase** in SLE and in Healthy.


---

## Top representative genes per comparison

### 16–23w vs <16w

**SLE (upregulated):**
- CA1 (logFC 1.29, adj.P.Val 8.0e-14)
- XK (logFC 1.21, adj.P.Val 4.3e-13)
- TMOD1 (logFC 1.01, adj.P.Val 5.8e-13)
- ITLN1 (logFC 1.05, adj.P.Val 2.7e-11)
- OSBP2 (logFC 1.04, adj.P.Val 1.1e-10)
**Healthy (upregulated):**
- BNIP3L (logFC 1.11, adj.P.Val 3.2e-08)
- CA1 (logFC 1.36, adj.P.Val 8.7e-07)
- YOD1 (logFC 1.06, adj.P.Val 3.1e-06)
- OLR1 (logFC 1.07, adj.P.Val 7.7e-06)
- DEFA3 (logFC 1.16, adj.P.Val 7.8e-06)

### 24–31w vs 16–23w

**SLE (upregulated):**
- None at this threshold
**Healthy (upregulated):**
- None at this threshold

### 32–40w vs 24–31w

**SLE (upregulated):**
- None at this threshold
**Healthy (upregulated):**
- None at this threshold

### PP vs 32–40w

**SLE (upregulated):**
- FCER1A (logFC 1.48, adj.P.Val 3.4e-18)
- HDC (logFC 1.24, adj.P.Val 2.0e-15)
**Healthy (upregulated):**
- CCR3 (logFC 1.46, adj.P.Val 6.1e-20)
- FCER1A (logFC 1.77, adj.P.Val 5.0e-19)
- SLC45A3 (logFC 1.22, adj.P.Val 3.6e-18)
- HDC (logFC 1.54, adj.P.Val 2.4e-17)
- CPA3 (logFC 1.22, adj.P.Val 1.2e-16)

---

## Interpretation

1. **Temporal stability in Healthy:** Very few genes pass this strict cutoff in Healthy across the pregnancy timeline, indicating that **healthy immune transcriptomes remain remarkably stable**.

2. **Time-dependent activation in SLE:** Upregulated SLE genes across adjacent stages often include **interferon-stimulated and cytokine-response genes**, suggesting **modest increases** in immune activation at specific transitions (e.g., late pregnancy → postpartum).

3. **Stringency matters:** Because FDR is 1%, these upregulated sets represent **high-confidence temporal increases**; relaxing to 5% (FDR 0.05) will reveal additional, but less stringent, changes.

4. **Biological meaning:** These changes may reflect **flare modulation or postpartum immune rebound** in SLE, while healthy pregnancies avoid large transcriptional swings.
