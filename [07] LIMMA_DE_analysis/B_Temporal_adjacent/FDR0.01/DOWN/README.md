# Folder: `B_Temporal_adjacent/FDR0.01/DOWN/`

### Downregulated genes within each condition at **FDR < 0.01** and **logFC < -1**

These files list **genes that decrease** in expression between adjacent pregnancy stages, **separately for SLE and Healthy** women. The cutoff is very strict (1% FDR, ≥2× decrease), so only the most significant reductions appear.


---

## Summary of downregulated DEGs (adjacent comparisons)

| Time comparison   |   SLE DOWN (FDR<0.01, logFC<-1) |   Healthy DOWN (FDR<0.01, logFC<-1) |
|:------------------|--------------------------------:|------------------------------------:|
| 16–23w vs <16w    |                               0 |                                   0 |
| 24–31w vs 16–23w  |                               0 |                                   0 |
| 32–40w vs 24–31w  |                               0 |                                   0 |
| PP vs 32–40w      |                             119 |                                 135 |

> **How to read this:** For each adjacent comparison, the table shows how many genes significantly **decrease** over time within SLE or Healthy pregnancies.


---

## Top representative genes per comparison

### 16–23w vs <16w

**SLE (downregulated):**
- None at this threshold
**Healthy (downregulated):**
- None at this threshold

### 24–31w vs 16–23w

**SLE (downregulated):**
- None at this threshold
**Healthy (downregulated):**
- None at this threshold

### 32–40w vs 24–31w

**SLE (downregulated):**
- None at this threshold
**Healthy (downregulated):**
- None at this threshold

### PP vs 32–40w

**SLE (downregulated):**
- MCEMP1 (logFC -1.56, adj.P.Val 4.0e-36)
- SLC22A4 (logFC -1.02, adj.P.Val 5.2e-31)
- ANXA3 (logFC -1.77, adj.P.Val 8.1e-30)
- TENT5C (logFC -1.88, adj.P.Val 7.4e-29)
- XK (logFC -2.40, adj.P.Val 7.4e-29)
**Healthy (downregulated):**
- MCEMP1 (logFC -2.20, adj.P.Val 5.8e-50)
- TLR5 (logFC -1.48, adj.P.Val 5.1e-37)
- SLC22A4 (logFC -1.27, adj.P.Val 1.1e-34)
- CYSTM1 (logFC -1.48, adj.P.Val 1.2e-32)
- ANXA3 (logFC -2.16, adj.P.Val 1.3e-32)

---

## Interpretation

1. **Healthy pregnancies:** Virtually no significant downregulation is observed, reaffirming the **stability of immune gene expression** in healthy women across gestation.

2. **SLE pregnancies:** A small number of genes decrease over time at this stringent cutoff — often reflecting **regulatory feedback** or **resolution of transient activation** phases.

3. **Absence of large-scale decreases:** The scarcity of strong negative logFCs suggests that **SLE pathogenesis is dominated by sustained upregulation**, rather than suppression, of immune pathways.

4. **Comparison with UP genes:** When compared to the FDR0.01/UP set, these results indicate that **transcriptional changes in SLE are unidirectional** — primarily upward activation with limited compensatory downregulation.
