# 03_Feature_Validation / 02_Results / Volcano_Plots

## Folder Overview

Volcano plots (log2 Fold Change vs. -log10 FDR) for each "after" gene panel, from the same `Method_All_Union_Metrics.tsv` source as the MA plots. Produced by Step 2 of `../../01_Feature_Validation.R`.

## Folder Structure & File Reference

| File | Genes | Description |
|---|---|---|
| `Volcano_DEG_Union.png` | 44 | Volcano plot for the lenient `DEG_Union` panel |
| `Volcano_Final_Panel.png` | 42 | Volcano plot for the strict `Final_Panel` |

---

## Results

**`Volcano_DEG_Union.png`**:

![Volcano plot, DEG Union panel](Volcano_DEG_Union.png)

Significance (-log10 FDR) ranges enormously across the panel — from single digits up to ~85 — reflecting the huge span of statistical confidence within an already-significant gene set: `SPATS2L`, `USP18`, and `ZBP1` sit at the very top (most significant), while `ORM1` sits at the bottom-left, both least significant *and* the panel's only negative fold-change (≈-1, -log10 FDR ≈12) — visually the clear outlier of the panel on both axes at once.

**`Volcano_Final_Panel.png`**:

![Volcano plot, Final Panel](Volcano_Final_Panel.png)

The same shape without `ORM1` — every remaining gene has a positive fold-change, and the same top genes (`SPATS2L`, `USP18`, `ZBP1`) anchor the most-significant end.

**Used downstream by:** nothing further in the pipeline — visual QC artifacts for the thesis discussion.
