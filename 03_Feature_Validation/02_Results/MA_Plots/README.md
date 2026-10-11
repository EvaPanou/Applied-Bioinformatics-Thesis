# 03_Feature_Validation / 02_Results / MA_Plots

## Folder Overview

MA plots (Average Expression vs. log2 Fold Change) for each "after" gene panel, restricted to genes already confirmed significant by `02_DGE_Analysis` — the plots show the spread of effect size across expression levels within an already-validated panel, not a search for new hits. Produced by Step 2 of `../../01_Feature_Validation.R`, from `Method_All_Union_Metrics.tsv`.

## Folder Structure & File Reference

| File | Genes | Description |
|---|---|---|
| `MA_DEG_Union.png` | 44 | MA plot for the lenient `DEG_Union` panel |
| `MA_Final_Panel.png` | 42 | MA plot for the strict `Final_Panel` (cross-method intersection) |

---

## Results

**`MA_DEG_Union.png`**:

![MA plot, DEG Union panel](MA_DEG_Union.png)

All 44 genes sit outside the dashed ±1 log2FC threshold lines, as expected for an already-filtered panel. Genes span roughly 6.5–13 on the average-expression axis and 1–3.4 on the fold-change axis, with the highest-expressed genes (`IFITM3`, `IFIT2`, `MX1`, `IFI6`, `ISG15` — average expression >11) clustering at moderate fold-changes (1.3–2.3), while `IFI27` and `RSAD2` — mid-range in average expression (~9–11) — show the largest fold-changes (~3.0–3.4). **`ORM1` is the one visible point below the lower dashed line** (average expression ~8.7, logFC ≈ -1) — the same single down-regulated gene identified throughout `02_DGE_Analysis`.

**`MA_Final_Panel.png`**:

![MA plot, Final Panel](MA_Final_Panel.png)

The same pattern, minus `ORM1` — since `ORM1` isn't in the 42-gene `Method_All_Intersection` panel (confirmed directly: it doesn't appear in `Method_All_Intersection.tsv`). Every point in this plot sits above the upper threshold line; there is no down-regulated gene in the Final Panel at all.

**Used downstream by:** nothing further in the pipeline — these are end-of-stage visual QC artifacts for the thesis discussion, not an input to `04_ML_Prerequisites`.
