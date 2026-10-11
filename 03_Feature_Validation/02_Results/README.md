# 03_Feature_Validation / 02_Results

## Folder Overview

All output from `../01_Feature_Validation.R`. Four subfolders, one per output type — this folder holds no files of its own directly.

## Folder Structure & File Reference

| Folder | Description |
|---|---|
| [`MA_Plots/`](./MA_Plots/README.md) | MA plots (Average Expression vs. logFC), one per after-panel (2 files) |
| [`PCA_Plots/`](./PCA_Plots/README.md) | PCA plots, one per panel including the before baseline (3 files) |
| [`PERMANOVA/`](./PERMANOVA/README.md) | The formal Condition/Time effect-size test, one combined summary table (1 file) |
| [`Volcano_Plots/`](./Volcano_Plots/README.md) | Volcano plots (logFC vs. -log10 FDR), one per after-panel (2 files) |

---

## Results — the headline comparison

The single most informative number this stage produces is in `PERMANOVA/PERMANOVA_Summary.tsv`: Condition's PERMANOVA R² at the donor level goes from 6.1% (`All_Genes`, before any gene selection) to 41.1% (`DEG_Union`, 44 genes) to 41.6% (`Final_Panel`, 42 genes) — see that subfolder's README for the full table and what each column means. Every other output in this folder — the PCA ellipse separation, the volcano/MA spread — is the visual counterpart to that same underlying result.
