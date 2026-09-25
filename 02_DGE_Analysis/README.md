# Differential Gene Expression Analysis

## Folder Overview

Differential gene expression (DGE) analysis stage of the SLE pregnancy biomarker discovery project. Takes the filtered, NP-donor-corrected expression data from `01_RawData_&_PCA` and identifies genes that reliably distinguish SLE from Healthy, using three independent statistical methods rather than a single test — the pipeline's own recurring theme of cross-method validation, echoed later in `03_Feature_Validation` and again in `01_RawData_&_PCA/04_ALASCA_Output`'s independent multivariate confirmation of the same signal.

## Folder Structure & File Reference

Per the repo-wide convention, `00_` files are inputs, odd numbers are code, even numbers are outputs:

| File / folder | Type | Description |
|---|---|---|
| `00_final_filtered_expression_non_NP_donors.tsv` | Input | The cleaned, filtered expression matrix from `01_RawData_&_PCA/02_Analysis_Results` — 9,055 genes × 475 samples, outliers already removed. |
| `00_non_NP_metadata.csv` | Input | Metadata for the 489-sample non-NP cohort (489, not 475 — the expression matrix has already dropped the 14 outlier samples that this metadata file still includes; the script aligns the two by sample ID rather than assuming equal row counts). |
| `01_limma_DEG_pipeline.R` | Code | Loads both inputs, fits three independent differential-expression models (Methods A, B, C), pools Method A's per-timepoint results via meta-analysis, cross-references all three methods into combined gene panels, and produces every downstream table/figure — 13 steps, wrapped in one function (`run_deg_pipeline()`), called once at the very end of the file. |
| `02_Results/` | Output (folder) | Everything the script produces: 5 subfolders plus 6 top-level files. See [`02_Results/README.md`](./02_Results/README.md). |

---

## Methodology, in detail

Three independent statistical methods test SLE vs. Healthy with the same significance rule (**FDR < 0.05, \|logFC\| > 1**, ~2-fold change), differing only in how the effect size and repeated-measures structure (multiple timepoints per donor) are estimated:

- **Method A** (Steps 2–4) — a sample-means design (`~0 + Condition:Time`, no shared intercept) tests SLE vs. Healthy separately within each of the 5 gestational timepoints, using limma's `duplicateCorrelation()` + donor blocking to correct for repeated measures. Run 5 times, this directly shows whether the signal is stable across gestation or concentrated at one stage, rather than assuming stability.
- **Method B** (Step 5) — a single pooled model across all timepoints (`~0 + Condition`), also via `duplicateCorrelation()`. This estimates *one* consensus within-donor correlation across the whole genome and applies it to every gene equally — it corrects the variance estimate for repeated measures, but doesn't rebalance each donor's influence on the group mean (a donor with 5 samples still pulls 5× harder on their Condition's average than a donor with 1 sample).
- **Method C** (Step 6) — the same pooled question as Method B, but via `variancePartition::dream()`, a true linear mixed model fit per gene (REML), estimating the within-donor correlation *separately for each gene* rather than one shared genome-wide number, and adjusting degrees of freedom accordingly.

**Method A pooling** (Step 9) uses random-effects meta-analysis (`metafor::rma()`, REML, inverse-variance weighted) across the 5 timepoint estimates, rather than simple vote-counting (a documented weak method — Light & Smith 1971; Hedges & Olkin 1980) or naive averaging. I² and Cochran's Q (QEp) are reported per gene to quantify heterogeneity across timepoints — a gene with a strong pooled effect *and* low heterogeneity is a good time-invariant candidate for a marker meant to work regardless of when in pregnancy blood is drawn.

**Final gene panels** (Step 10) are built by cross-referencing which genes reach significance in each of the three methods — union, intersection, and combinations thereof — rather than trusting any single method alone. See [`02_Results/DE_Genes/README.md`](./02_Results/DE_Genes/README.md) for the full naming convention and every resulting file.

**A formal interaction test** (Step 12) directly tests whether the Condition effect itself changes shape across gestation — a joint F-test across all `Condition:Time` interaction coefficients simultaneously, following the temporal-DEG axis from TiSA (Lefol et al. 2023). This is a different, more rigorous question than anything Methods A–C individually answer, since Method A's per-timepoint estimates and the meta-analysis's per-gene heterogeneity stats don't add up to one combined significance test across the whole interaction term. **Result: 0 of 9,055 genes reach FDR < 0.05** — see `Global_Time_by_Condition_Ftest.tsv`.

**Two real bugs were caught and fixed during this script's development** (both documented directly in code comments, and confirmed against the actual `_Rhistory` console log uploaded alongside this stage):
1. `dream()`'s `L=` contrast argument keeps the raw per-group cell-means coefficients (`ConditionHealthy`, `ConditionSLE`) alongside the requested contrast column, rather than replacing them. An earlier version of the script looped over all 3 resulting columns as if they were all meaningful contrasts — but `ConditionHealthy`/`ConditionSLE` each test whether that group's raw average log-expression differs from *zero*, which is true for almost every gene on a log scale. This silently inflated Method C's "significant" gene count to nearly the entire dataset (confirmed: ~9,055 of 9,055). Fixed by explicitly selecting only the `SLE_vs_Healthy_Dream` column.
2. The global F-test's output table was originally saved without a `Gene` column — `topTable()`'s row names (the gene symbols) aren't written to a TSV by `fwrite()` unless explicitly added as a column first, unlike every other `topTable()` call in the script, which already did this. Fixed by adding `interaction_top_table$Gene <- rownames(interaction_top_table)` before saving.

**How to run it:** Open in RStudio. **Before running**, edit the hardcoded absolute path at the top (`setwd("C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/02_DGE_Analysis")`) and the two `EXPR_FILE`/`META_FILE` paths built from it — same local-machine-specific path issue as `01_RawData_&_PCA/03_ALASCA_Analysis.R`, and again referencing an older folder name (`Experiment/02_DGE_Analysis`) that predates this repo's current structure. Source the whole file — this only *defines* `run_deg_pipeline()` and loads packages, it doesn't run the analysis until the final line (`run_deg_pipeline()`) executes, which happens automatically when you source the file top to bottom. Comment out that last line if you want the function available without immediately running all 13 steps.

**Environment:** R, via RStudio. Packages: `limma`, `ComplexHeatmap`, `circlize`, `UpSetR`, `tidyverse`, `data.table`, `matrixStats`, `metafor`, `variancePartition` (for `dream()`), `BiocParallel`. No R or package version numbers are pinned in the script. `BiocParallel` is deliberately registered as `SerialParam` (one gene at a time, with a progress bar) rather than the faster `SnowParam` parallel backend — the script's own comment notes that `SnowParam` hung indefinitely with no warning on a prior run (suspected Windows/RStudio-specific issue), and recommends only trying it after a serial run has completed successfully at least once.

**Runtime note:** `dream()` fits one mixed model per gene (9,055 genes) and is explicitly flagged in the script as taking substantially longer than Methods A or B.

---

## Key Result

**42 genes** are found significant by all three methods simultaneously (`Method_All_Intersection` in `02_Results/DE_Genes/`) — a dominant, overwhelmingly up-regulated interferon-stimulated gene (ISG) signature (`IFI44L`, `MX1`, `OAS1/2/3`, `ISG15`, `RSAD2`, and others), confirmed stable across gestation by the genome-wide interaction test (0 of 9,055 genes significant).

## Next Stage

`03_Feature_Validation` validates this DEG feature set (union vs. intersection) via PCA, PERMANOVA, and volcano/MA plots — see that stage's own README for which panel it carries forward.
