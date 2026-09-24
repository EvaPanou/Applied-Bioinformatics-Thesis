# 01_RawData_&_PCA / 02_Analysis_Results

## Folder Overview

Every table and figure produced by `../01_RawData_Processing.ipynb`: sample-level outlier detection (two independent methods), gene-level filtering, and exploratory PCA on the raw microarray expression matrix. This folder is pure output — no code, no inputs — numbered in the order the notebook generates them.

## Folder Structure

As with the equivalent metadata-stage folder, this doesn't follow the 00/odd/even convention (there's no code or input here, only outputs from `../01_RawData_Processing.ipynb`), and files are numbered in generation order:

```text
02_Analysis_Results/
├── 01_boxplot_raw_expression_all_donors.png
├── 02_metadata_NP_only_locked.tsv
├── 03_expression_NP_only_locked.tsv
├── 04_outlier_samples_non_NP_donors.tsv
├── 05_extreme_outlier_samples_non_NP_donors.tsv
├── 06_pca_outlier_check_non_NP_donors.png
├── 07_samples_flagged_by_both_methods_non_NP_donors.tsv
├── 08_intersection_outliers_non_NP_donors.png
├── 09_expression_no_outliers_non_NP_donors.tsv
├── 10_boxplot_no_outliers_non_NP_donors.png
├── 11_gene_variance_distribution_non_NP_donors.png
├── 12_all_genes_backup_no_expression_filter_non_NP_donors.tsv
├── 13_low_expressed_genes_removed_non_NP_donors.tsv
├── 14_genes_kept_after_expression_filter_non_NP_donors.tsv
├── 15_final_filtered_expression_non_NP_donors.tsv
├── 16_histogram_final_expression_non_NP_donors.png
├── 17_pca_sample_level_by_condition_non_NP_donors.png
├── 18_pca_sample_level_by_timepoint_non_NP_donors.png
└── 19_metadata_for_alasca_non_NP_donors.tsv
```

## File Reference

| File | Description |
|---|---|
| `01_boxplot_raw_expression_all_donors.png` | Per-sample expression boxplot, raw data, all 512 samples, before any cleaning |
| `02_metadata_NP_only_locked.tsv` | Metadata for the 23 never-pregnant donors, locked aside |
| `03_expression_NP_only_locked.tsv` | Expression data for the same 23 NP donors, locked aside |
| `04_outlier_samples_non_NP_donors.tsv` | Samples beyond 1.5× IQR of the sample-median distribution |
| `05_extreme_outlier_samples_non_NP_donors.tsv` | Samples beyond 3× IQR (the "extreme" tier) |
| `06_pca_outlier_check_non_NP_donors.png` | PCA + 95% Mahalanobis/χ² ellipse outlier plot |
| `07_samples_flagged_by_both_methods_non_NP_donors.tsv` | The intersection: extreme-IQR **and** PCA-flagged — the actual removal set |
| `08_intersection_outliers_non_NP_donors.png` | Same PCA plot, highlighting only the intersection set |
| `09_expression_no_outliers_non_NP_donors.tsv` | Expression matrix with the 14 intersection outliers removed |
| `10_boxplot_no_outliers_non_NP_donors.png` | Per-sample boxplot after outlier removal |
| `11_gene_variance_distribution_non_NP_donors.png` | Gene variance histogram with the 10th-percentile filter cutoff marked |
| `12_all_genes_backup_no_expression_filter_non_NP_donors.tsv` | Variance-filtered matrix, before the low-expression filter — kept as a comparison baseline |
| `13_low_expressed_genes_removed_non_NP_donors.tsv` | The 3,019 genes removed by the low-expression filter (≤25th percentile mean) |
| `14_genes_kept_after_expression_filter_non_NP_donors.tsv` | The 9,055 genes retained |
| `15_final_filtered_expression_non_NP_donors.tsv` | The final cleaned, filtered expression matrix |
| `16_histogram_final_expression_non_NP_donors.png` | Value distribution of the final matrix |
| `17_pca_sample_level_by_condition_non_NP_donors.png` | Sample PCA, filtered vs. unfiltered, colored by Condition |
| `18_pca_sample_level_by_timepoint_non_NP_donors.png` | Same comparison, colored by Timepoint |
| `19_metadata_for_alasca_non_NP_donors.tsv` | Minimal metadata (`Donor_id`, `time_point`, `Condition`) for the surviving 475 samples, prepared for the R/ALASCA script |

---

## Results

### Raw data check (Step 4)

**`01_boxplot_raw_expression_all_donors.png`** — all 512 samples, before any cleaning:

![Boxplot raw expression](01_boxplot_raw_expression_all_donors.png)

Values sit in roughly the 4–14 range across all samples, consistent with data that's already been log2-transformed and normalized upstream by ADEx — no obvious gross technical failures visible at this stage (no sample sitting wildly off from the rest).

### NP donors locked aside (Step 5)

**`02_metadata_NP_only_locked.tsv`** — 23 rows, same 29-column schema as the parent metadata files. Header and first row:

```text
Sample,GSE,Experimental Strategy,GPL,Condition,Tissue,Gender,Age Group,Age,Race,Ethnicity,grp_p_tp,Sample_name,Donor_id,sle,apl,lac,tp,pe,fd,nnd,pl_insuff,iugr,sga,batch,time_point,ga_at_collection,ga_at_end_of_pregnancy,if_pe_before_or_after_36_weeks
GSM2901826,GSE108497,Expression,GPL10558,Healthy,Whole blood,Female,21-30,25,C,Not Hispanic or Latino,HC_NP_5,HC2013_1,106346,0,0,0,5,0,0,0,0,0,0,3,NP,,,
```

**`03_expression_NP_only_locked.tsv`** — 13,416 genes × 23 samples, unfiltered (variance/expression filtering downstream only applies to the non-NP working set). This and the file above are held for the specificity checks used later in `06_Final_Panel_Validation`, not touched again in this notebook.

### Outlier detection (Step 6)

**`04_outlier_samples_non_NP_donors.tsv`** — 53 samples flagged by the loose 1.5×IQR rule (out of 489 non-NP samples — about 11%).

**`05_extreme_outlier_samples_non_NP_donors.tsv`** — 20 samples flagged by the stricter 3×IQR "extreme" rule.

**`06_pca_outlier_check_non_NP_donors.png`** — the PCA/Mahalanobis view of the same question:

![PCA-based outlier check](06_pca_outlier_check_non_NP_donors.png)

A cluster of tightly-grouped samples (within the 95% ellipse) with a scattered set of clear outliers extending out along PC1 — the flagged points aren't borderline, they sit well outside the main cloud.

**`07_samples_flagged_by_both_methods_non_NP_donors.tsv`** — the actual removal set: only samples flagged as *extreme* by the IQR method **and** flagged by PCA. **14 samples** — a deliberately conservative intersection rather than the union of either method alone (53 or 20 samples individually).

**`08_intersection_outliers_non_NP_donors.png`** — same PCA view, now highlighting only those 14:

![Intersection of IQR and PCA outliers](08_intersection_outliers_non_NP_donors.png)

All 14 are the most extreme points from the plot above — visually, the intersection rule doesn't appear to be dropping anything borderline.

**`09_expression_no_outliers_non_NP_donors.tsv`** — 475 samples (489 − 14) × 13,416 genes, unfiltered at the gene level. **This sample count (475) is the one that persists through every later stage** — it's the final analytic sample size referenced in `02_DGE_Analysis`, `03_Feature_Validation`, and `04_ML_Prerequisites`.

### Post-removal check & gene filtering (Steps 7–10)

**`10_boxplot_no_outliers_non_NP_donors.png`**:

![Boxplot after outlier removal](10_boxplot_no_outliers_non_NP_donors.png)

Visually near-identical to the pre-removal boxplot — expected, since only 14 of 489 samples were removed and none were egregious enough to visibly skew the overall per-sample distributions in the first plot either.

**`11_gene_variance_distribution_non_NP_donors.png`** — variance filter cutoff (10th percentile):

![Gene variance distribution](11_gene_variance_distribution_non_NP_donors.png)

A heavily right-skewed distribution — most genes have low variance, a long tail of a few genes with much higher variance. Removing the bottom 10% takes the gene count from 13,416 to 12,074.

**`12_all_genes_backup_no_expression_filter_non_NP_donors.tsv`** — 475 samples × 12,074 genes. Kept as a comparison baseline against the more heavily filtered matrix below; referenced again in Step 12's filtered-vs-unfiltered PCA comparison.

**`13_low_expressed_genes_removed_non_NP_donors.tsv`** — 3,019 genes removed for low mean expression (≤25th percentile). Lowest-expressed removed genes include `CRISP3`, `APOBEC3B`, `SNORD3C`.

**`14_genes_kept_after_expression_filter_non_NP_donors.tsv`** — 9,055 genes retained. Highest-expressed kept genes include `HBA1`, `HBG2`, `SRGN` — hemoglobin and blood-cell-granule genes, as expected for whole-blood microarray data. **3,019 + 9,055 = 12,074**, confirming the two lists partition the variance-filtered gene set exactly.

**`15_final_filtered_expression_non_NP_donors.tsv`** — 475 samples × 9,055 genes. **This is the expression matrix `02_DGE_Analysis` builds its differential expression analysis on directly.**

### Final checks & exploratory PCA (Steps 11–12)

**`16_histogram_final_expression_non_NP_donors.png`**:

![Histogram of final filtered expression values](16_histogram_final_expression_non_NP_donors.png)

A right-skewed but unimodal distribution peaking around 5.5–6, tailing out toward 14 — the expected shape for filtered log2 microarray intensities.

**`17_pca_sample_level_by_condition_non_NP_donors.png`** — filtered vs. unfiltered matrix, colored by Condition:

![Sample-level PCA by Condition](17_pca_sample_level_by_condition_non_NP_donors.png)

SLE and Healthy samples show heavy overlap in both panels, with SLE (orange) spreading out somewhat further along PC1 in both — expected at this stage, since this is whole-transcriptome variance, not the DEG-restricted signal that later stages isolate. The filtered and unfiltered panels look broadly similar, suggesting the low-expression filter isn't distorting the overall sample-level structure.

**`18_pca_sample_level_by_timepoint_non_NP_donors.png`** — same comparison, colored by Timepoint:

![Sample-level PCA by Timepoint](18_pca_sample_level_by_timepoint_non_NP_donors.png)

No visually obvious timepoint clustering or trajectory in either panel — timepoints are thoroughly intermixed, consistent with `02_DGE_Analysis`'s later finding of no significant Condition × Time interaction. This whole-transcriptome PCA isn't really the tool to detect a subtle DEG-level effect either way, so this is a sanity check rather than a definitive result.

### ALASCA handoff (Step 13)

**`19_metadata_for_alasca_non_NP_donors.tsv`** — 475 rows, 4 columns (`Sample`, `Donor_id`, `time_point`, `Condition`). Header and first row:

```text
Sample,Donor_id,time_point,Condition
GSM2901849,JC12,<16 weeks,Healthy
```

**Used directly by `../03_ALASCA_Analysis.R`**, alongside `15_final_filtered_expression_non_NP_donors.tsv` — together, these are the two files the R script reads to build its long-format input.
