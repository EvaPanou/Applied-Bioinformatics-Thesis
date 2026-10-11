# 04_ML_Prerequisites

## Folder Overview

Builds the machine-learning-ready input table from the 44-gene `Method_All_Union` DEG panel (from `02_DGE_Analysis`) and the corrected sample metadata. Pure merge/formatting only — no model training, no CV strategy, no scaling, no derived/composite features. Those belong to the separate, later `05_ML_Pipeline` stage.

## Pipeline context

```text
02_DGE_Analysis
      ↓
03_Feature_Validation
      ↓
04_ML_Prerequisites   <- this stage
      ↓
05_ML_Pipeline
```

## Folder Structure & File Reference

Per the repo-wide convention, `00_` files are inputs, odd numbers are code, even numbers are outputs:

| File | Type | Description |
|---|---|---|
| `00_ML_FeatureMatrix_Method_All_Union_log2.tsv` | Input | Log2 expression matrix, 475 samples × 44 genes (`Method_All_Union` panel, ≥1 of 3 DE methods) — copied from `02_DGE_Analysis/02_Results/DE_Genes/`. |
| `00_non_NP_metadata.csv` | Input | Corrected sample metadata, 489 non-NP samples (pre-outlier-removal) — used for `SampleID`, `DonorID`, `Condition`, `Timepoint`, and their readable labels. |
| `01_ML_Prerequisites_Notebook.ipynb` | Code | Loads both inputs, merges, validates, and saves the ML-ready table — 5 short steps, no functions, no custom logic beyond a rename/merge/reorder. |
| `02_ML_Input_File.tsv` | Output | The final merged table — 475 samples × 50 columns (6 metadata columns + 44 gene columns) — the input to `05_ML_Pipeline`. |

---

## Methodology

**What the notebook does**, step by step:

1. **Load** both inputs; print shapes to confirm before doing anything else (475×44 for the expression matrix, 489×29 for the metadata).
2. **Select and rename metadata columns** — `Sample`→`SampleID`, `Donor_id`→`DonorID`, `sle`→`Condition` (0/1), `Condition`→`Condition_label`, `tp`→`Timepoint` (ordinal 1–5), `time_point`→`Time_label`. Only these 6 columns are kept from the 29-column metadata file.
3. **Merge** metadata (489 rows) with the expression matrix (475 rows) on `SampleID`, via an **inner join**, then reorder columns (metadata first, then the 44 gene columns unmodified). The inner join is what actually drops the 14 samples flagged as outliers in `01_RawData_&_PCA` — those 14 `SampleID`s exist in the metadata but have no matching row in the already-outlier-removed expression matrix, so the join keeps only IDs present in both tables. No separate filtering step is needed or written.
4. **Validate**: prints final sample count (475), gene count (44), the `Condition_label` distribution, and unique donor count — a quick sanity check before saving, not a formal test.
5. **Save** as `02_ML_Input_File.tsv`.

**How to run it:** All four input/output paths in the notebook are plain filenames (no hardcoded absolute paths, unlike the R scripts in earlier stages) — it assumes both `00_` input files sit in the same directory as the notebook. Run top to bottom in Jupyter.

**Environment:** Python (Anaconda-managed local kernel, consistent with the rest of the thesis). Packages: `pandas` only. No version pinned in the notebook itself.

---

## A note on this stage's simplicity

Unlike every other code artifact in this repo, this notebook has no bugs documented, no methodological caveats, and no custom functions — it's a straightforward rename/merge/reorder/save, five steps, no branching logic. That's a deliberate design choice stated in the folder's own original description: this stage is explicitly "pure merge/formatting only," with modeling decisions (splitting strategy, feature scaling, derived features) all deferred to `05_ML_Pipeline`, which is a separate, later stage built to read this file's output directly.

## Important Consideration

This is a longitudinal dataset — multiple samples per donor across gestational timepoints. `DonorID` is retained specifically so `05_ML_Pipeline` can split train/test **by donor**, not by sample, to avoid data leakage from repeated measures of the same individual — the `DonorID` column carried through here is what makes that possible downstream.

## Verified Output

`02_ML_Input_File.tsv` — 475 rows × 50 columns confirmed directly. Header and first row:

```text
SampleID,DonorID,Condition,Condition_label,Timepoint,Time_label,HES4,HERC6,TRIM6,RSAD2,...
GSM2901849,JC12,0,Healthy,1,<16 weeks,7.81,6.88,4.95,5.94,...
```

**Used downstream by:** `05_ML_Pipeline`, which reads this exact file (copied forward as `00_ML_Input_File.tsv` in that stage's own folder) as its sole starting point — every gene available to the donor-aware benchmarking, feature selection, and final-panel derivation in that stage comes from these 44 columns.
