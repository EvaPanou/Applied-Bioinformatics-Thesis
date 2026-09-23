# 04_ML_Prerequisites

## Notes

This stage builds the machine-learning-ready input table from the 44-gene
`Method_All_Union` DEG panel (from `02_DGE_Analysis`) and the corrected
sample metadata. Pure merge/formatting only — no model training, no CV
strategy, no scaling, no derived/composite features. Those belong to a
separate, later ML pipeline stage.

---

## Pipeline context

    02_DGE_Analysis
          ↓
    03_Feature_Validation
          ↓
    04_ML_Prerequisites   <- this stage
          ↓
    (later, separate) Model training / comparison

---

## Folder contents

**Notebook**

- `04_ML_Prerequisites_Notebook.ipynb` — loads inputs, merges, validates,
  and saves the ML-ready table.

**Input files**

- `00_ML_FeatureMatrix_Method_All_Union_log2.tsv` — log2 expression matrix,
  475 samples × 44 genes (`Method_All_Union` panel, ≥1 of 3 DE methods).
- `00_non_NP_metadata.csv` — corrected sample metadata, 489 non-NP samples
  (pre-outlier-removal), used for `SampleID`, `DonorID`, `Condition`,
  `Timepoint`, and their readable labels.

**Output file**

- `ML_Input_File.tsv` — final merged table, 475 samples × 50 columns
  (6 metadata columns + 44 gene columns), used as the input to the next,
  separate ML modeling stage.

---

## Workflow

1. **Load** both input files; confirm shapes.
2. **Select & rename metadata columns** — `Sample`→`SampleID`,
   `Donor_id`→`DonorID`, `sle`→`Condition` (0/1), `Condition`→
   `Condition_label`, `tp`→`Timepoint` (ordinal 1–5), `time_point`→
   `Time_label`.
3. **Merge** metadata with the expression matrix on `SampleID` (inner join)
   and reorder columns (metadata first, then genes). The inner join is
   what drops the 14 samples flagged as outliers in `01_RawData_&_PCA` —
   they exist in the metadata (489 rows) but not in the expression matrix
   (475 rows), so no separate filtering step is needed.
4. **Validate**: sample count, gene count, `Condition_label` distribution,
   unique donor count.
5. **Save** as `ML_Input_File.tsv`.

---

## Important consideration

This is a longitudinal dataset — multiple samples per donor across
gestational timepoints. `DonorID` is retained specifically so the later
modeling stage can split train/test **by donor**, not by sample, to avoid
data leakage from repeated measures of the same individual.
