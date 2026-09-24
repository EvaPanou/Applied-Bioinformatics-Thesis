# 00_Metadata

## Folder Overview

This is the first stage of the pipeline. It takes the two independent metadata sources for GSE108497 — the ADEx-processed sample table and the raw NCBI GEO series matrix — and merges them into one assembled table (`01_Metadata_Assembly.ipynb`), then explores, corrects, and characterizes that table (`03_Metadata_EDA.ipynb`), producing the corrected sample-level metadata file used by every subsequent stage of the pipeline. This is also where the critical **NP-donor correction** happens: 23 never-pregnant donors, whose `time_point` field was blank in the source data (rather than explicitly marked), are identified and relabeled, then held out as their own group rather than being silently folded into the postpartum (`PP`) timepoint they would otherwise collide with.

## Folder Structure

Per the repo-wide convention, `00_` files are inputs, odd numbers are code, even numbers (from `02_`) are outputs:

```text
00_Metadata/
├── 00_GSE108497_metadata.tsv              [input]  ADEx-processed sample metadata
├── 00_GSE108497_series_matrix.txt         [input]  raw NCBI GEO series matrix
├── 01_Metadata_Assembly.ipynb             [code]   merges the two inputs above
├── 02_GSE108497_assembled_metadata.csv    [output] merged table (output of 01)
├── 03_Metadata_EDA.ipynb                  [code]   corrects, explores, characterizes, splits
├── 04_GSE108497_processed_metadata.csv    [output] corrected full table (output of 03)
└── 04_Analysis_Results/                   [output] all figures/tables from 03 — see its own README
```

## File Reference

| File | Type | Description |
|---|---|---|
| `00_GSE108497_metadata.tsv` | Input | ADEx's curated, tab-separated sample table for GSE108497 — 512 samples × 10 columns (`Sample`, `GSE`, `Experimental Strategy`, `GPL`, `Condition`, `Tissue`, `Cell Type`, `Gender`, `Age`, `Ethnicity`). |
| `00_GSE108497_series_matrix.txt` | Input | Raw NCBI GEO series matrix file for GSE108497 — series-level title/summary/design lines plus per-sample `!Sample_characteristics_ch1` fields (donor ID, disease flags, timepoint, gestational age, etc.). |
| `01_Metadata_Assembly.ipynb` | Code | Merges the two inputs into one table, using the ADEx file as the anchor and pulling additional per-sample characteristics out of the raw GEO file's `!Sample_characteristics_ch1` lines. |
| `02_GSE108497_assembled_metadata.csv` | Output | The merged result — 512 rows × 29 columns. |
| `03_Metadata_EDA.ipynb` | Code | Loads the assembled table, fixes the NP/PP timepoint collision, builds donor-level demographic/clinical tables and figures, runs statistical tests, splits NP donors out, and saves the corrected full table plus the non-NP analytic cohort. |
| `04_GSE108497_processed_metadata.csv` | Output | The corrected full table (512 rows, NP included) — this is the file saved *before* the NP/non-NP split at the end of `03`. |
| `04_Analysis_Results/` | Output (folder) | Every figure and summary table `03` produces, plus the two split files (`GSE108497_NP_stage.csv`, `non_NP_metadata.csv`). See [`04_Analysis_Results/README.md`](./04_Analysis_Results/README.md). |

---

## Methodology

Both notebooks follow the thesis's established code style: plain procedural top-to-bottom code, no custom function definitions, heavily commented.

### `01_Metadata_Assembly.ipynb`

**What it does:** ADEx's metadata file is loaded first (with the `Age` column forced to `str` on read — a value like `"11-20"` was otherwise being silently parsed as a date by pandas). The raw GEO series matrix is then parsed by hand: it isn't tabular, so the notebook scans it line-by-line, pulls the sample ID list from the `!Sample_geo_accession` line, and builds a `{sample_id: {characteristic: value}}` dictionary from every `!Sample_characteristics_ch1` line (handling samples with duplicate characteristic keys by joining them). That dictionary becomes a DataFrame and is left-merged onto the ADEx table on `Sample`, with an assertion that the merge doesn't change the row count (guards against duplicate sample IDs). Columns are then renamed for clarity (ADEx's bucketed `Age` → `Age Group`; GEO's continuous `age` → `Age`; `donor_id` → `Donor_id`, etc.), a handful of redundant/empty columns are dropped, and everything is reordered into a fixed column sequence. Two sanity checks close the notebook out: `Sample` values must be unique, and the row count must be exactly 512.

**How to run it:** Open in Jupyter with `00_GSE108497_metadata.tsv` and `00_GSE108497_series_matrix.txt` in the same working directory, and run top to bottom. Produces `02_GSE108497_assembled_metadata.csv`.

### `03_Metadata_EDA.ipynb`

**What it does**, in the order the notebook runs:

1. **Loads** the assembled table from step 1 above.
2. **Fixes the NP/PP collision** — donors in the IVF/never-pregnant cohort (identifiable via `_NP_` in the `grp_p_tp` field) had a blank `time_point`, which would otherwise be indistinguishable from a missing value; a lookup-table cross-check first confirms `tp` and `time_point` agree for every pregnant-group row, then `time_point` is explicitly set to `"NP"` for the NP donors.
3. **Builds a donor-level view** (one row per `Donor_id`, via `drop_duplicates`) for demographic tables. Five columns (`Ethnicity`, `Age Group`, `Age`, `Race`, `apl`) sometimes disagree across a donor's own repeated samples — for each, the majority value across that donor's samples is kept; a genuine tie is left as `NA` rather than guessed.
4. **Generates donor-level demographic tables/figures** (age, age group, race, ethnicity, condition) and **sample-level structure tables** (timepoint counts, per-donor sample counts, the longitudinal donor × timepoint heatmap, batch structure) — see `04_Analysis_Results/README.md` for each output individually.
5. **Runs correlation and statistical-test analyses** — a Spearman correlation matrix (chosen over Pearson since it doesn't assume linear relationships) across clinical/pregnancy-outcome variables, plus chi-squared tests (categorical flags) and one-way ANOVA (continuous variables) against Age Group and against Race/Ethnicity.
6. **Saves the corrected full table** (`04_GSE108497_processed_metadata.csv`) — this happens *before* the NP split below, so it still contains all 512 samples.
7. **Splits NP donors out**: the 23 never-pregnant donors (`time_point == "NP"`) are separated into their own file (`GSE108497_NP_stage.csv`, saved into `04_Analysis_Results/`) and excluded from `metadata` going forward.
8. **Re-runs the demographic/correlation breakdown three ways** — "All" (158 donors), "Pregnancy Cohort" (135), and "NP" (23) — as a direct three-way comparison, saving each as its own table/figure.
9. **Saves the final non-NP analytic cohort** (`non_NP_metadata.csv`, into `04_Analysis_Results/`) — this is the file that feeds `01_RawData_&_PCA` next.

**How to run it:** Requires `02_GSE108497_assembled_metadata.csv` (output of the assembly notebook) in the same directory. Creates the `04_Analysis_Results/` folder itself if it doesn't exist. Run top to bottom in Jupyter.

**Environment:** Python (Anaconda-managed local kernel). Packages imported: `pandas`, `matplotlib.pyplot`, `seaborn`, `networkx`, `matplotlib.patches`, `scipy.stats` (`chi2_contingency`, `f_oneway`, `ttest_ind`), `pathlib`. Exact package version numbers aren't pinned anywhere in either notebook — if you want them recorded here, the easiest source is a `pip freeze` (or `conda list`) taken from the environment you actually ran these in; happy to fold that in once you have it.

---

## Results

### `02_GSE108497_assembled_metadata.csv`

512 rows × 29 columns. Header and first two rows:

```text
Sample,GSE,Experimental Strategy,GPL,Condition,Tissue,Gender,Age Group,Age,Race,Ethnicity,grp_p_tp,Sample_name,Donor_id,sle,apl,lac,tp,pe,fd,nnd,pl_insuff,iugr,sga,batch,time_point,ga_at_collection,ga_at_end_of_pregnancy,if_pe_before_or_after_36_weeks
GSM2901826,GSE108497,Expression,GPL10558,Healthy,Whole blood,Female,21-30,25,C,Not Hispanic or Latino,HC_NP_5,HC2013_1,106346,0,0,0,5,0,0,0,0,0,0,3,,,,
GSM2901827,GSE108497,Expression,GPL10558,Healthy,Whole blood,Female,21-30,25,C,Not Hispanic or Latino,HC_NP_5,HC2013-15,134642,0,0,0,5,0,0,0,0,0,0,3,,,,
```

This is the raw merge, before the NP `time_point` fix — note the two rows above (both NP donors, `grp_p_tp = HC_NP_5`) still have a blank `time_point`. **Feeds directly into `03_Metadata_EDA.ipynb`.**

### `04_GSE108497_processed_metadata.csv`

512 rows × 29 columns, same schema as above, but with the NP/PP fix applied — the same two sample rows now read `time_point = NP` instead of blank. This is the full corrected dataset *before* the NP/non-NP split; the two post-split files (`GSE108497_NP_stage.csv`, `non_NP_metadata.csv`) live in `04_Analysis_Results/` and are documented there.

**Used downstream by:** `01_RawData_&_PCA`, `02_DGE_Analysis`, and every later stage that needs sample metadata — either this file directly or its `04_Analysis_Results/non_NP_metadata.csv` derivative, depending on whether NP donors are wanted for that step.
