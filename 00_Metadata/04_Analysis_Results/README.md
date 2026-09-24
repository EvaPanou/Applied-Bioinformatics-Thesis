# 04_Analysis_Results

## Folder Overview

Every figure and summary table produced by `../03_Metadata_EDA.ipynb`, plus the two files that notebook saves when it splits the corrected metadata into the never-pregnant (NP) cohort and the pregnancy+postpartum analytic cohort. This folder is pure output — there is no code here, only what the parent notebook generated, numbered in the order it generates them.

## Folder Structure

This folder doesn't follow the input/code/output numbering convention used elsewhere in the repo (there's no code or input here — everything is an output of `../03_Metadata_EDA.ipynb`). Instead, files are numbered in the order the notebook produces them; a `.csv` and `.png` sharing a number are usually the data and the plot for the same result, but the two numbering sequences drift apart partway through (e.g. `11_conditions_per_batch_prop.csv` pairs with `09_proportional_condition_per_batch.png`, not `11_...png`).

```text
04_Analysis_Results/
├── 01_AgeGroup_counts.csv
├── 01_donor_age_distribution.png
├── 02_donor_AgeGroup_distribution.png
├── 02_samples_per_timepoint.csv
├── 03_samplesVSdonors_per_Condition.png
├── 03_samples_per_condition_timepoint.csv
├── 04_samples_per_patient.csv
├── 04_samples_per_timepoint.png
├── 05_donors_with_5_total_samples.csv
├── 05_samples_per_condition_tp.png
├── 06_donors_with_1_sample.csv
├── 06_longitudinal_graph.png
├── 07_batch_counts.csv
├── 07_samples_per_batch.png
├── 08_proportional_tp_per_batch.png
├── 08_timepoints_per_batch.csv
├── 09_conditions_per_batch.csv
├── 09_proportional_condition_per_batch.png
├── 10_race_and_ethnicity.png
├── 10_timepoints_per_batch_prop.csv
├── 11_conditions_per_batch_prop.csv
├── 11_proportional_condition_per_race.png
├── 12_condition_per_race_ethnicity.png
├── 12_race_ethnicity_counts.csv
├── 13_proportional_condition_per_race_ethnicity.png
├── 13_race_per_condition_proportions.csv
├── 14_correlated_clinical_variables.png
├── 14_race_ethnicity_combined_counts.csv
├── 15_network_clinical_variables.png
├── 15_race_ethnicity_condition_proportions.csv
├── 16_age_by_group.png
├── 16_clinical_variables_correlation_matrix.csv
├── 17_agegroup_by_group.png
├── 17_sorted_strong_correlations.csv
├── 18_age_correlations.csv
├── 18_race_by_group.png
├── 19_ageGroup_clinical_stat_results.csv
├── 19_ethnicity_by_group.png
├── 20_correlation_All.png
├── 20_race_clinical_stat_results.csv
├── 21_correlation_PregnancyCohort.png
├── 21_final_sample_condition_counts.csv
├── 22_correlation_NP.png
├── 22_final_sample_agegroup_counts.csv
├── 23_final_sample_race_counts.csv
├── 24_age_by_group.csv
├── 25_agegroup_by_group.csv
├── 26_race_by_group.csv
├── 27_ethnicity_by_group.csv
├── 28_correlation_matrix_All.csv
├── 29_correlation_matrix_PregnancyCohort.csv
├── 30_correlation_matrix_NP.csv
├── GSE108497_NP_stage.csv
└── non_NP_metadata.csv
```

## File Reference

| File | Description |
|---|---|
| `01_AgeGroup_counts.csv` / `01_donor_age_distribution.png` / `02_donor_AgeGroup_distribution.png` | Donor age distribution, full 158-donor dataset |
| `02_samples_per_timepoint.csv` / `04_samples_per_timepoint.png` | Sample counts per timepoint, full dataset (n=512) |
| `03_samplesVSdonors_per_Condition.png` | Donors vs. samples per Condition, full dataset |
| `03_samples_per_condition_timepoint.csv` / `05_samples_per_condition_tp.png` | Sample counts per Condition × timepoint |
| `04_samples_per_patient.csv` / `06_longitudinal_graph.png` | Per-donor sample presence across timepoints (donor × timepoint matrix) |
| `05_donors_with_5_total_samples.csv` | The 36 donors with a complete 5-timepoint series |
| `06_donors_with_1_sample.csv` | The 28 donors with only 1 sample total |
| `07_batch_counts.csv` / `07_samples_per_batch.png` | Sample counts per processing batch |
| `08_timepoints_per_batch.csv` / `08_proportional_tp_per_batch.png` | Timepoint composition within each batch |
| `09_conditions_per_batch.csv` / `09_proportional_condition_per_batch.png` | Condition composition within each batch |
| `10_race_and_ethnicity.png` | Ethnicity distribution within each Race category, donor-level |
| `10_timepoints_per_batch_prop.csv` | Same as `08_timepoints_per_batch.csv`, as proportions |
| `11_conditions_per_batch_prop.csv` | Proportion of Healthy vs. SLE samples within each of the 4 processing batches |
| `11_proportional_condition_per_race.png` / `13_race_per_condition_proportions.csv` | Proportion of Healthy vs. SLE within each donor Race category |
| `12_condition_per_race_ethnicity.png` | Raw counts of Healthy vs. SLE within each Race–Ethnicity combination |
| `12_race_ethnicity_counts.csv` | Cross-tabulation of Race × Ethnicity (donor-level), with row/column totals |
| `13_proportional_condition_per_race_ethnicity.png` / `15_race_ethnicity_condition_proportions.csv` | Proportion of Healthy vs. SLE within each Race–Ethnicity combination |
| `14_correlated_clinical_variables.png` / `16_clinical_variables_correlation_matrix.csv` | Spearman correlation matrix across 13 clinical/pregnancy-outcome variables, full donor-level dataset (n=158) |
| `14_race_ethnicity_combined_counts.csv` | Donor counts per combined Race–Ethnicity label |
| `15_network_clinical_variables.png` | Network-graph view of the same correlation matrix, edges filtered to \|r\| > 0.3 |
| `17_sorted_strong_correlations.csv` | The same \|r\| > 0.3 pairs as a flat, sorted, de-duplicated list |
| `18_age_correlations.csv` | Spearman correlation of donor Age against each clinical variable |
| `19_ageGroup_clinical_stat_results.csv` | Chi-squared (categorical) / ANOVA (continuous) tests of each clinical variable against Age Group |
| `20_race_clinical_stat_results.csv` | Same tests, against Race–Ethnicity |
| `16_age_by_group.png` / `24_age_by_group.csv` | Age distribution/summary stats: All vs. Pregnancy Cohort vs. NP |
| `17_agegroup_by_group.png` / `25_agegroup_by_group.csv` | Age Group counts: All vs. Pregnancy Cohort vs. NP |
| `18_race_by_group.png` / `26_race_by_group.csv` | Race counts: All vs. Pregnancy Cohort vs. NP |
| `19_ethnicity_by_group.png` / `27_ethnicity_by_group.csv` | Ethnicity counts: All vs. Pregnancy Cohort vs. NP |
| `20_correlation_All.png` / `28_correlation_matrix_All.csv` | Clinical-variable correlation matrix, "All" group only (mathematically identical to item 14/16 — see note below) |
| `21_correlation_PregnancyCohort.png` / `29_correlation_matrix_PregnancyCohort.csv` | Clinical-variable correlation matrix, Pregnancy Cohort only (n=135) |
| `22_correlation_NP.png` / `30_correlation_matrix_NP.csv` | Clinical-variable correlation matrix, NP donors only (n=23) — mostly undefined, see below |
| `21_final_sample_condition_counts.csv` | Final analytic-cohort donor counts by Condition |
| `22_final_sample_agegroup_counts.csv` | Final analytic-cohort donor counts by Age Group |
| `23_final_sample_race_counts.csv` | Final analytic-cohort donor counts by Race |
| `GSE108497_NP_stage.csv` | The 23 NP donors, split out of the corrected metadata, held for specificity checks |
| `non_NP_metadata.csv` | The final analytic cohort — corrected metadata with NP donors removed |

---

## Results

### Donor age & sample structure (Section 5–6 of the notebook)

**`01_AgeGroup_counts.csv`** — donor-level (n=158) age group counts: `21-30`: 81, `31-40`: 71, `11-20`: 3, `41-50`: 3. Donors are overwhelmingly in their 20s–30s.

**`01_donor_age_distribution.png`**:

![Age Distribution of Donors](01_donor_age_distribution.png)

A roughly bell-shaped distribution centered around 29–30, ranging from 19 to 43, with a secondary bump around 36–37.

**`02_donor_AgeGroup_distribution.png`** — the same counts as a bar chart:

![Age Group Distribution of Donors](02_donor_AgeGroup_distribution.png)

**`03_samplesVSdonors_per_Condition.png`** — donors vs. samples, full 158-donor dataset (NP included):

![Comparison of Samples vs. Unique Donors per Condition](03_samplesVSdonors_per_Condition.png)

Healthy: 66 donors / 187 samples (2.8 samples/donor). SLE: 92 donors / 325 samples (3.5 samples/donor) — SLE donors are sampled somewhat more densely on average than Healthy donors. These 187/325 figures match the series-level sample counts stated directly in the raw GEO series matrix (`512 total samples; 325 SLE samples; 187 without SLE`), which is a useful cross-check that the assembled metadata hasn't dropped or duplicated anything. Note this "66 Healthy donors" figure includes all 23 NP donors (who are Healthy by definition) — it's not the same denominator as the 43-Healthy-donor figure used later for the final non-NP analytic cohort (item 21).

**`02_samples_per_timepoint.csv`** / **`04_samples_per_timepoint.png`** — sample counts per timepoint, full dataset (n=512): `<16 weeks`: 122, `16-23 weeks`: 112, `24-31 weeks`: 107, `32-40 weeks`: 90, `PP`: 58, `NP`: 23 (sums to 512).

![Number of Samples per Timepoint](04_samples_per_timepoint.png)

A steady decline in sample count across pregnancy — expected, since not every donor has a sample at every timepoint (some pregnancies end early, some donors weren't sampled at all points).

**`03_samples_per_condition_timepoint.csv`** / **`05_samples_per_condition_tp.png`** — the same breakdown split by Condition:

![Number of Samples per Condition and Timepoint](05_samples_per_condition_tp.png)

SLE consistently outnumbers Healthy at every timepoint (e.g. `<16 weeks`: 84 SLE vs. 38 Healthy), and only Healthy has any `NP` samples (23) — confirming NP is Healthy-only, consistent with the study design. Row sums here reproduce the 187/325 Healthy/SLE totals from item 03 exactly.

**`04_samples_per_patient.csv`**, **`05_donors_with_5_total_samples.csv`**, **`06_donors_with_1_sample.csv`**, **`06_longitudinal_graph.png`** — per-donor sample presence across timepoints. 36 donors have a complete 5-timepoint series; 28 donors have only a single sample (12 of those are the numeric-ID NP donors sampled only once by design, plus several pregnancy-cohort donors like `H59`, `N95`, `T139`, `T75`, `U91` who simply weren't sampled at other timepoints).

![Sample Type per Donor Across Timepoints](06_longitudinal_graph.png)

The heatmap makes the dataset's actual shape visible at a glance: NP donors (top block) only ever have a value in the `NP` column; pregnancy-cohort donors below show a mostly-but-not-fully-populated grid across `<16 weeks` → `PP`, color-coded red (Healthy) or blue (SLE) — this is the visual confirmation of why donor-aware splitting matters everywhere downstream: most donors contribute multiple, but not identical, sets of timepoints.

### Batch structure (Section 6.IV–VI)

**`07_batch_counts.csv`** / **`07_samples_per_batch.png`** — sample counts per batch: batch 4: 325, batch 1: 86, batch 2: 56, batch 3: 45.

![Samples per Batch](07_samples_per_batch.png)

Batch 4 alone accounts for 63% of all samples (325/512) — the batches are far from evenly sized.

**`08_timepoints_per_batch.csv`** / **`08_proportional_tp_per_batch.png`**:

![Proportion of Timepoints per Batch](08_proportional_tp_per_batch.png)

Batch 3 stands out: 26.7% of its samples are NP (12/45), versus 11.6% in batch 1, 1.8% in batch 2, and **0%** in batch 4 — NP samples are concentrated almost entirely in batches 1 and 3.

**`09_conditions_per_batch.csv`** / **`09_proportional_condition_per_batch.png`**:

![Proportion of Conditions per Batch](09_proportional_condition_per_batch.png)

This is the most important finding in this folder: **batch 3 is 97.8% Healthy (44/45 samples)**, while batches 1, 2, and 4 are all majority-SLE (58%, 64%, and 73% SLE respectively). Batch is heavily confounded with Condition in this dataset — batch 3 functions almost as a dedicated "healthy batch." **This matters directly for `02_DGE_Analysis` and `06_Final_Panel_Validation`**: any differential expression or downstream model that doesn't explicitly account for batch risks mistaking a batch effect for a disease effect, particularly for whatever's driving batch 3's near-total Healthy composition. (This is the same concern the `06_Final_Panel_Validation` README references when discussing why batch is adjusted for directly rather than tested separately — this table is the underlying evidence for that decision.)

**`10_timepoints_per_batch_prop.csv`** — the proportional version of item 08, same numbers as percentages.

**`10_race_and_ethnicity.png`** — Ethnicity distribution within each Race category, donor-level:

![Distribution of Ethnicity within Race Categories](10_race_and_ethnicity.png)

Consistent with `12_race_ethnicity_counts.csv` covered below: the `C` (Caucasian) group is the only one with a substantial Hispanic-or-Latino subgroup (~12 of ~80); `AA` and `As` donors are almost entirely Not-Hispanic-or-Latino, and `H` donors are Hispanic-or-Latino by definition of the label.

### Donor-level demographics & cross-tabulations (Section 5–7 continued)

**`12_race_ethnicity_counts.csv`** — Race × Ethnicity cross-tab, donor-level (n=158):

```text
Race_filled,Hispanic or Latino,NA,Not Hispanic or Latino,Total
AA,0,0,23,23
As,0,1,16,17
```

The cohort is predominantly Caucasian (`C`) with smaller African-American (`AA`), Asian (`As`), Hispanic (`H`), and unclassified groups.

**`14_race_ethnicity_combined_counts.csv`** — the same, collapsed into a single combined label per donor (e.g. `"C - Not Hispanic or Latino"`: 80 donors, the largest single group; `"AA - Not Hispanic or Latino"`: 23).

**`11_proportional_condition_per_race.png`**:

![Proportion of Condition within Race Categories](11_proportional_condition_per_race.png)

Disease proportion varies noticeably by Race: African-American donors are ~70% SLE, Hispanic (`H`) donors are ~80% Healthy, and the small `nan`/unclassified-race group is ~82% SLE. **Worth flagging for later stages:** if Race correlates with SLE status this strongly, it's a potential confounder for any downstream analysis that doesn't account for it.

**`13_race_per_condition_proportions.csv`** is the numeric backing for the plot above.

**`12_condition_per_race_ethnicity.png`** (raw counts) and **`13_proportional_condition_per_race_ethnicity.png`** / **`15_race_ethnicity_condition_proportions.csv`** (proportions) break the same relationship down further by combined Race–Ethnicity label — several small subgroups (e.g. `"As - NA"`, `"C - NA"`) are 100% one condition, but those are single- or few-donor cells, not a reliable signal.

![Condition Distribution within Race-Ethnicity Combinations](12_condition_per_race_ethnicity.png)

![Proportion of Condition within Race-Ethnicity Combinations](13_proportional_condition_per_race_ethnicity.png)

### Clinical variable correlations & statistical tests (Section 8–11)

**`14_correlated_clinical_variables.png`** / **`16_clinical_variables_correlation_matrix.csv`** — Spearman correlation across 13 clinical/pregnancy-outcome variables (`sle`, `apl`, `lac`, `tp`, `pe`, `fd`, `nnd`, `pl_insuff`, `iugr`, `sga`, `ga_at_collection`, `ga_at_end_of_pregnancy`, `if_pe_before_or_after_36_weeks`), full donor-level dataset (n=158):

![Spearman Correlation Matrix for Clinical Variables](14_correlated_clinical_variables.png)

The strongest relationship by far is `apl`–`lac` (r=0.83) — antiphospholipid antibodies and lupus anticoagulant are related lab markers, so this is expected rather than a red flag. `sle` correlates negatively with `ga_at_end_of_pregnancy` (r=-0.54) and `tp` (r=-0.42), and `pe` (preeclampsia) also correlates negatively with `ga_at_end_of_pregnancy` (r=-0.55) — consistent with complicated pregnancies ending earlier. `pl_insuff`–`iugr` (r=0.65) reflects that placental insufficiency and growth restriction commonly co-occur clinically.

**`17_sorted_strong_correlations.csv`** — the same matrix flattened to just the \|r\| > 0.3 pairs, sorted:

```text
,,0
apl,lac,0.8345914531829285
ga_at_end_of_pregnancy,if_pe_before_or_after_36_weeks,0.6461987757371916
pl_insuff,iugr,0.6461513141187033
tp,ga_at_collection,0.5109553283018766
```
*(16 pairs total — full list in the file.)*

**`15_network_clinical_variables.png`** — the same \|r\| > 0.3 relationships as a network graph (orange = positive, blue = negative, edge width = |r|):

![Network Graph of Strong Clinical Variable Correlations](15_network_clinical_variables.png)

`sle` sits at a hub connecting to `apl`, `lac`, `tp`, and `ga_at_end_of_pregnancy` — visually confirming it's the most-connected clinical variable in the dataset, consistent with the correlation matrix above.

**`18_age_correlations.csv`** — none of the clinical variables correlate meaningfully with Age (all \|r\| < 0.18); the strongest is `ga_at_end_of_pregnancy` at r=0.098. **Conclusion: Age is not a clinical confounder here.**

**`19_ageGroup_clinical_stat_results.csv`** — only `apl` reaches significance against Age Group (χ²=9.70, p=0.021); `sle` itself is borderline (p=0.052) but not significant.

**`20_race_clinical_stat_results.csv`** — by contrast, Race–Ethnicity shows several significant associations: `sle` (χ²=23.68, p=0.014), `lac` (p=0.031), `pl_insuff` (p<0.001), and `iugr` (p=0.010). **This statistically confirms the Race-vs-Condition imbalance flagged above under item 11** — Race (like batch) is a real potential confounder in this cohort and should be kept in mind for any analysis stage that doesn't explicitly account for it.

### Three-way comparison: All vs. Pregnancy Cohort vs. NP (Section 13)

**`21_final_sample_condition_counts.csv`**, **`22_final_sample_agegroup_counts.csv`**, **`23_final_sample_race_counts.csv`** — donor counts for the final analytic (non-NP) cohort:

```text
Condition,count        Age Group,count      Race,count
SLE,92                 21-30,68              C,84
Healthy,43             31-40,62              AA,21
                        41-50,3               As,14
                        11-20,2               (unclassified),10
                                               Other,4
                                               H,2
```

**`16_age_by_group.png`** / **`24_age_by_group.csv`**:

![Age Distribution: All vs. Pregnancy Cohort vs. NP](16_age_by_group.png)

The NP group skews slightly younger (mean 28.7) than the Pregnancy Cohort (mean 30.6) — a modest difference, not a large one.

**`17_agegroup_by_group.png`** / **`25_agegroup_by_group.csv`**:

![Age Group: All vs. Pregnancy Cohort vs. NP](17_agegroup_by_group.png)

**`18_race_by_group.png`** / **`26_race_by_group.csv`**:

![Race: All vs. Pregnancy Cohort vs. NP](18_race_by_group.png)

**`19_ethnicity_by_group.png`** / **`27_ethnicity_by_group.csv`**:

![Ethnicity: All vs. Pregnancy Cohort vs. NP](19_ethnicity_by_group.png)

Race and Ethnicity proportions look broadly similar between the Pregnancy Cohort and NP group, though NP is a small sample (n=23) so individual-category swings (e.g. `H` donors: 8/23 in NP vs. 2/135 in the Pregnancy Cohort) carry a lot of sampling noise.

**`20_correlation_All.png`** / **`28_correlation_matrix_All.csv`**, **`21_correlation_PregnancyCohort.png`** / **`29_correlation_matrix_PregnancyCohort.csv`**, **`22_correlation_NP.png`** / **`30_correlation_matrix_NP.csv`** — the same clinical-variable correlation matrix, computed independently within each of the three groups:

![Spearman Correlation Matrix — Pregnancy Cohort](21_correlation_PregnancyCohort.png)

![Spearman Correlation Matrix — NP](22_correlation_NP.png)

**Note:** `20_correlation_All.png`/`28_correlation_matrix_All.csv` are mathematically identical to `14_correlated_clinical_variables.png`/`16_clinical_variables_correlation_matrix.csv` above — both are computed on the same full 158-donor `donor_level` table, just at two different points in the notebook (Section 8 and Section 13.VI). Not a bug, just a redundant recomputation worth knowing about if you're looking for where a number comes from.

The NP correlation matrix (n=23) is almost entirely blank — this is expected and called out directly in the notebook's own print statement: NP donors have no pregnancy, so `pe`, `fd`, `nnd`, `pl_insuff`, `iugr`, `sga`, `ga_at_collection`, `ga_at_end_of_pregnancy`, and `if_pe_before_or_after_36_weeks` are either constant (0) or entirely undefined for every one of them, and a correlation is mathematically undefined for a constant column.

### Final split files (Section 12 & 14)

**`GSE108497_NP_stage.csv`** — the 23 never-pregnant donors, held out here rather than discarded. Same 29-column schema as the parent stage's metadata files.

**`non_NP_metadata.csv`** — the final analytic cohort: 489 samples (512 minus the 23 NP-donor rows), same schema. Header and first row:

```text
Sample,GSE,Experimental Strategy,GPL,Condition,Tissue,Gender,Age Group,Age,Race,Ethnicity,grp_p_tp,Sample_name,Donor_id,sle,apl,lac,tp,pe,fd,nnd,pl_insuff,iugr,sga,batch,time_point,ga_at_collection,ga_at_end_of_pregnancy,if_pe_before_or_after_36_weeks
GSM2901849,GSE108497,Expression,GPL10558,Healthy,Whole blood,Female,11-20,20,AA,Not Hispanic or Latino,HC_P_1,JC12(392)6-11w,JC12,0,0,0,1,0,0,0,0,0,0,4,<16 weeks,10,40.3,
```

**This is the file that feeds every later stage of the pipeline** — `01_RawData_&_PCA` and `02_DGE_Analysis` both build on the non-NP cohort defined here, while `GSE108497_NP_stage.csv` is held in reserve as the specificity-check cohort used later in `06_Final_Panel_Validation`.

---

## Summary of confounders identified in this folder

Two independent variables show a statistically meaningful association with Condition (SLE vs. Healthy) in this cohort, both worth carrying forward as caveats into later analysis stages:

- **Batch** — batch 3 is 97.8% Healthy, while the other three batches are all majority-SLE (item 09).
- **Race** — significantly associated with `sle` (χ²=23.68, p=0.014) and several other clinical variables (item 20).

Age, by contrast, shows no meaningful association with Condition or any clinical variable (items 18–19).
