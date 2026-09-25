# 03_Feature_Validation / 02_Results / PERMANOVA

## Folder Overview

The formal, quantitative counterpart to the PCA plots: a PERMANOVA (`vegan::adonis2()`) effect-size test for Condition and for Time, run separately on each of the 3 gene panels. This is the single most information-dense result in this stage — a number for the "how much" question the PCA plots can only show visually. Produced by Step 4 of `../../01_Feature_Validation.R`.

## Folder Structure & File Reference

| File | Description |
|---|---|
| `PERMANOVA_Summary.tsv` | One row per panel × test (6 rows total: 3 panels × {Condition, Time}) |

---

## Results

Full table:

```text
Panel        Stage                Test                                      N    N_genes  R2      F_statistic  P_value
All_Genes    Before (unfiltered)  Condition (donor-level, 1 point/donor)    135  9055     0.0609  8.624        0.001
All_Genes    Before (unfiltered)  Time (donor-blocked, sample-level)        475  9055     0.0659  8.287        0.001
DEG_Union    After (DEG-selected) Condition (donor-level, 1 point/donor)    135  44       0.4112  92.885       0.001
DEG_Union    After (DEG-selected) Time (donor-blocked, sample-level)        475  44       0.0209  2.506        0.001
Final_Panel  After (DEG-selected) Condition (donor-level, 1 point/donor)    135  42       0.4158  94.676       0.001
Final_Panel  After (DEG-selected) Time (donor-blocked, sample-level)        475  42       0.0158  1.883        0.002
```

**Condition R² triples-then-some from before to after**: 6.1% → 41.1% (`DEG_Union`) → 41.6% (`Final_Panel`) — roughly a 7× concentration of Condition-related variance once the expression matrix is restricted to a validated panel. All three reach p=0.001, the smallest p-value obtainable with 999 permutations (`p = 1/(999+1)`) — every panel's Condition effect is as significant as this test can report, before vs. after included, so the *R²* is the number that actually differentiates the panels, not the p-value.

**`Final_Panel` edges out `DEG_Union`** (41.6% vs. 41.1%) despite having 2 fewer genes — consistent with the PCA plots showing no visible loss of separation, and with those 2 genes being the ones Method A's meta-analysis doesn't independently confirm (see `02_DGE_Analysis/02_Results/DE_Genes/README.md`).

**Time R² moves in the opposite direction**, and this is arguably the more informative half of the table: 6.6% (`All_Genes`) → 2.1% (`DEG_Union`) → 1.6% (`Final_Panel`). Gestational timepoint explains noticeably *less* variance once the matrix is restricted to Condition-selected genes — exactly what you'd want to see if the panel is capturing a genuinely Condition-specific signal rather than some general, non-specific difference between samples. This is the PERMANOVA-level echo of `02_DGE_Analysis`'s global interaction F-test (0 of 9,055 genes show significant Condition × Time interaction) and `01_RawData_&_PCA/04_ALASCA_Output`'s multivariate decomposition — three independent analyses, three different statistical tools, all agreeing that Time and Condition are cleanly separable effects in this dataset.

**Used downstream by:** nothing further in the pipeline computationally — this table is the headline quantitative result for the thesis discussion's feature-validation section, referenced from the parent stage's README and the root README's methodology summary.
