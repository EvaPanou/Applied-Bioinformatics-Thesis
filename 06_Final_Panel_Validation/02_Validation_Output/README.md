# 06_Final_Panel_Validation / 02_Validation_Output

## Folder Overview

Everything `../01_Panel_Validation.ipynb` produces — 3 CSVs (one per model, one combined) and 3 diagnostic plots, one per grouping variable tested.

## Folder Structure & File Reference

| File | Description |
|---|---|
| `model1_complication_specificity.csv` | Model 1 results — complication specificity (SLE donors only), 3 genes |
| `model2_pregnancy_specificity.csv` | Model 2 results — pregnancy specificity (Healthy donors only), 3 genes |
| `combined_results.csv` | Both models stacked, 6 rows total — the single table that answers this stage's question |
| `panel_by_batch.png` | Per-gene expression by batch, colored by Condition — all 498 samples |
| `panel_by_complication_group.png` | Per-gene expression by complication group (NC/PE/OC) — SLE donors only |
| `panel_by_pregnancy_status.png` | Per-gene expression by pregnancy status (NP/Pregnant) — Healthy donors only |

---

## Results

### `combined_results.csv` — the headline table

```text
gene,test,likelihood_ratio_statistic,degrees_of_freedom,raw_p_value,fdr_corrected_p_value
BATF2,complication_specificity,2.438,2,0.296,0.443
BATF2,pregnancy_specificity,1.894,1,0.169,0.253
HERC5,complication_specificity,0.443,2,0.801,0.801
HERC5,pregnancy_specificity,0.185,1,0.667,0.667
SPATS2L,complication_specificity,5.445,2,0.066,0.197
SPATS2L,pregnancy_specificity,2.409,1,0.121,0.253
```

**All six FDR-corrected p-values sit above 0.05** — none of the three genes show evidence of tracking complication severity or pregnancy status once SLE status and batch are accounted for. The panel's SLE-vs-Healthy signal isn't attributable to either confound in this cohort.

**Worth reading carefully rather than glossing over**: `SPATS2L`'s complication-specificity row has the smallest raw p-value in the table (0.066) — close enough to 0.05 that the honest description is "no evidence of an effect after FDR correction," not a comfortably settled null the way the other five rows are.

### `model1_complication_specificity.csv` — Model 1 alone

```text
gene,likelihood_ratio_statistic,degrees_of_freedom,raw_p_value,fdr_corrected_p_value
HERC5,0.443,2,0.801,0.801
BATF2,2.438,2,0.296,0.443
SPATS2L,5.445,2,0.066,0.197
```

Whether expression differs across NC/PE/OC within SLE donors, batch-adjusted. Nothing survives FDR correction, so the panel doesn't appear to be tracking complication severity rather than SLE itself. `SPATS2L` is the one gene whose raw p-value would have been nominally noteworthy before correction — worth a one-line mention in a thesis discussion as a direction to watch in a larger cohort, not as a finding on its own.

### `model2_pregnancy_specificity.csv` — Model 2 alone

```text
gene,likelihood_ratio_statistic,degrees_of_freedom,raw_p_value,fdr_corrected_p_value
HERC5,0.185,1,0.667,0.667
BATF2,1.894,1,0.169,0.253
SPATS2L,2.409,1,0.121,0.253
```

Whether expression differs between Healthy-Pregnant and Healthy-Never-Pregnant donors, batch-adjusted. All three raw p-values sit comfortably above 0.05, with no borderline case this time — `BATF2`'s and `SPATS2L`'s FDR values are tied at 0.253 simply because BH-FDR ties the two smallest p-values when they're close together, not because of anything special about either gene. This is the cleaner of the two specificity checks: nothing here suggests the panel is picking up pregnancy status rather than SLE.

### `panel_by_batch.png`

![Final 3-gene panel, expression by batch, colored by Condition](panel_by_batch.png)

This plot is what makes the batch-adjustment rationale concrete rather than hypothetical. Batch 3 is visibly almost all blue (Healthy) across all three genes — consistent with the rest of the repo's own finding that batch 3 is 44 of 45 Healthy samples, with only a single SLE sample in it. That's a real, visible batch/condition imbalance, exactly the kind that could let a processing artifact masquerade as a disease signal if left unadjusted. The fact that both specificity models still come out non-significant *after* adjusting for this confound is more reassuring than a null result would have been from a dataset where batch and Condition were already independent of each other.

### `panel_by_complication_group.png`

![Final 3-gene panel, SLE donors only, by complication group](panel_by_complication_group.png)

The direct visual counterpart to Model 1. Same story as the p-values: no gene shows a clean separation across NC/PE/OC. `SPATS2L`'s boxes shift slightly upward from NC (median ≈7.5) through PE (≈8.1) to OC (≈8.2) — consistent with it being the one borderline result in the table above — but the spread within each group is large enough relative to that shift that it isn't a visually obvious effect either. `HERC5` and `BATF2` show essentially flat medians across all three groups.

### `panel_by_pregnancy_status.png`

![Final 3-gene panel, Healthy donors only, by pregnancy status](panel_by_pregnancy_status.png)

The direct visual counterpart to Model 2. Medians sit close between NP and Pregnant for all three genes, with a somewhat wider spread and a handful of high outliers on the Pregnant side — nothing that reads as a systematic shift.

---

## Are the three genes expected to move together across complication groups?

Not automatically — if the three genes were independent, unrelated measurements, there would be no particular reason for all three to shift the same direction from NC through PE to OC. Seeing them agree is not a given.

It is not surprising for this specific panel, though. `HERC5`, `BATF2`, and `SPATS2L` were not selected independently — they were selected together specifically because they behave similarly across SLE vs. Healthy, and (per `05_ML_Pipeline`'s own methodology) `HERC5` and `BATF2` in particular are general type I interferon markers rather than narrowly SLE-specific ones. Genes on the same biological axis tend to move together: if any shared factor nudges expression up slightly in PE/OC donors relative to NC — not necessarily anything complication-specific — all three genes would plausibly drift the same direction together, simply because they are reading out the same underlying signal.

A concrete candidate for that shared factor comes directly from the source paper: Hong et al. report that baseline disease activity (SLEDAI) differs by complication group — SLE-PE averaged 4.43, SLE-OC 2.86, SLE-NC 2.35 (their own P=0.026). Disease activity tracks the type I interferon response, which is exactly what these three genes are reading. So a small, same-direction shift across all three genes — with none of them individually clearing significance — is consistent with "PE/OC donors happened to have slightly higher baseline disease activity in this cohort," rather than with any complication-specific biology.

This is a useful framing for a thesis discussion section: the consistent direction across genes doesn't strengthen the case that complication status itself is driving expression (the FDR-corrected test already says it isn't), but it does suggest the borderline `SPATS2L` result is more likely a faint echo of general disease activity than noise in one gene alone — worth a sentence in the discussion, not a finding, and not something that changes any result in this stage.
