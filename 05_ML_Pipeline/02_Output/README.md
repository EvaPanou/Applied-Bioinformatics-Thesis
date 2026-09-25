# 05_ML_Pipeline / 02_Output

## Folder Overview

Everything `../01_Code/sle_pipeline.py` produces from a single run — 32 files across 4 sections plus the run log. All real values below are read directly from these files (JSON contents, CSV rows, the pipeline log's own timestamps), not inferred from the code.

## Folder Structure & File Reference

| File | Section | Description |
|---|---|---|
| `01_split_summary.json` | 1 | Donor/sample counts per side of the development/held-out split |
| `02_benchmark_results.csv` | 2 | Every outer-fold AUC — 180 rows (5 folds × 6 FS methods × 6 classifiers) |
| `02_benchmark_summary.csv` | 2 | Mean/SD/count per FS+classifier combination — 36 rows |
| `02_benchmark_boxplot.png` / `02_benchmark_barplot.png` / `02_benchmark_heatmap.png` / `02_benchmark_leaderboard.png` | 2 | Benchmark visualizations |
| `03_final_panel.csv` | 3 | The final 3-gene panel with rank and bootstrap frequency |
| `03_full_stability.csv` | 3 | Bootstrap selection frequency for the 24 genes Boruta confirmed on the full development set |
| `03_high_freq_ranked_genes.csv` | 3 | The 15 genes that passed the ≥80% stability threshold, with FS rank and final-panel membership |
| `03_ranked_stable_genes_pre_pearson.csv` | 3 | Same 15 stable genes, as a plain list, before Pearson pruning |
| `03_final_panel_meta.json` | 3 | Winner, full 44-gene ranking, final panel, thresholds used |
| `03_final_panel_pearson_full_matrix.csv` | 3 | Full 15×15 Pearson correlation matrix over the pre-pruning stable genes |
| `03_final_panel_correlation_matrix.csv` | 3 | 3×3 correlation matrix, final panel only |
| `03_final_panel_pearson_high_pairs.csv` | 3 | Every gene pair at or above the pruning threshold — 83 pairs |
| `03_panel_stability_barplot.png` / `03_full_stability_top20.png` / `03_final_panel_correlation_heatmap.png` | 3 | Final-panel visualizations |
| `04_holdout_metrics.json` | 4 | All held-out performance metrics |
| `04_predictions.csv` | 4 | Per-sample true label, prediction, probability — 96 rows |
| `04_roc_data.csv` | 4 | ROC curve coordinates — 31 points |
| `04_permutation_null.csv` / `04_permutation_meta.json` | 4 | Permutation-test null distribution (200 values) and summary |
| `04_lr_comparator.json` | 4 | Logistic-regression baseline results |
| `04_shap_values.csv` / `04_shap_summary.csv` | 4 | Per-sample (96×3) and mean-absolute (3-row) SHAP values |
| `04_roc_curve_with_lr.png` / `04_confusion_matrix.png` / `04_permutation_histogram.png` / `04_predicted_probability_violin.png` / `04_shap_barplot.png` / `04_shap_beeswarm.png` | 4 | Held-out evaluation visualizations |
| `03_pipeline_log.txt` | Throughout | Full run log — real timestamps used for the runtime breakdown below |

---

## Results

### Section 1 — Split

**`01_split_summary.json`**:

```json
{
  "n_samples_dev": 379, "n_samples_hold": 96,
  "n_donors_dev": 108, "n_donors_hold": 27,
  "n_sle_dev": 250, "n_sle_hold": 63,
  "n_gene_columns": 44
}
```

108 development donors (379 samples), 27 held-out donors (96 samples) — a clean 80/20 split at the donor level. The pipeline log confirms the same numbers plus the donor-level class balance before splitting: **135 donors total, 92 SLE / 43 Healthy** — note this is the *donor*-level count, distinct from the 475 *sample*-level count reported everywhere upstream (162 Healthy samples / 313 SLE samples, per the log's `Sample class counts` line) — donors contribute unequal numbers of repeated samples, which is exactly why every split and CV fold in this pipeline operates at donor level rather than sample level.

### Section 2 — Outer benchmark

**`02_benchmark_summary.csv`** — top of the table, sorted by mean AUC:

```text
fs_method,classifier,mean,std,count
boruta,svm,0.9168,0.0379,5
rf,svm,0.9128,0.0389,5
elasticnet,svm,0.9104,0.0351,5
```

**`02_benchmark_leaderboard.png`** confirms **Boruta + SVM wins** at mean AUC 0.917 ± 0.038 (n=5 folds):

![Benchmark leaderboard, top 15](02_benchmark_leaderboard.png)

**`02_benchmark_heatmap.png`** — all 36 FS×classifier combinations:

![Benchmark heatmap](02_benchmark_heatmap.png)

SVM is the strongest classifier across nearly every feature selector (rightmost-but-one column, consistently dark), and Boruta the strongest feature selector paired with it — the combination isn't a fluke of one classifier or one selector alone.

**`02_benchmark_boxplot.png`** — full fold-level spread for all 36 combinations:

![Benchmark boxplot](02_benchmark_boxplot.png)

**`02_benchmark_barplot.png`** — mean ± SD as a bar chart, same ranking:

![Benchmark barplot](02_benchmark_barplot.png)

### Section 3 — Final panel

**`03_final_panel_meta.json`** — the winner and panel:

```json
{
  "winner": {"fs_method": "boruta", "classifier": "svm", "mean_auc": 0.9168, "std_auc": 0.0379, "n_folds": 5},
  "panel": ["HERC5", "BATF2", "SPATS2L"],
  "final_stability_threshold": 0.8,
  "final_pearson_threshold": 0.9
}
```

The full `ranking` field lists all 44 genes in Boruta's own full-development-set order; `high_freq_ranked_genes` lists the 15 genes that passed the stability bootstrap.

**The funnel, from the pipeline log directly**: 44 candidate genes → Boruta confirms 24 on the full development set → 200-iteration donor-bootstrap stability keeps 15 at ≥80% confirmation frequency → Pearson pruning (`|r|<0.90`) keeps 3.

**`03_final_panel.csv`**:

```text
gene,rank_in_ranking,bootstrap_selection_freq
HERC5,1,0.835
BATF2,12,0.945
SPATS2L,14,0.875
```

`HERC5` ranks #1 in Boruta's own ordering despite having the *lowest* bootstrap stability of the three (83.5% vs. 94.5% and 87.5%) — rank and stability aren't the same thing, and the panel keeps both columns specifically so that distinction stays visible rather than collapsing into a single score.

**`03_panel_stability_barplot.png`**:

![Final panel stability](03_panel_stability_barplot.png)

**`03_full_stability_top20.png`** — the top of the full 24-gene stability list (not just the final 3):

![Top 20 most stable genes](03_full_stability_top20.png)

`OASL` tops the full stability list at ~99% — higher than any of the 3 genes that made the final panel — but gets pruned for redundancy with panel members, not for instability. This is directly confirmed by `03_final_panel_pearson_high_pairs.csv`, which lists 83 gene pairs at `|r|≥0.90` among the 15 pre-pruning stable genes — the interferon-stimulated gene set is highly intercorrelated, so most of the pruning step's work is redundancy removal, not a stability judgment.

**`03_final_panel_correlation_heatmap.png`** — the final 3 genes' own correlations, all still substantial (0.82–0.89) despite surviving pruning against each other:

![Final panel Pearson correlation heatmap](03_final_panel_correlation_heatmap.png)

```text
,HERC5,BATF2,SPATS2L
HERC5,1.00,0.815,0.892
BATF2,0.815,1.00,0.828
SPATS2L,0.892,0.828,1.00
```

**Worth noting explicitly**: the final 3 genes are still correlated with *each other* at 0.82–0.89 — well above what might intuitively seem "independent." They survived pruning only because the threshold (0.90) is a hard cutoff and none of these three pairs quite crosses it; this is a panel of the *least redundant* genes among a highly co-regulated set, not a panel of mutually independent markers.

### Section 4 — Held-out evaluation

**`04_holdout_metrics.json`** — full results:

```json
{
  "auc": 0.8605, "auc_ci_lower": 0.7735, "auc_ci_upper": 0.9209,
  "sensitivity": 0.8095, "specificity": 0.7879, "accuracy": 0.8021, "f1": 0.8430,
  "tp": 51, "fp": 7, "fn": 12, "tn": 26,
  "n_positive": 63, "n_negative": 33, "n_samples": 96,
  "panel": ["HERC5", "BATF2", "SPATS2L"], "classifier": "svm", "fs_method": "boruta"
}
```

**`04_confusion_matrix.png`**:

![Held-out confusion matrix](04_confusion_matrix.png)

51 true positives, 26 true negatives, 12 false negatives, 7 false positives, out of 96 held-out samples (63 SLE, 33 Healthy) — sensitivity trades off against specificity in the direction typical of an ISG-driven panel: the model catches most SLE cases (81%) but is somewhat less reliable at correctly clearing Healthy donors (79%).

**`04_lr_comparator.json`**:

```json
{"auc": 0.8836, "ci_lo": 0.8032, "ci_hi": 0.9393, "best_params": {"C": 10.0}}
```

**The logistic-regression baseline actually outperforms the winning SVM on held-out AUC** (0.884 vs. 0.861) — reported plainly rather than smoothed over. With only 3 features, a simple linear model is a fair competitor to SVM, and this is exactly the kind of transparent baseline comparison the pipeline was built to surface rather than hide.

**`04_roc_curve_with_lr.png`**:

![Held-out ROC curve with LR comparator](04_roc_curve_with_lr.png)

**`04_permutation_meta.json`**:

```json
{"observed_auc": 0.9233, "p_value": 0.005, "n_permutations": 200, "null_mean": 0.4968}
```

**`04_permutation_histogram.png`** — the observed development-set CV AUC (0.923) sits far outside the null distribution (centered at ~0.497, essentially chance):

![Permutation-test null distribution](04_permutation_histogram.png)

p=0.005 is the smallest p-value obtainable from 200 permutations (`1/(200+1)`) — the observed result is as significant as this test can report.

**`04_predicted_probability_violin.png`** — held-out predicted probabilities, split by true class:

![Held-out predicted probabilities by true class](04_predicted_probability_violin.png)

SLE samples cluster tightly near probability 1.0; Healthy samples spread much more widely across the range, with a real tail sitting above the 0.5 decision line — this is the same asymmetry the confusion matrix shows numerically (more false positives among Healthy than false negatives among SLE, proportionally).

**SHAP interpretation** — `04_shap_summary.csv`:

```text
gene,mean_abs_shap
HERC5,0.1445
BATF2,0.1138
SPATS2L,0.0835
```

**`04_shap_barplot.png`**:

![Top SHAP feature importances](04_shap_barplot.png)

**`04_shap_beeswarm.png`**:

![SHAP beeswarm plot](04_shap_beeswarm.png)

High expression (red) of all three genes pushes the model toward an SLE prediction (positive SHAP), low expression (blue) pushes toward Healthy — consistent with all three genes being up-regulated in SLE throughout `02_DGE_Analysis`. `HERC5` has both the highest mean importance and the widest spread of SHAP values, meaning it's driving the largest and most variable share of individual predictions.

---

## Actual runtime, from the log's own timestamps

The "Known open items" note in the parent stage's README flags the Boruta-based stability check's runtime as untested — the log resolves that directly:

| Section | Start | End | Duration |
|---|---|---|---|
| 1 (split) | 05:55 | 05:55 | instant |
| 2 (outer benchmark, 180 fits) | 05:55 | 06:04 | ~9 min |
| 3 (200-iteration bootstrap Boruta stability) | 06:04 | 09:12 | **~3h 8min** |
| 4 (held-out eval + LR comparator + 200-permutation test) | 09:12 | 14:04 | **~4h 52min** |
| **Total** | 05:55 | 14:04 | **~8h 9min** |

The full run took just over 8 hours, dominated almost entirely by the two 200-iteration resampling steps (bootstrap stability and the permutation test) — the actual outer-fold benchmark itself is fast (~9 minutes). This confirms the concern the top-level README raises about Boruta's stability check specifically, and adds that the permutation test is comparably expensive on this dataset.
