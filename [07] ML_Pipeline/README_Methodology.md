The goal of this comparison of ML algorithms is to see which method results in the most reliable gene signature for the differentiation of Lupus vs Healthy

# ML PIPELINE METHODOLOGY

This repository implements a donor-aware machine-learning pipeline for
building a diagnostic 6‑gene signature for systemic lupus erythematosus
(SLE) from whole-blood transcriptomic data, starting from an upstream
candidate gene set and ending with evaluation on a sealed held-out donor
set.[1–4]

Key implementation details, including donor-aware cross-validation,
nested model selection, bootstrap stability, Pearson correlation pruning,
and permutation-based significance testing, follow current guidance for
reproducible clinical prediction models.[1,3]

---

## Study design and donor-aware splits

We build a **sample-level diagnostic classifier** (SLE vs healthy) from a
candidate gene set obtained upstream (e.g. limma `duplicateCorrelation`
on repeated visits).[4] The cohort comprises 158 donors (92 SLE, 66
healthy), 507 samples across up to 5 timepoints, and 212 input genes.

To avoid leakage from repeated measures, all splitting is done at the
**donor level**:

- An 80/20 donor split defines **development** vs **held-out** sets
  (`01_split_summary.json`).
- All cross-validation uses `StratifiedGroupKFold` with donor IDs passed
  via `groups=` so no donor contributes samples to both train and test
  folds.[5]

In the recorded run, this produced 404 development samples from 126
donors and 103 held-out samples from 32 donors
(see `01_split_summary.json` for exact counts).

---

## Feature-selection and classifier benchmark

On the **development set only**, we benchmark multiple feature selectors
and classifiers in nested donor-aware cross-validation:[3,6]

- Feature selectors: LASSO, elastic net, random forest importance,
  mRMR, SVM‑RFE, and Boruta.
- Classifiers: logistic regression, Naïve Bayes, SVM, k‑nearest
  neighbours (kNN), random forest (RF), and XGBoost.

For each selector–classifier pair, an inner donor-stratified CV tunes
hyperparameters, and an outer donor-stratified CV estimates performance.
We record mean and standard deviation of outer-fold AUC for every
combination in `02_benchmark_summary.csv`.

The benchmark is visualised via:

- `06_benchmark_barplot.png` – mean AUC per selector–classifier pair.  
- `06_benchmark_boxplot.png` – distribution of outer-fold AUCs.  
- `06_benchmark_heatmap.png` – heatmap of mean AUC across combinations.

All combinations achieve high AUC (≈0.86–0.90), with Boruta + RF ranked
among the best-performing pairs by mean outer-fold AUC.

---

## Final panel construction

### Stability-based filtering

Using the **Boruta + RF** winner as the base model, we perform a
bootstrap stability analysis on the development cohort.[7] We generate
400 donor-level bootstrap resamples, re-run the Boruta + RF pipeline on
each, and record how often each gene is selected. The resulting
bootstrap frequencies are stored in `03_full_stability.csv` and
summarised in `03_full_stability_top20.png`.

Genes with bootstrap selection frequency ≥0.80 define a set of 29
high-stability candidates. These genes, along with their global ranking
under Boruta + RF and an indicator of whether they enter the final
panel, are listed in `03_high_freq_ranked_genes.csv` and
`03_ranked_stable_genes_per_pearson.csv`.

### Elbow analysis and Pearson pruning

We order the 29 stable genes by their development-set importance under
Boruta + RF and probe panel size \(k\) with an inner donor-aware
AUC‑vs‑\(k\) curve (top \(k = 2,\dots,30\) genes). Mean inner‑CV AUC
values are stored in `03_elbow_curve.csv` and rendered as
`03_elbow_curve.png`, which shows a high plateau (≈0.88–0.89) with only
modest variation across \(k\).

To reduce redundancy, we compute the full Pearson correlation matrix for
the 29 stable genes (`03_final_panel_pearson_full_matrix.csv`) and
collect high-correlation pairs (|r| ≥ 0.90) in
`03_final_panel_high_pairs.csv`. Starting from the top-ranked gene, we
iteratively drop later-ranked genes that are too strongly correlated
with genes already kept, until no pair exceeds |r| = 0.90.

Correlations for the final 6‑gene subset are stored in
`03_final_panel_correlation_matrix.csv` and summarised visually in
`03_final_panel_correlation_heatmap.png`.

This yields a 6‑gene panel:

- IFI44L, USP18, EIF2AK2, PARP14, GALM, BATF2.

All six genes are positively correlated (r > 0.75) but by design no pair
has |r| ≥ 0.90. Metadata for the final panel, including the winning
selector/classifier, full ranked list, stability thresholds, and
correlation thresholds, is stored in `03_final_panel_meta.json`.

---

## Held-out evaluation and interpretation

The final random forest classifier (Boruta + RF) is retrained on the
development set using only the 6‑gene panel, with hyperparameters
tuned by inner donor-aware CV (`04_holdout_metrics.json`). We then
evaluate the trained model once on the sealed held-out set
(70 SLE, 33 healthy), saving all predictions to `04_predictions.csv`.

From these, we compute AUC, sensitivity, specificity, accuracy,
F1‑score, and bootstrap CIs, as summarised in `04_holdout_metrics.json`.
ROC coordinates are stored in `04_roc_data.csv` and plotted in
`04_roc_curve.png`.

The confusion matrix and predicted-probability distributions are
visualised in `04_confusion_matrix.png` and
`04_predicted_probability_violin.png`.

To test whether the observed dev-set AUC could arise by chance, we run a
donor-aware permutation test on the development set, following Ojala &
Garriga and Phipson & Smyth.[8,9] The null AUC distribution is stored in
`04_permutation_null.csv`, and summary statistics in
`04_permutation_meta.json`, visualised as
`04_permutation_histogram.png`.

As a transparent comparator recommended by TRIPOD+AI, we fit an
L2-penalised logistic regression on the same 6‑gene panel, tuned by
nested CV.[1] Its held-out AUC and CI are stored in
`04_lr_comparison.json`, and its ROC is overlaid with the RF ROC in
`04_roc_curve_with_lr.png`.

Finally, we compute SHAP values for the random forest on the held-out
set. Per-sample SHAP values are saved in `04_shap_summary.csv`, and mean
absolute SHAP values per gene in `04_sghap_summary.csv`. These are
visualised in `04_shap_barplot.png` and `04_shap_beeswarm.png`.

---

## References

1. Collins GS, Moons KGM, Reitsma JB, et al. TRIPOD+AI: reporting guidelines for studies developing and validating prediction models using artificial intelligence and machine learning. *BMJ*. 2024.  
2. TRIPOD Steering Group. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD) statement. *Ann Intern Med*. 2015.  
3. Cawley GC, Talbot NLC. On over-fitting in model selection and subsequent selection bias in performance evaluation. *JMLR*. 2010;11:2079–2107.  
4. Yang X et al. Longitudinal limma and machine-learning diagnostic signatures in blood transcriptomics. *Front Endocrinol*. 2025.  
5. Saeb S, Lonini L, Jayaraman A, et al. The need to report trial design and statistical methods in biomedical machine learning. *GigaScience*. 2017;6(5):gix020.  
6. Ambroise C, McLachlan GJ. Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS*. 2002;99(10):6562–6566.  
7. Meinshausen N, Bühlmann P. Stability selection. *JRSS B*. 2010;72(4):417–473.  
8. Ojala M, Garriga GC. Permutation tests for studying classifier performance. *JMLR*. 2010;11:1833–1863.  
9. Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Stat Appl Genet Mol Biol*. 2010;9(1):Article39.  
