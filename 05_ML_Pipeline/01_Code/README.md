# 05_ML_Pipeline / 01_Code

## Folder Overview

Everything `sle_pipeline.py` needs to run — 5 Python files that must stay together (the imports between them are relative to this folder). `sle_pipeline.py` is the only one you actually run; the other four are imported by it (`configurations` for settings, `requirements` for third-party imports, `helper_functions` for the ML logic, `plot_creation` for every figure).

## Folder Structure & File Reference

| File | Role |
|---|---|
| `sle_pipeline.py` | Workflow controller — the only script you run. Coordinates every section in order and writes every output table, JSON file, and log entry. One function, `main()`, called at the bottom of the file. |
| `configurations.py` | All paths, metadata column names, the random seed, CV/split settings, feature-selector and classifier grids, and thresholds. |
| `requirements.py` | Imports and re-exports every third-party package the rest of the code needs, so the other scripts write `from requirements import *` instead of repeating import lists. |
| `helper_functions.py` | All the actual ML logic — 29 functions. |
| `plot_creation.py` | Every figure the pipeline produces — 17 functions (13 save a figure; the rest are small shared helpers). |

---

## `configurations.py`

Derives every path from `Path(__file__).resolve().parent` — the one script in this repo that **doesn't** have a hardcoded-absolute-path problem, since it resolves relative to wherever the repo actually sits. Key settings:

| Setting | Value |
|---|---|
| `SEED` | 42 |
| `HOLDOUT_FRACTION` | 0.20 |
| `N_OUTER_FOLDS` / `N_INNER_FOLDS` / `N_REPEATS` | 5 / 5 / 1 |
| `FS_METHODS` | `mrmr`, `lasso`, `elasticnet`, `svm_rfe`, `rf`, `boruta` |
| `CLASSIFIERS` | `svm`, `logreg`, `rf`, `xgboost`, `nb`, `knn` |
| `N_BOOTSTRAP` / `FINAL_STABILITY_THRESHOLD` | 200 / 0.80 |
| `FINAL_PEARSON_THRESHOLD` | 0.90 |
| `N_PERMUTATIONS` / `N_HOLDOUT_BOOTSTRAP` | 200 / 200 |

Also defines the colorblind-safe palette used throughout the plots (`CB_BLUE = "#0072B2"` for SLE/primary, `CB_ORANGE = "#E69F00"` for Healthy/comparator, `CB_GREY = "#999999"` for null distributions).

## `requirements.py`

Imports: standard library (`json`, `logging`, `os`, `random`, `collections.Counter`, `dataclasses.dataclass`, `pathlib.Path`, `typing`); plotting (`matplotlib`, forced to the non-interactive `Agg` backend before `pyplot` is imported, so the pipeline can run headless); `numpy`, `pandas`, `boruta.BorutaPy`, `xgboost.XGBClassifier`, `kneed.KneeLocator`; scikit-learn (`RandomForestClassifier`, `RFECV`, `LogisticRegression`/`LogisticRegressionCV`, the 5 metrics functions used in evaluation, `GridSearchCV`/`StratifiedGroupKFold`/`StratifiedShuffleSplit`, `GaussianNB`, `StandardScaler`, `SVC`). SHAP is imported inside a `try/except`, setting a module-level `HAS_SHAP` flag (deliberately with no leading underscore, so `import *` doesn't hide it from the other scripts) — if SHAP isn't installed, the flag is `False` and the pipeline later skips the SHAP plots gracefully rather than crashing.

## `helper_functions.py` — function reference

**Setup & splitting**

| Function | What it does |
|---|---|
| `setup_logging()` | Configures a logger that writes to both the console and `03_pipeline_log.txt`. |
| `donor_stratified_cv()` | Returns a `StratifiedGroupKFold` splitter — donor-aware CV, reused for both outer and inner loops. |
| `assert_no_donor_leakage()` | Hard-fails with a clear error if any donor ID appears in both a train and test split — called after every split, not assumed safe. |
| `summarise_split()` | Logs sample/donor/class counts for a named split. |
| `load_dataset()` | Reads `00_ML_Input_File.tsv`, splits it into metadata, the gene matrix (`X`), labels (`y`), donor IDs, and gene column names. |
| `make_sealed_donor_split()` | Builds the development/held-out split at donor level via `StratifiedShuffleSplit` on one row per donor; raises `RuntimeError` if any donor overlap is detected between the two sets. |

**Feature selection**

| Function | What it does |
|---|---|
| `run_mrmr_selector()` | Tries the external `mrmr` package first (ranks all genes, knee-detects the cutoff on the relevance curve); falls back to a hand-written mutual-information-relevance-minus-correlation-redundancy loop if the package isn't available. |
| `run_lasso_selector()` | L1-regularized logistic regression; genes with non-zero coefficients are selected. |
| `run_elasticnet_selector()` | Elastic Net logistic regression (L1 + L2), same non-zero-coefficient selection logic. |
| `run_svm_rfe_selector()` | Recursive feature elimination with a linear SVM, donor-aware inner CV (`RFECV`). |
| `run_rf_selector()` | Random Forest importances, knee-detected cutoff (via `_knee_cutoff()`). |
| `run_boruta_selector()` | `BorutaPy` wrapped around a class-weighted Random Forest — compares real features against shuffled "shadow" features and keeps everything Boruta confirms relevant. This is the method that ultimately wins the outer benchmark. |
| `run_selector()` | Dispatches to one of the six functions above by name string. |
| `_scale()` / `_knee_cutoff()` | Shared helpers — standard scaling, and `kneed.KneeLocator`-based elbow detection used by mRMR and RF importance. |

**Benchmarking**

| Function | What it does |
|---|---|
| `_classifier_grid()` | Builds the six classifiers with their hyperparameter search grids, applying class-imbalance weighting to every one that supports it (including XGBoost's `scale_pos_weight`, added in this pass — see the top-level README's "Fixes applied"). |
| `_benchmark_classifiers_one_panel()` | For one outer fold and one feature-selector's chosen panel: tunes and evaluates all six classifiers via inner CV, returns one AUC per classifier. |
| `run_outer_benchmark()` | The full Section 2 loop — 5 outer folds × 6 feature selectors × 6 classifiers, collecting one AUC per combination per fold. |
| `summarise_benchmark()` | Collapses fold-level results into mean/SD/count per FS+classifier pair. |
| `pick_winner()` | Selects the FS+classifier combination with the highest mean outer-fold AUC. |

**Final panel**

| Function | What it does |
|---|---|
| `build_final_signature()` | Reruns the winning feature selector on the *full* development set (not just one fold), keeping its real selected genes. |
| `run_bootstrap_boruta_stability()` | Resamples **donors** (not rows) with replacement, reruns Boruta on each resample, and keeps genes confirmed at or above the stability threshold across iterations — with a fallback to the top-N most-frequent genes if too few pass. This is the method described in the top-level README's Section 3, grounded in Kursa (2014). |
| `prune_correlated_genes()` | Walks a ranked gene list in order, dropping any gene whose correlation with an already-kept gene exceeds the Pearson threshold — returns the pruned list and the full correlation matrix (computed once, reused by the caller rather than recomputed). |
| `finalise_panel()` | Orchestrates all of the above: pick winner → rebuild signature on full dev set → bootstrap stability (restricted to the winner's own genes) → Pearson pruning → final panel table. |

**Model training & held-out evaluation**

| Function | What it does |
|---|---|
| `make_classifier()` | Instantiates one classifier by name with given hyperparameters. |
| `train_final_model()` | Fits the winning classifier on the full development set, using only the final panel's genes, with the winning hyperparameters from the outer benchmark. |
| `evaluate_on_holdout()` | Applies the trained model to the held-out set once — classification threshold fixed at probability ≥0.5 — returns AUC, accuracy, sensitivity, specificity, F1, and the raw confusion-matrix counts. |
| `save_json()` | Small shared helper — writes a dict to a JSON file with consistent formatting. |
| `bootstrap_auc_ci()` | Bootstraps the held-out predictions (200 resamples) to get a confidence interval on AUC. |
| `permutation()` | Donor-stratified label permutation on the development set (200 permutations) — builds a null AUC distribution to compute a p-value for the observed result. |

---

## `plot_creation.py` — function reference

| Function | Produces |
|---|---|
| `get_log()` / `ensure_figures_dir()` / `save_figure()` | Shared helpers — logger fallback, output-directory creation, consistent figure saving. |
| `plot_benchmark_boxplot()` | `02_benchmark_boxplot.png` |
| `plot_benchmark_barplot()` | `02_benchmark_barplot.png` |
| `plot_benchmark_heatmap()` | `02_benchmark_heatmap.png` |
| `plot_benchmark_leaderboard()` | `02_benchmark_leaderboard.png` |
| `plot_stability_barplot()` | `03_panel_stability_barplot.png` (had a real column-detection bug fixed this pass — see the top-level README) |
| `plot_full_stability_top20()` | `03_full_stability_top20.png` |
| `plot_final_panel_correlation_heatmap()` | `03_final_panel_correlation_heatmap.png` |
| `plot_roc_curve()` / `plot_roc_curve_with_lr()` | `04_roc_curve.png` or `04_roc_curve_with_lr.png` (the `_with_lr` variant is used whenever the LR comparator actually ran) |
| `plot_confusion_matrices()` | `04_confusion_matrix.png` |
| `plot_permutation_histogram()` | `04_permutation_histogram.png` |
| `plot_predicted_probability_violin()` | `04_predicted_probability_violin.png` |
| `plot_shap_barplot()` | `04_shap_barplot.png` |
| `plot_shap_beeswarm()` | `04_shap_beeswarm.png` |

See [`../02_Output/README.md`](../02_Output/README.md) for every one of these plots with the actual result shown.
