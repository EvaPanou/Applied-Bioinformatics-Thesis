# Code architecture

This repository contains all python script files used in this 
pipeline.The pipeline involves phases from data & packages loading, 
donor-strattified splitting tofeature-selection benchmarking, final 
panel construction, held-out test subset evaluation, and figure generation.

At a high level:

- `sle_pipeline.py` is the workflow controller, the only script that you need to run, 
                    it includes the main pipeline that depends on all the other 
                    scripts and calls on them itself in the appropriate 
                    phase. It also creates all the csv output tables.
--------------------------------------------------------------------------------------------
- `requirements.py` manages package installation, if they do not already exist 
                    and imports them for future use in the rest of the scripts
- `configurations.py` defines all input and output directory paths, necessary 
                      data input columns, parameters, and thresholds.
- `helper_functions.py` defines all the necessary main pipeline functions, the donor 
                        stratification, the train-test splitting and implements 
                        the core ML logic, every feature selection method & 
                        classifier benchmark. It also creates all the json output files.
- `plot_creation.py` defines all plot functions and implements the generation
                      of all PNG figures from the saved results.

All scripts must live in the same code directory for the imports to
work correctly.

---

## Script overview

### `requirements.py`

Purpose: environment/bootstrap script that makes sure the required
Python packages are available and then exposes shared imports used by
the rest of the codebase.

It:

- Tries to import each required package; if missing, runs
  `pip install` in the current Python interpreter (via `sys.executable`).
- Uses `matplotlib.use("Agg")` for safe non-interactive plotting.
- Imports and re-exports:
  - Scientific stack: `numpy`, `pandas`, `matplotlib.pyplot`.
  - ML stack: `RandomForestClassifier`, `RFECV`,
    `LogisticRegression`, `LogisticRegressionCV`,
    `GridSearchCV`, `StratifiedGroupKFold`,
    `StratifiedShuffleSplit`, `GaussianNB`,
    `StandardScaler`, `SVC`, `XGBClassifier`, `BorutaPy`.
  - Metrics: `accuracy_score`, `confusion_matrix`,
    `roc_auc_score`, `roc_curve`.
- Handles SHAP as optional, as it sometimes had issues while installing:
  - Tries `import shap`; if successful sets `HAS_SHAP = True` and
    silences SHAP’s own logging noise.
  - If install/import fails, sets `HAS_SHAP = False` so downstream
    code can skip SHAP plots gracefully.

All other scripts import these objects from `requirements.py` rather
than importing them directly.

---

### `configurations.py`

Purpose: central configuration file for paths, metadata columns, random
seeds, cross-validation settings, feature selectors, classifier grids,
bootstrap settings, and final-panel thresholds.

Key elements:

- Paths:
  - `CODE_DIR` = the directory the script is in
  - `PROJECT_DIR` = the parent directory were the Code folder + the data & output folders are within
  - `INPUT_DIR = PROJECT_DIR / "[01] Input"`
  - `OUTPUT_DIR = PROJECT_DIR / "[03] Output"`
  - `DATA_PATH`, `LOG_PATH`, `FIGURES_DIR` = specific filepaths either for input data files, or output files.
- Metadata columns (some existing in more than one column in the original file):
  - Sample & Donor ID columns: `SAMPLE_ID_COL`, `DONOR_ID_COL`
  - Condition Label Columns: `LABEL_COL` (`Condition`) = in 0/1, `LABEL_TEXT_COL`= in text labels "SLE" & "HEALTHY"
  - Timepoint Columns: `TIMEPOINT_COL` = numerical, `TIMEPOINT_LABEL_COL` = categorical.
  - Collected into `METADATA_COLUMNS` for helper functions.
- Reproducibility and splits:
  - `SEED = 42`
  - `HOLDOUT_FRAC = 0.20` (TRAIN - TEST SPLIT, with 20% of samples held-out for testing in the end of the pipeline)
  - Cross Validation Folds Parameters: `N_OUTER_FOLDS = 5`, `N_INNER_FOLDS = 5`
- Feature-selection settings:
  - `FS_METHODS = ("mrmr", "lasso", "elasticnet", "svm_rfe", "rf", "boruta")`
  - Method-specific grids and constants:
    `MRMR_TOP_K`, `LASSO_C_VALUES`, `ELASTICNET_L1_RATIOS`,
    `RF_IMPORTANCE_TREES`, `SVM_RFE_MIN_FEATURES`,
    `BORUTA_MAX_ITER`.
- Classifier benchmark settings:
  - `CLASSIFIERS = ("svm", "logreg", "rf", "xgboost", "nb", "knn")`
  - Grids: `LR_C_VALUES`, `SVM_C_VALUES`, `SVM_GAMMA_VALUES`,
    `RF_DEPTH_VALUES`, `RF_LEAF_VALUES`,
    `XGB_DEPTH_VALUES`, `XGB_LR_VALUES`, `XGB_N_ESTIMATORS`,
    `NB_VAR_SMOOTHING`, `KNN_NEIGHBORS`.
- Stability and final panel:
  - Fold-level bootstrap: `FOLD_BOOTSTRAP_ITERATIONS`,
    `FOLD_BOOTSTRAP_TREES`, `FOLD_BOOTSTRAP_TOP_PERCENTILE`,
    `FOLD_STABILITY_THRESHOLD`,
    `FOLD_STABLE_FALLBACK_TOP_N`, `FOLD_MIN_STABLE_GENES`.
  - Final Gene signature panel parameters: `N_BOOTSTRAP=200`,
    `FINAL_STABILITY_THRESHOLD = 0.80`,
    `FINAL_PEARSON_THRESHOLD = 0.90`.
- Held-out test set evaluation parameters:
  - `N_PERMUTATIONS = 200`
  - `N_HOLDOUT_BOOTSTRAP = 200`.
- Plot colours:
  - `CB_BLUE` = `#0072B2`, `CB_ORANGE` = `#E69F00`, `CB_GREY` = `#999999`
    (neutral palette re-used in all figures).

You only edit this file when changing dataset paths, column names,
model grids, or thresholds.

----
### `helper_functions.py`

Purpose: foundation of the core analysis. This file contains the reusable
functions for data loading, donor-aware splitting, feature selection,
bootstrap stability, nested CV benchmarking, final signature
construction, model training, held-out evaluation, and JSON/CSV
export.

Key functional blocks:

- Basic setup:
  - `setup_logging()`: sets up a logger that shows the pipeline progress and writes both to
    `pipeline_log.txt` and the terminal console.
  - `set_global_seed()`: fixes Python, NumPy, and hash seeds.
  - `ensure_output_dir()`: creates the output folder if it doesn't already exist.
- Donor-stratified cross-validation:
  - `make_donor_stratified_cv()`: returns a `StratifiedGroupKFold`
    object that keeps all samples from a donor together, so that some do not spill over in the testing set.
  - `assert_no_donor_leakage()`: sanity-check that train/test donors
    are disjoint.
  - `summarise_split()`: logs basic counts for each split.
- Data loading and sealed split:
  - `load_dataset()`: reads the TSV/CSV matrix at `DATA_PATH`, checks
    `METADATA_COLUMNS`, separates metadata from gene columns, and
    returns `df`, `X` (gene expression), `y` (0/1 labels), and donor IDs.
  - `make_sealed_donor_split()`: creates an 80/20 (configurable)
    development vs held-out split at the donor level using
    `StratifiedShuffleSplit`, then logs the split.
- Stability filtering:
  - `run_bootstrap_rf_stability()`: donor-level bootstrap:
    resamples donors, fits a Random Forest, tracks which genes fall in
    the top importance percentile each iteration, and computes
    selection frequencies. Applies thresholds and fallback logic to
    define per-fold stable genes.
- Feature selection:
  - `run_mrmr_selector()`: runs mRMR via the external `mrmr` package
    when available, or falls back to mutual-information relevance plus
    correlation-based redundancy. Returns both selected genes and a
    ranking.
  - Additional selectors (lasso, elastic net, SVM-RFE, RF importance,
    Boruta) are implemented and called from higher-level helpers
    (e.g. `run_selector` and the outer benchmark functions).
- Benchmark and winner selection:
  - Functions to run outer nested CV on the development set:
    multiple selectors × classifiers, donor-aware, recording fold-level AUC.
  - `pick_winner()`: summarises `results_df` and picks the
    `fs_method` + `classifier` combination with the highest mean
    AUC.
- Final gene signature construction:
  - `build_final_signature()`: reruns the winning feature selector on
    the full development set, then calls the elbow routine to choose
    the optimal panel size based on inner-CV AUC.
  - Since the elbow panel size was verified to be extremely small 
    and the elbow curve not optimal, the pipeline is built to create 
    the final gene panel another way
  - The ouput of the feature selector on the full development set 
    is bootstrapped for stability:
  - `finalise_panel(...)` (called from `sle_pipeline.py`): wraps
    final stability, ranking, and Pearson-based pruning to produce the
    final 6-gene panel and associated metadata.
- Final model training:
  - `make_classifier()`: creates an estimator object by name
    (`svm`, `logreg`, `rf`, `xgboost`, `nb`, `knn`).
  - `train_final_model()`: takes the final signature, fits a scaler,
    runs donor-aware inner CV grid-search, and returns the best
    classifier with its hyperparameters and inner-CV AUC.
- Held-out evaluation:
  - `evaluate_on_holdout()`: applies the final model to the sealed
    held-out subset, returning AUC, accuracy, sensitivity, specificity,
    F1, confusion matrix, individual predictions, and probabilities.
- Utilities:
  - `save_json()`: small helper to save structured outputs as JSON.

This script is where all donor-aware modelling and stability logic
actually live; the main pipeline simply orchestrates these functions
in a defined order.

---

### `plot_creation.py`

Purpose: figure-generation module. It contains all plotting functions
for benchmark performance, final panel properties, held-out
evaluation, permutation null, and SHAP-based interpretation.

Internal helpers:

- `get_log()`: returns a progress log file.
- `ensure_figures_dir()`: ensures an `output/figures` subfolder exists.
- `save_figure()`: standardised Matplotlib saving routine (tight
  layout, 300 dpi, logging).

Main plot functions (selection):

- Benchmark:
  - `plot_benchmark_boxplot()`: outer-fold AUC distributions per
    feature-selector + classifier combination.
  - `plot_benchmark_barplot()`: mean outer-fold AUC ± SD per
    combination.
  - `plot_benchmark_heatmap()`: mean AUC heatmap (rows = selectors,
    columns = classifiers).
- Final gene signature panel:
  - `plot_elbow_curve()`: inner-CV AUC vs. panel size, marking
    the selected `elbow_k`.
  - `plot_stability_barplot()`: stability values for genes in the
    final panel.
  - `plot_full_stability_top20()`: top-20 most stable genes from the
    full bootstrap frequency series.
  - `plot_final_panel_correlation_heatmap()`: Pearson correlation
    matrix for the final panel genes.
- Held-out test set evaluation:
  - `plot_roc_curve()`: held-out ROC curve with AUC and optional CI.
  - Additional functions (called from the main pipeline) generate
    confusion matrices, permutation histograms, and violin plots of
    predicted probabilities.
- SHAP:
  - `plot_shap_barplot()`: mean absolute SHAP values per gene.
  - `plot_shap_beeswarm()`: SHAP beeswarm plot for per-sample effects
    (when SHAP is installed and available).

All functions accept an `out_dir` and write PNG files into the
`figures` subdirectory of the main output folder.

---

### `sle_pipeline.py`

Purpose: top-level **workflow controller**. This is the script you run
to execute the full pipeline from raw input file to final tables,
JSON files, and figures.

It:

1. Imports:
   - `requirements` (environment/tools).
   - `configurations` (paths, parameters).
   - `helper_functions` (ML logic).
   - `plot_creation` (visualization).
2. Sets up output folder, logging, and global random seeds.
3. Section 1 – Data loading and donor split:
   - Checks that `DATA_PATH` exists.
   - Calls `helpers.load_dataset(...)`.
   - Calls `helpers.make_sealed_donor_split(...)`.
   - Saves a split summary to `01_split_summary.json`.
4. Section 2 – Outer benchmarking (development set only):
   - Calls `helpers.run_outer_benchmark(...)`.
   - Saves `02_benchmark_results.csv` (fold-level AUCs).
   - Summarises to `02_benchmark_summary.csv`.
   - Produces benchmark boxplot, barplot, and heatmap via
     `plot_creation.py`.
5. Section 3 – Final panel and stability:
   - Calls `helpers.finalise_panel(...)` to build the final panel with
     bootstrap stability and Pearson correlation pruning.
   - Saves:
     - `03_final_panel.csv`
     - `03_full_stability.csv`
     - `03_elbow_curve.csv`
     - `03_final_panel_meta.json`
     - `03_high_freq_ranked_genes.csv`
     - Pearson correlation matrices and high-correlation pairs.
   - Generates elbow, stability, and correlation plots.
6. Section 4 – Held-out test set evaluation:
   - Trains the final model on the development set using the final
     panel (`helpers.train_final_model`).
   - Evaluates on the sealed held-out set
     (`helpers.evaluate_on_holdout`).
   - Computes ROC coordinates, bootstrap AUC CI, and a permutation null
     distribution on the development set.
   - Optionally trains a logistic regression comparator on the same
     panel and evaluates it on the held-out set.
   - Computes SHAP values if SHAP is installed and produces SHAP-based
     plots.
   - Saves:
     - `04_predictions.csv`
     - `04_roc_data.csv`
     - `04_holdout_metrics.json`
     - `04_permutation_null.csv` + `04_permutation_meta.json`
     - Optional `04_lr_comparator.json`
     - `04_shap_values.csv` and `04_shap_summary.csv`.
   - Generates ROC, confusion-matrix, permutation, probability, and
     SHAP plots.

You can run the full pipeline with:

```bash
python sle_pipeline.py
```

from the project’s code directory.

---

## Execution flow (high level)

1. Load packages and configuration.
2. Load labelled donor-level expression matrix.
3. Create an 80/20 donor-level development/held-out split.
4. Run nested donor-aware outer CV benchmark on the development set
   for multiple selectors × classifiers.
5. Select the best-performing selector/classifier combination.
6. Build a stable, low-redundancy final gene panel using bootstrap
   stability and Pearson correlation pruning.
7. Train a final tuned classifier on the development set using only
   the final panel.
8. Evaluate once on the sealed held-out donors.
9. Generate all tables, JSON summaries, and figures.
