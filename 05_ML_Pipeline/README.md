# 05_ML_Pipeline

## Overview

This stage runs the full donor-aware machine-learning pipeline: starting from the 44-gene ML-ready input table (`00_ML_Input_File.tsv`, produced by `04_ML_Prerequisites`), it benchmarks six feature-selection methods against six classifiers, builds a final compact gene-signature panel, and evaluates it once on a sealed held-out set of donors.

This is a from-scratch rerun of the pipeline against the corrected, 44-gene `Method_All_Union` panel — not a continuation of any earlier run. It also carries a full line-by-line review and fix pass across every script (see "Fixes applied in this pass" below).

**A note on this folder's documentation history:** two files named `README.md` were uploaded for this stage in the same batch. Only one survived (the one this document is built from) — the platform keeps just one file when two uploads share an exact filename. The other described an older `[01] Input/ [02] Code/ [03] Output/ [04] Flowchart/` folder-naming scheme with PDF/PNG flowchart diagrams, which doesn't match the current, verified folder structure below (no `Flowchart/` subfolder exists among the files actually confirmed here). Given this stage's own README explicitly describes itself as a from-scratch rerun with "a full line-by-line review and fix pass," the flowchart README is almost certainly a stale draft from before that rebuild — the same pattern as the stale root README and the duplicate ALASCA README caught earlier in this repo. If a `[04] Flowchart/` folder with those diagrams does still exist in the real repo, it wasn't part of what was uploaded here and isn't covered by this documentation pass.

## Pipeline context

```text
04_ML_Prerequisites
      |
      v
05_ML_Pipeline   <- this stage
      |
      v
06_Final_Panel_Validation
```

## Folder Structure & File Reference

| File / folder | Description |
|---|---|
| `00_ML_Input_File.tsv` | Foundational input — 475 samples × 50 columns (6 metadata + 44 genes), from `04_ML_Prerequisites`. |
| `00_requirements.txt` | Pinned dependency list — `pip install -r 00_requirements.txt` from within this folder. |
| [`01_Code/`](./01_Code/README.md) | Everything `sle_pipeline.py` needs to run — 5 Python files. |
| [`02_Output/`](./02_Output/README.md) | Everything running `01_Code` produces — tables, JSON, figures, log (23 files). |

Figures save directly into `02_Output` alongside the tables and JSON files — there is no separate `figures` subfolder.

---

## Input

**`00_ML_Input_File.tsv`** — 475 samples × 50 columns: `SampleID`, `DonorID`, `Condition` (0/1), `Condition_label`, `Timepoint` (ordinal 1–5), `Time_label`, then the 44 `Method_All_Union` gene columns (log2 expression), unmodified. Produced by `04_ML_Prerequisites`.

---

## Code architecture

| Script | Role |
|---|---|
| `sle_pipeline.py` | Workflow controller — the only script you run. Coordinates every stage in order and writes every output table, JSON file, and log entry. |
| `configurations.py` | All paths, metadata column names, the random seed, CV/split settings, feature-selector and classifier grids, and thresholds. Edit this file to change any parameter. |
| `requirements.py` | Imports and re-exports every third-party package the rest of the code needs (numpy, pandas, scikit-learn, XGBoost, Boruta, kneed, SHAP), so the other scripts just write `from requirements import *`. |
| `helper_functions.py` | All the actual ML logic — donor-stratified splitting, the six feature selectors, bootstrap stability filtering, the outer benchmark, final-panel construction, Pearson pruning, final model training, held-out evaluation, permutation testing. 29 functions. |
| `plot_creation.py` | Every figure the pipeline produces, one function per plot — 17 functions (13 produce a saved figure; the rest are small shared helpers). |

All five files must live together in `01_Code` — the imports between them depend on that. See [`01_Code/README.md`](./01_Code/README.md) for a function-by-function breakdown of the two logic-heavy files.

---

## Methodology

### Section 1 — Input & donor split

Loads `00_ML_Input_File.tsv`, then splits at the **donor** level (not sample level) into a development cohort and a sealed held-out cohort (`HOLDOUT_FRACTION` in `configurations.py`, default 20%, via `StratifiedShuffleSplit` on one row per donor). A donor's repeated timepoint samples always stay together in one side of the split — checked explicitly with a hard `RuntimeError` if any donor overlap is found, not assumed. The held-out cohort is not touched again until Section 4.

### Section 2 — Outer benchmark (development set only)

Donor-aware 5-fold outer cross-validation (`StratifiedGroupKFold`). Within each outer training fold:
1. Each of six feature-selection methods (mRMR, LASSO, Elastic Net, SVM-RFE, RF importance, Boruta) runs directly on the fold's full 44-gene set and decides its own panel size independently — mRMR and RF importance use knee detection (`kneed.KneeLocator`) on their own score curves; LASSO/Elastic Net/SVM-RFE/Boruta are naturally self-determining.
2. Each of six classifiers (Logistic Regression, Naive Bayes, SVM, kNN, Random Forest, XGBoost) is tuned via donor-aware inner 5-fold CV (`GridSearchCV` with a `StratifiedGroupKFold` splitter) and evaluated on that outer fold's held-out donors.
3. Every outer-fold AUC is recorded; the combination with the highest mean AUC across folds wins.

(An earlier version of this stage also ran a bootstrap Random Forest stability filter here, narrowing the candidate pool before feature selection. Removed after reviewing real results on this 44-gene panel, where the filter's threshold rarely passed naturally and a fallback mechanism was doing most of the work instead of the intended stability logic — see "Fixes applied in this pass".)

### Section 3 — Final panel

The winning feature selector is rerun on the full development set, keeping its own real selected genes (not a separate re-sizing step). Those genes go through a second, stricter stability pass (`run_bootstrap_boruta_stability()` — 200 resamples, ≥80% threshold) restricted to that selected set. This resamples **donors** (not rows) with replacement and reruns Boruta itself on each resample, keeping genes it confirms consistently — not a generic RF-importance percentile cutoff (see "Fixes applied in this pass"). Then Pearson correlation pruning (`prune_correlated_genes()`, `|r| >= 0.90`) walks the genes in ranked order and drops any gene too correlated with an already-kept one, producing the final compact panel.

### Section 4 — Held-out evaluation

The winning classifier is trained once on the full development set using only the final panel (`train_final_model()`), then applied once to the sealed held-out donors (`evaluate_on_holdout()` — classification threshold fixed at probability ≥0.5). Reports ROC AUC, sensitivity, specificity, accuracy, F1, and a confusion matrix, plus a bootstrap confidence interval on held-out AUC (`bootstrap_auc_ci()`, 200 resamples) and a donor-stratified permutation test on the development set (`permutation()`, 200 permutations, labels shuffled at the donor level) to confirm the result isn't attributable to chance. A logistic-regression comparator is trained on the same panel as a transparent baseline **only if the winning classifier isn't already logreg** (it wasn't, here — the winner was SVM). SHAP values are computed on the held-out donors: `TreeExplainer` if the winning classifier is RF or XGBoost, otherwise a `KernelExplainer` with a 10-point k-means background (falls back gracefully — wrapped in its own `try/except` — if SHAP isn't installed or errors).

---

## Key parameters (`configurations.py`)

| Parameter | Value | What it controls |
|---|---|---|
| `SEED` | 42 | Every random operation in the pipeline |
| `HOLDOUT_FRACTION` | 0.20 | Donor-level development/held-out split |
| `N_OUTER_FOLDS` / `N_INNER_FOLDS` | 5 / 5 | Nested CV structure |
| `FS_METHODS` | 6 methods | Which feature selectors are benchmarked |
| `CLASSIFIERS` | 6 classifiers | Which classifiers are benchmarked |
| `FOLD_BOOTSTRAP_TREES` / `FOLD_BOOTSTRAP_TOP_PERCENTILE` / `FOLD_STABLE_FALLBACK_TOP_N` / `FOLD_MIN_STABLE_GENES` | 200 / 20 / 20 / 5 | Reused as Section 3's final-panel stability settings (fallback size, minimum stable count) — no longer used for a fold-level pre-filter, see Section 2 |
| `N_BOOTSTRAP` / `FINAL_STABILITY_THRESHOLD` | 200 / 0.80 | Final-panel stability pass (Section 3) |
| `FINAL_PEARSON_THRESHOLD` | 0.90 | Redundancy pruning cutoff |
| `N_PERMUTATIONS` | 200 | Development-set permutation test |
| `N_HOLDOUT_BOOTSTRAP` | 200 | Held-out AUC confidence interval |

---

## Fixes applied in this pass

This pass took the pipeline through a full line-by-line review against the 44-gene input. Highlights: mRMR and RF-importance now use knee detection on their own score curves instead of a fixed top-30 cap that no longer made sense at this scale; the now-redundant elbow-curve re-sizing step was removed from both the outer-fold loop and final-panel construction; the final panel's bootstrap stability and Pearson pruning now correctly run only on the winning method's own selected genes, instead of the full candidate pool; XGBoost now receives class-imbalance weighting (`scale_pos_weight`) like every other classifier; several crash-causing references to functions and fields that no longer exist were caught and fixed; a genuine bug in `plot_stability_barplot`'s column detection (which meant that plot had likely never once generated successfully) was fixed; the fold-level bootstrap stability pre-filter was removed from the outer-fold loop after real results on this data showed its threshold rarely passed naturally; a separate bug was caught in the same area where the final-panel stability step was silently running 400 bootstrap iterations instead of the intended 200 (`N_BOOTSTRAP`), due to a default parameter pointing at the wrong config constant; the final-panel stability check (`run_bootstrap_boruta_stability`, replacing `run_bootstrap_rf_stability`) now resamples donors and reruns Boruta itself each time, keeping genes it confirms consistently, instead of a generic RF-importance top-percentile cutoff — grounded in Kursa (2014), the paper this project's own stability methodology already cites, which found Boruta the most self-consistent of the RF-based gene selection approaches it compared, and motivated by real evidence from a run on this data where the old mechanism passed zero genes and the resulting panel was anchored on a gene (`HERC5`) that was itself one of the least stable by the old measure; and `00_requirements.txt` was fixed twice — first because its header used Python-style `"""` comments, which pip can't parse at all, and then because its package versions were unpinned, which let one install silently jump numpy/pandas/scikit-learn to incompatible newer major versions (now pinned to a verified-compatible set).

## Known open items

- **Runtime of the Boruta-based final-panel stability check was untested at full scale when this pass was written — the actual run's log resolves it.** The full pipeline run took **~8h 9min** total, broken down as: Section 2 (outer benchmark) ~9 min, **Section 3 (the 200-iteration bootstrap Boruta stability check) ~3h 8min**, Section 4 (held-out eval + the 200-permutation test) ~4h 52min. The concern was justified — see [`02_Output/README.md`](./02_Output/README.md) for the full timing breakdown from the log's own timestamps.
- **mRMR's second return value** stays as its own short list rather than the full gene ranking (unlike RF importance), specifically so this pass doesn't silently change what the outer-fold loop's panel-size logic saw for mRMR.
- **Per-fold gene identities are not logged or saved anywhere** — `panels_per_fold` (built during the outer benchmark) is computed but never written to a file or logged by gene name, only by count. If you want to trace exactly which genes each fold/method selected, this needs a small addition (a log line was drafted during review but not yet added).
- **`FOLD_BOOTSTRAP_ITERATIONS`, `FOLD_BOOTSTRAP_TREES`, `FOLD_BOOTSTRAP_TOP_PERCENTILE`, and `FOLD_STABILITY_THRESHOLD`** in `configurations.py` are now unused in practice, left over from the removed fold-level pre-filter — kept in place as documented-but-unused rather than deleted.

---

## How to run in VS Code

**One-time setup (only needed the first time, or after changing packages):**

1. Open the `05_ML_Pipeline` folder in VS Code (`File > Open Folder...`).
2. Open a terminal in VS Code (`Terminal > New Terminal`, or `` Ctrl+` ``).
3. Activate your Anaconda environment if it isn't already active: `conda activate <your-environment-name>`
4. Install the required packages: `pip install -r 00_requirements.txt`
5. In the bottom-right corner of VS Code, confirm the selected Python interpreter is your Anaconda environment (Python 3.11.5) — click it to change if it's showing something else.

**Running the pipeline:**

1. Open `01_Code/sle_pipeline.py` in the editor.
2. Click the **Run** button (the triangle/play icon in the top-right of the editor), or right-click anywhere in the file and choose **Run Python File in Terminal**.
3. Progress prints to the integrated terminal as each section runs, and the same output is written to `02_Output/03_pipeline_log.txt` as it goes — open that file at any point (during or after the run) to check progress or review a finished run later.
4. When it finishes, every table, JSON file, and figure listed in `02_Output/README.md` will be sitting in `02_Output/`.

**Before your first real run**, it's worth doing one quick sanity check: open `configurations.py` and confirm `DATA_PATH` resolves correctly — the easiest way is to run the pipeline once and check the very first lines of `03_pipeline_log.txt`, which log the resolved paths and settings before anything else happens. Notably, this is the one script in the whole repo that **doesn't** have a hardcoded absolute path problem — `configurations.py` derives every path from `Path(__file__).resolve().parent`, so it works correctly regardless of where the repo is cloned.

**A note on runtime**: this pipeline does real work — 200-iteration bootstrap resampling for the final-panel stability check, six feature selectors and six classifiers each with their own hyperparameter grid search across 5 outer folds, then a second 200-iteration bootstrap and a 200-permutation test at the end. Expect this to take a while to complete; let it run rather than assuming it has stalled if the terminal is quiet for a stretch between logged section headers.

---

## Key Result

**Winning combination:** Boruta + SVM (mean outer-fold AUC 0.917). **Final panel: HERC5, BATF2, SPATS2L** — held-out AUC 0.861 (95% CI 0.774–0.921), sensitivity 0.810, specificity 0.788, permutation p=0.005. See [`02_Output/README.md`](./02_Output/README.md) for every metric and plot in full.

## Next Stage

`06_Final_Panel_Validation` confirms this 3-gene panel tracks SLE status specifically, not pregnancy status, complication severity, or batch.
