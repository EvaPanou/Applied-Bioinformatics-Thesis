# --------------------------------------------------------------------
# Helper functions for the SLE diagnostic gene-signature pipeline.
# --------------------------------------------------------------------

# I use this file to keep the main pipeline script cleaner and easier to follow.
# It contains the functions for loading the data, making donor-stratified splits,
# running the feature-selection methods, benchmarking models, building the
# final gene signature, and testing the final model on the held-out donors.

# The overall aim is to support the main ML pipeline for
# classifying SLE versus healthy controls from gene-expression data.
# The main script (sle_pipeline.py) only orchestrates these steps.

# IT WILL NEED TO BE LOCATED IN THE SAME DIRECTORY AS THE MAIN PIPELINE SCRIPT "sle_pipeline.py"
# in order to be directly called into action

# --------------------------
# IMPORT STEP
# --------------------------

from __future__ import annotations # special Python import that lets a future statement for enabling features before they become standard behavior.

from collections import Counter
from pathlib import Path
from typing import Iterable, Tuple

import logging
import os
import random

import numpy as np
import pandas as pd

from collections import Counter
from pathlib import Path
from typing import Iterable, Tuple
import json
import logging
import os
import random

from requirements import (
    StandardScaler,
    RandomForestClassifier,
    StratifiedGroupKFold,
    StratifiedShuffleSplit,
    LogisticRegression,
    LogisticRegressionCV,
    RFECV,
    SVC,
    GaussianNB,
    GridSearchCV,
    accuracy_score,
    roc_auc_score,
    confusion_matrix,
    BorutaPy,
    XGBClassifier,
)

from configurations import (
    LOG_PATH,
    OUTPUT_DIR,
    DATA_PATH,
    METADATA_COLUMNS,
    DONOR_ID_COL,
    LABEL_COL,
    HOLDOUT_FRAC,
    N_OUTER_FOLDS,
    N_INNER_FOLDS,
    N_REPEATS,
    SEED,
    FS_METHODS,
    CLASSIFIERS,
    MRMR_TOP_K,
    LASSO_C_VALUES,
    ELASTICNET_C_VALUES,
    ELASTICNET_L1_RATIOS,
    RF_IMPORTANCE_TREES,
    SVM_RFE_STEP,
    SVM_RFE_MIN_FEATURES,
    SVM_RFE_INNER_FOLDS,
    BORUTA_MAX_ITER,
    LR_C_VALUES,
    SVM_C_VALUES,
    SVM_GAMMA_VALUES,
    RF_DEPTH_VALUES,
    RF_LEAF_VALUES,
    XGB_DEPTH_VALUES,
    XGB_LR_VALUES,
    XGB_N_ESTIMATORS,
    NB_VAR_SMOOTHING,
    KNN_NEIGHBORS,
    FOLD_BOOTSTRAP_ITERATIONS,
    FOLD_BOOTSTRAP_TREES,
    FOLD_BOOTSTRAP_TOP_PERCENTILE,
    FOLD_STABILITY_THRESHOLD,
    FOLD_STABLE_FALLBACK_TOP_N,
    FOLD_MIN_STABLE_GENES,
    FINAL_STABILITY_THRESHOLD,     # ← add this
    FINAL_PEARSON_THRESHOLD,       # ← add this
)

# ------------------------------------------------------------------
# Basic setup functions
# ------------------------------------------------------------------

def setup_logging(log_path: Path = LOG_PATH) -> logging.Logger:
    """
    Set up logging of all actions straight to both a file and the terminal.

    This helps keep a LOG FILE of what happened during the run.
    If the logger already has handlers, they are removed first
    so the same message is not printed twice.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Remove old handlers if they already exist
    for handler in list(logger.handlers):
        logger.removeHandler(handler)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_path, mode="w")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logger.info("Logger started. Writing log to %s", log_path)
    return logger


def set_global_seed(seed: int = SEED) -> None:
    """
    This function sets the random SEED for reproducibility.

    This affects Python's random module, NumPy, and hash behaviour.
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


def ensure_output_dir(path: Path = OUTPUT_DIR) -> Path:
    """
    This function creates the output folder if it does not already exist.
    """
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# Donor-stratification CROSS VALIDATION functions
# ---------------------------------------------------------------------------

def make_donor_stratified_cv(
    n_splits: int = N_OUTER_FOLDS,
    seed: int = SEED,
) -> StratifiedGroupKFold:
    """
    This function creates donor-stratified cross-validation folds.

    Samples from the same donor always stay in the same fold.
    This is important because the dataset contains repeated
    samples from the same donor.
    """
    return StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=seed,
    )


def assert_no_donor_leakage(
    train_donors: Iterable,
    test_donors: Iterable,
    fold_label: str = "",
) -> None:
    """
    This function checks that no donor appears in both train and test.

    If the same donor is found in both sets, the function stops
    the run with an error.
    """
    overlap = set(train_donors) & set(test_donors)

    if overlap:
        example_donor = next(iter(overlap))
        raise RuntimeError(
            f"[{fold_label}] Donor leakage found: {len(overlap)} overlapping donor(s). "
            f"Example donor: {example_donor}"
        )


def summarise_split(
    name: str,
    y: np.ndarray,
    donors: np.ndarray,
    logger: logging.Logger,
) -> None:
    """
    Print a short summary of one dataset split.

    This includes:
    - number of donors
    - number of samples
    - number of SLE samples
    - number of healthy samples
    """
    n_sle = int(np.sum(y == 1))
    n_healthy = int(np.sum(y == 0))

    logger.info(
        "%-12s | donors=%3d | samples=%4d | SLE=%3d | Healthy=%3d",
        name,
        len(np.unique(donors)),
        len(y),
        n_sle,
        n_healthy,
    )

# ---------------------------------------------------------------------------
# Data loading + sealed donor split
# ---------------------------------------------------------------------------

def load_dataset(
    data_path: Path = DATA_PATH,
    logger: logging.Logger | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, list[str]]:
    """
    Load the input file and separate metadata from gene columns.

    Returns:
    - df: full dataframe
    - X: gene expression matrix
    - y: labels (1 = SLE, 0 = Healthy)
    - donors: donor IDs
    - gene_columns: list of gene column names
    """
    log = logger or logging.getLogger()
    path = Path(data_path)

    if not path.exists():
        raise FileNotFoundError(f"Dataset not found: {data_path}")

    # Use tab for TSV files, comma otherwise
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    df = pd.read_csv(path, sep=sep)

    # Check that all expected metadata columns exist
    missing_columns = [col for col in METADATA_COLUMNS if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required metadata columns: {missing_columns}")

    # Everything that is not metadata is treated as a gene column
    gene_columns = [col for col in df.columns if col not in METADATA_COLUMNS]

    X = df[gene_columns].astype(float).copy()
    y = df[LABEL_COL].astype(int).copy()
    donors = df[DONOR_ID_COL].astype(str).copy()

    # Labels should be only 0 and 1
    if sorted(y.unique().tolist()) != [0, 1]:
        raise ValueError(
            f"{LABEL_COL} must contain only 0 and 1. Found: {sorted(y.unique())}"
        )

    log.info("Input file: %s", path)
    log.info(
        "Loaded %d samples, %d donors, %d gene columns",
        len(df),
        donors.nunique(),
        len(gene_columns),
    )
    log.info("Sample class counts: %s", dict(Counter(y.tolist())))

    donor_class_counts = (
        df[[DONOR_ID_COL, LABEL_COL]]
        .drop_duplicates()[LABEL_COL]
        .value_counts()
        .to_dict()
    )
    log.info("Donor class counts: %s", donor_class_counts)

    return df, X, y, donors, gene_columns


def make_sealed_donor_split(
    df: pd.DataFrame,
    X: pd.DataFrame,
    y: pd.Series,
    donors: pd.Series,
    holdout_frac: float = HOLDOUT_FRAC,
    seed: int = SEED,
    logger: logging.Logger | None = None,
) -> dict:
    """
    Create a donor-level development / held-out split.

    The split is done at donor level, not sample level.
    This means a donor can only appear in one set.
    """
    log = logger or logging.getLogger()
    log.info(
        "Creating donor-level held-out split (holdout_frac=%.2f, seed=%d)",
        holdout_frac,
        seed,
    )

    # Build one row per donor for stratified splitting
    donor_table = (
        df[[DONOR_ID_COL, LABEL_COL]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    donor_y = donor_table[LABEL_COL].astype(int)

    splitter = StratifiedShuffleSplit(
        n_splits=1,
        test_size=holdout_frac,
        random_state=seed,
    )

    dev_idx, hold_idx = next(
        splitter.split(np.arange(len(donor_table)), donor_y)
    )

    dev_donors = set(donor_table.iloc[dev_idx][DONOR_ID_COL].astype(str))
    hold_donors = set(donor_table.iloc[hold_idx][DONOR_ID_COL].astype(str))

    # Safety check
    overlap = dev_donors & hold_donors
    if overlap:
        raise RuntimeError(
            f"Donor overlap found between development and held-out sets: {len(overlap)}"
        )

    dev_mask = donors.isin(dev_donors)
    hold_mask = donors.isin(hold_donors)

    split = {
        "X_dev": X.loc[dev_mask].reset_index(drop=True),
        "y_dev": y.loc[dev_mask].reset_index(drop=True),
        "donors_dev": donors.loc[dev_mask].reset_index(drop=True),
        "X_hold": X.loc[hold_mask].reset_index(drop=True),
        "y_hold": y.loc[hold_mask].reset_index(drop=True),
        "donors_hold": donors.loc[hold_mask].reset_index(drop=True),
        "dev_donors": sorted(dev_donors),
        "hold_donors": sorted(hold_donors),
    }

    summarise_split(
        "development",
        split["y_dev"].values,
        split["donors_dev"].values,
        log,
    )
    summarise_split(
        "held-out",
        split["y_hold"].values,
        split["donors_hold"].values,
        log,
    )

    log.info("Held-out donor set is now sealed.")
    return split


# ---------------------------------------------------------------------------
# Stability filter
# ---------------------------------------------------------------------------

def run_bootstrap_rf_stability(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    donors_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    n_iter: int = FOLD_BOOTSTRAP_ITERATIONS,
    n_trees: int = FOLD_BOOTSTRAP_TREES,
    top_percentile: int = FOLD_BOOTSTRAP_TOP_PERCENTILE,
    stability_threshold: float = FOLD_STABILITY_THRESHOLD,
    fallback_top_n: int = FOLD_STABLE_FALLBACK_TOP_N,
    min_stable: int = FOLD_MIN_STABLE_GENES,
    seed: int = SEED,
) -> Tuple[pd.Series, pd.Series]:
    """
    Run a donor-level bootstrap stability filter with Random Forest.

    In each bootstrap iteration, I resample donors with replacement,
    fit a random forest, and record which genes are in the top
    importance group. At the end, I keep genes that appear often enough
    across iterations.

    Returns:
    - stable_genes: genes that passed the stability rule
    - frequency: selection frequency for all genes
    """
    log = logger or logging.getLogger()
    log.info("%s | Running bootstrap RF stability filter (%d iterations)", fold_id, n_iter)

    gene_names = X_train.columns.tolist()
    counts = pd.Series(0.0, index=gene_names)

    unique_donors = donors_train.unique()
    rng = np.random.RandomState(seed)
    valid_iterations = 0

    # Scale once at the beginning
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(
        scaler.fit_transform(X_train),
        columns=gene_names,
        index=X_train.index,
    )

    for i in range(n_iter):
        # Resample donors, not rows
        sampled_donors = rng.choice(
            unique_donors,
            size=len(unique_donors),
            replace=True,
        )

        mask = donors_train.isin(sampled_donors)
        X_boot = X_scaled.loc[mask]
        y_boot = y_train.loc[mask]

        # Skip if one class disappears
        if y_boot.nunique() < 2:
            continue

        rf = RandomForestClassifier(
            n_estimators=n_trees,
            random_state=seed + i,
            class_weight="balanced",
            n_jobs=-1,
        )
        rf.fit(X_boot, y_boot)

        importances = pd.Series(rf.feature_importances_, index=gene_names)

        # Example: top_percentile = 20 means keep genes above the 80th percentile
        cutoff = np.percentile(importances, 100 - top_percentile)
        selected_now = importances[importances >= cutoff].index

        counts.loc[selected_now] += 1
        valid_iterations += 1

        if (i + 1) % 50 == 0 or (i + 1) == n_iter:
            log.info(
                "%s | Bootstrap progress: %d/%d (valid=%d)",
                fold_id, i + 1, n_iter, valid_iterations,
            )

    if valid_iterations == 0:
        raise RuntimeError(f"{fold_id}: no valid bootstrap iterations completed")

    frequency = (counts / valid_iterations).sort_values(ascending=False)
    stable_genes = frequency[frequency >= stability_threshold]

    # Fallback if very few genes survive the threshold
    if len(stable_genes) < min_stable:
        log.info(
            "%s | Too few stable genes (%d). Using fallback top %d genes.",
            fold_id, len(stable_genes), fallback_top_n,
        )
        stable_genes = frequency.head(fallback_top_n)

    log.info("%s | Stable genes kept: %d", fold_id, len(stable_genes))
    log.info(
        "%s | Top stable genes: %s",
        fold_id, ", ".join(stable_genes.head(10).index.tolist()),
    )

    return stable_genes, frequency


# ------------------------------------------------------------------
# Small helper for scaling
# ------------------------------------------------------------------

def _scale(X: pd.DataFrame) -> pd.DataFrame:
    """
    Standardise a dataframe and keep the same row/column labels.
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return pd.DataFrame(X_scaled, columns=X.columns, index=X.index)



# ------------------------------------------------------------------
# Feature selectors
# ------------------------------------------------------------------

def run_mrmr_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    top_k: int = MRMR_TOP_K,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run mRMR feature selection.

    First I try to use the external mrmr package. If that is not available,
    I use a simple fallback based on:
    - mutual information = relevance
    - average absolute correlation = redundancy

    Returns:
    - selected genes
    - ranking of genes from best to worst
    """
    log = logger or logging.getLogger()
    K = min(top_k, X_train.shape[1])
    log.info("%s | Running mRMR (top_k=%d)", fold_id, K)

    try:
        from mrmr import mrmr_classif  # type: ignore

        selected = list(
            mrmr_classif(X=X_train, y=y_train, K=K, show_progress=False)
        )
        log.info("%s | mRMR selected %d genes using the package", fold_id, len(selected))
        return selected, selected

    except Exception as e:
        log.info("%s | mRMR package not available, using fallback (%s)", fold_id, e)

    from sklearn.feature_selection import mutual_info_classif

    X_scaled = _scale(X_train)

    relevance = pd.Series(
        mutual_info_classif(X_scaled, y_train, random_state=seed),
        index=X_train.columns,
    )

    corr_abs = X_scaled.corr().abs().fillna(0.0)

    selected: list[str] = []
    remaining = list(X_train.columns)

    # Start with the most relevant gene
    first_gene = relevance.sort_values(ascending=False).index[0]
    selected.append(first_gene)
    remaining.remove(first_gene)

    while len(selected) < K and remaining:
        redundancy = corr_abs.loc[remaining, selected].mean(axis=1)
        score = relevance.loc[remaining] - redundancy
        best_gene = score.idxmax()

        selected.append(best_gene)
        remaining.remove(best_gene)

    log.info("%s | Fallback mRMR selected %d genes", fold_id, len(selected))
    return selected, selected


def run_lasso_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    C_values: tuple = LASSO_C_VALUES,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run LASSO logistic regression for feature selection.

    Genes with non-zero coefficients are kept.
    All genes are also ranked by absolute coefficient size.
    """
    log = logger or logging.getLogger()
    log.info("%s | Running LASSO", fold_id)

    X_scaled = _scale(X_train).values

    model = LogisticRegressionCV(
        Cs=list(C_values),
        cv=3,
        penalty="l1",
        solver="saga",
        scoring="roc_auc",
        max_iter=5000,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )
    model.fit(X_scaled, y_train)

    coef = pd.Series(model.coef_[0], index=X_train.columns)

    selected = coef[coef != 0].index.tolist()
    ranking = coef.abs().sort_values(ascending=False).index.tolist()

    log.info(
        "%s | LASSO selected %d genes (best C = %.4g)",
        fold_id, len(selected), float(model.C_[0]),
    )

    return selected, ranking


def run_elasticnet_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    C_values: tuple = ELASTICNET_C_VALUES,
    l1_ratios: tuple = ELASTICNET_L1_RATIOS,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run Elastic Net logistic regression for feature selection.

    This is similar to LASSO, but it mixes L1 and L2 penalty.
    Genes with non-zero coefficients are selected.
    """
    log = logger or logging.getLogger()
    log.info("%s | Running Elastic Net", fold_id)

    X_scaled = _scale(X_train).values

    model = LogisticRegressionCV(
        Cs=list(C_values),
        cv=3,
        penalty="elasticnet",
        solver="saga",
        scoring="roc_auc",
        max_iter=5000,
        class_weight="balanced",
        random_state=seed,
        l1_ratios=list(l1_ratios),
        n_jobs=-1,
    )
    model.fit(X_scaled, y_train)

    coef = pd.Series(model.coef_[0], index=X_train.columns)

    selected = coef[coef != 0].index.tolist()
    ranking = coef.abs().sort_values(ascending=False).index.tolist()

    log.info(
        "%s | Elastic Net selected %d genes (best C = %.4g, best l1_ratio = %.2f)",
        fold_id,
        len(selected),
        float(model.C_[0]),
        float(model.l1_ratio_[0]),
    )

    return selected, ranking


def run_svm_rfe_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    donors_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    step: int = SVM_RFE_STEP,
    min_features: int = SVM_RFE_MIN_FEATURES,
    inner_folds: int = SVM_RFE_INNER_FOLDS,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run SVM-RFE with donor-stratified inner cross-validation.

    RFECV removes features step by step and tries to find a good subset size.
    Because the data has repeated samples per donor, I use donor-stratified splits.
    """
    log = logger or logging.getLogger()
    log.info(
        "%s | Running SVM-RFE (step=%d, min_features=%d)",
        fold_id, step, min_features,
    )

    X_scaled = _scale(X_train)

    inner_cv = StratifiedGroupKFold(
        n_splits=inner_folds,
        shuffle=True,
        random_state=seed,
    )
    splits = list(inner_cv.split(X_scaled.values, y_train, groups=donors_train))

    selector = RFECV(
        estimator=SVC(kernel="linear"),
        step=step,
        cv=splits,
        scoring="roc_auc",
        min_features_to_select=min_features,
        n_jobs=-1,
    )
    selector.fit(X_scaled.values, y_train)

    selected = X_train.columns[selector.support_].tolist()

    ranking_series = pd.Series(selector.ranking_, index=X_train.columns)
    ranking = ranking_series.sort_values(ascending=True).index.tolist()

    log.info(
        "%s | SVM-RFE selected %d genes",
        fold_id, len(selected),
    )

    return selected, ranking


def run_rf_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    n_trees: int = RF_IMPORTANCE_TREES,
    top_k: int = MRMR_TOP_K,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run Random Forest feature selection based on feature importance.

    The genes are ranked by RF importance, and the top_k genes are selected.
    """
    log = logger or logging.getLogger()
    log.info("%s | Running RF feature importance", fold_id)

    rf = RandomForestClassifier(
        n_estimators=n_trees,
        random_state=seed,
        class_weight="balanced",
        n_jobs=-1,
    )
    rf.fit(X_train, y_train)

    importance = pd.Series(rf.feature_importances_, index=X_train.columns)
    ranking = importance.sort_values(ascending=False).index.tolist()
    selected = ranking[: min(top_k, len(ranking))]

    log.info("%s | RF selected top %d genes", fold_id, len(selected))
    return selected, ranking


def run_boruta_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    max_iter: int = BORUTA_MAX_ITER,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run Boruta feature selection.

    Boruta compares real features against shuffled shadow features
    and tries to keep all relevant variables.
    """
    log = logger or logging.getLogger()
    log.info("%s | Running Boruta (max_iter=%d)", fold_id, max_iter)

    X_scaled = _scale(X_train).values

    rf = RandomForestClassifier(
        n_estimators=RF_IMPORTANCE_TREES,
        random_state=seed,
        class_weight="balanced",
        n_jobs=-1,
    )

    selector = BorutaPy(
        estimator=rf,
        n_estimators="auto",
        max_iter=max_iter,
        random_state=seed,
        verbose=0,
    )
    selector.fit(X_scaled, y_train.values)

    support = pd.Series(selector.support_, index=X_train.columns)
    ranking_series = pd.Series(selector.ranking_, index=X_train.columns)

    selected = support[support].index.tolist()
    ranking = ranking_series.sort_values().index.tolist()

    log.info("%s | Boruta selected %d genes", fold_id, len(selected))
    return selected, ranking


def run_selector(
    name: str,
    X_train: pd.DataFrame,
    y_train: pd.Series,
    donors_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
) -> Tuple[list[str], list[str]]:
    """
    Run one feature selector by name.

    Returns:
    - selected genes
    - full ranking of genes
    """
    name = name.lower()

    if name == "mrmr":
        return run_mrmr_selector(X_train, y_train, fold_id, logger)

    if name == "lasso":
        return run_lasso_selector(X_train, y_train, fold_id, logger)

    if name == "elasticnet":
        return run_elasticnet_selector(X_train, y_train, fold_id, logger)

    if name == "svm_rfe":
        return run_svm_rfe_selector(
            X_train, y_train, donors_train, fold_id, logger
        )

    if name == "rf":
        return run_rf_selector(X_train, y_train, fold_id, logger)

    if name == "boruta":
        return run_boruta_selector(X_train, y_train, fold_id, logger)

    raise ValueError(
        f"Unknown feature selector: {name!r}. Known values are: {list(FS_METHODS)}"
    )

# ------------------------------------------------------------------
# Elbow step
# ------------------------------------------------------------------

def _auc_curve_donor_cv(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    donors_train: pd.Series,
    ranked_genes: list[str],
    max_k: int,
    n_folds: int,
    seed: int,
) -> np.ndarray:
    """
    Calculate the mean donor-stratified CV AUC for different numbers of genes.

    For each k from 2 up to max_k:
    - keep the top-k genes from the ranking,
    - fit a simple logistic regression model,
    - evaluate mean AUC across donor-stratified folds.

    This function is used only to decide panel size.
    """
    auc_values = []

    cv = StratifiedGroupKFold(
        n_splits=n_folds,
        shuffle=True,
        random_state=seed,
    )

    for k in range(2, max_k + 1):
        selected_genes = ranked_genes[:k]
        X_sub = X_train[selected_genes]
        X_scaled = _scale(X_sub).values

        fold_aucs = []

        for train_idx, valid_idx in cv.split(X_scaled, y_train, groups=donors_train):
            model = LogisticRegression(
                class_weight="balanced",
                max_iter=5000,
                random_state=seed,
            )

            model.fit(X_scaled[train_idx], y_train.iloc[train_idx])
            proba = model.predict_proba(X_scaled[valid_idx])[:, 1]

            fold_auc = roc_auc_score(y_train.iloc[valid_idx], proba)
            fold_aucs.append(fold_auc)

        auc_values.append(float(np.mean(fold_aucs)))

    return np.asarray(auc_values)


def elbow_optimal_subset(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    donors_train: pd.Series,
    ranked_genes: list[str],
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    max_k: int = MRMR_TOP_K,
    min_improvement: float = 0.01,
    n_folds: int = N_INNER_FOLDS,
    seed: int = SEED,
) -> Tuple[list[str], np.ndarray, int]:
    """
    Choose a panel size using a simple elbow rule on the AUC curve.

    I calculate mean donor-astratifiedware CV AUC for k = 2, 3, ..., max_k.
    Then I look for the first point where adding one more gene improves
    AUC by less than min_improvement.

    Returns:
    - chosen genes
    - full AUC curve
    - chosen k
    """
    log = logger or logging.getLogger()

    if len(ranked_genes) <= 2:
        log.info(
            "%s | Elbow skipped because the ranking has only %d genes",
            fold_id, len(ranked_genes),
        )
        return list(ranked_genes), np.array([]), len(ranked_genes)

    max_k_eff = min(max_k, len(ranked_genes))

    log.info(
        "%s | Running elbow step on k = 2 to %d",
        fold_id, max_k_eff,
    )

    auc_curve = _auc_curve_donor_cv(
        X_train=X_train,
        y_train=y_train,
        donors_train=donors_train,
        ranked_genes=ranked_genes,
        max_k=max_k_eff,
        n_folds=n_folds,
        seed=seed,
    )

    improvements = np.diff(auc_curve)
    elbow_k = max_k_eff

    for i, gain in enumerate(improvements):
        if gain < min_improvement:
            elbow_k = i + 2
            break

    chosen_genes = list(ranked_genes[:elbow_k])

    log.info(
        "%s | Elbow chose %d genes",
        fold_id, elbow_k,
    )

    return chosen_genes, auc_curve, elbow_k


# ------------------------------------------------------------------
# Classifier benchmark helpers
# ------------------------------------------------------------------

def _classifier_grid() -> dict:
    """
    Create the classifier grid used in the benchmark step.

    Each entry contains:
    - the classifier object
    - the parameter grid for GridSearchCV
    """
    grid = {}

    if "svm" in CLASSIFIERS:
        grid["svm"] = (
            SVC(
                probability=True,
                class_weight="balanced",
                random_state=SEED,
            ),
            {
                "C": list(SVM_C_VALUES),
                "gamma": list(SVM_GAMMA_VALUES),
            },
        )

    if "logreg" in CLASSIFIERS:
        grid["logreg"] = (
            LogisticRegression(
                class_weight="balanced",
                max_iter=5000,
                random_state=SEED,
            ),
            {
                "C": list(LR_C_VALUES),
            },
        )

    if "rf" in CLASSIFIERS:
        grid["rf"] = (
            RandomForestClassifier(
                class_weight="balanced",
                random_state=SEED,
                n_jobs=-1,
            ),
            {
                "max_depth": list(RF_DEPTH_VALUES),
                "min_samples_leaf": list(RF_LEAF_VALUES),
            },
        )

    if "xgboost" in CLASSIFIERS:
        grid["xgboost"] = (
            XGBClassifier(
                objective="binary:logistic",
                eval_metric="logloss",
                random_state=SEED,
                n_jobs=-1,
            ),
            {
                "max_depth": list(XGB_DEPTH_VALUES),
                "learning_rate": list(XGB_LR_VALUES),
                "n_estimators": list(XGB_N_ESTIMATORS),
            },
        )

    if "nb" in CLASSIFIERS:
        grid["nb"] = (
            GaussianNB(),
            {
                "var_smoothing": list(NB_VAR_SMOOTHING),
            },
        )

    if "knn" in CLASSIFIERS:
        from sklearn.neighbors import KNeighborsClassifier

        grid["knn"] = (
            KNeighborsClassifier(),
            {
                "n_neighbors": list(KNN_NEIGHBORS),
                "weights": ["uniform", "distance"],
            },
        )

    return grid


def _benchmark_classifiers_one_panel(
    X_tr: pd.DataFrame,
    y_tr: pd.Series,
    donors_tr: pd.Series,
    X_va: pd.DataFrame,
    y_va: pd.Series,
    fold_id: str,
    n_inner_folds: int,
    seed: int,
    logger: logging.Logger,
) -> list[dict]:
    """
    Tune and score all classifiers on one fixed gene panel.

    Steps:
    - scale training and validation data,
    - create donor-stratified inner CV folds,
    - run GridSearchCV for each classifier,
    - evaluate the best model on the outer validation fold.

    Returns a list of result rows.
    """
    scaler = StandardScaler()
    X_tr_scaled = scaler.fit_transform(X_tr)
    X_va_scaled = scaler.transform(X_va)

    inner_cv = StratifiedGroupKFold(
        n_splits=n_inner_folds,
        shuffle=True,
        random_state=seed,
    )
    inner_splits = list(inner_cv.split(X_tr_scaled, y_tr, groups=donors_tr))

    results = []

    for clf_name, (estimator, param_grid) in _classifier_grid().items():
        search = GridSearchCV(
            estimator=estimator,
            param_grid=param_grid,
            scoring="roc_auc",
            cv=inner_splits,
            n_jobs=-1,
        )
        search.fit(X_tr_scaled, y_tr)

        best_model = search.best_estimator_
        proba = best_model.predict_proba(X_va_scaled)[:, 1]
        pred = (proba >= 0.5).astype(int)

        auc = float(roc_auc_score(y_va, proba))
        acc = float(accuracy_score(y_va, pred))

        logger.info(
            "%s | %-8s | AUC=%.4f | ACC=%.3f | best=%s",
            fold_id,
            clf_name,
            auc,
            acc,
            search.best_params_,
        )

        results.append(
            {
                "classifier": clf_name,
                "auc": auc,
                "accuracy": acc,
                "best_params": search.best_params_,
            }
        )

    return results


def run_outer_benchmark(
    X_dev: pd.DataFrame,
    y_dev: pd.Series,
    donors_dev: pd.Series,
    logger: logging.Logger | None = None,
    *,
    n_outer_folds: int = N_OUTER_FOLDS,
    n_inner_folds: int = N_INNER_FOLDS,
    n_repeats: int = N_REPEATS,
    fs_methods: tuple = FS_METHODS,
    seed: int = SEED,
) -> Tuple[pd.DataFrame, dict]:
    """
    Run the outer benchmark of (feature selector x classifier).

    For each outer fold:
    1. apply the bootstrap stability filter on the training part,
    2. run each feature selector,
    3. use the elbow step to choose panel size,
    4. benchmark all classifiers on that panel,
    5. store the outer-fold results.

    Returns:
    - results dataframe
    - dictionary with the selected panel in each fold
    """
    log = logger or logging.getLogger()
    log.info(
        "Starting outer benchmark | outer_folds=%d | repeats=%d",
        n_outer_folds,
        n_repeats,
    )

    all_rows = []
    panels_per_fold = {}

    for repeat_idx in range(n_repeats):
        repeat_seed = seed + 1000 * repeat_idx

        outer_cv = StratifiedGroupKFold(
            n_splits=n_outer_folds,
            shuffle=True,
            random_state=repeat_seed,
        )

        for fold_idx, (train_idx, valid_idx) in enumerate(
            outer_cv.split(X_dev, y_dev, groups=donors_dev),
            start=1,
        ):
            fold_id = f"r{repeat_idx + 1}f{fold_idx}"

            X_tr = X_dev.iloc[train_idx].reset_index(drop=True)
            y_tr = y_dev.iloc[train_idx].reset_index(drop=True)
            d_tr = donors_dev.iloc[train_idx].reset_index(drop=True)

            X_va = X_dev.iloc[valid_idx].reset_index(drop=True)
            y_va = y_dev.iloc[valid_idx].reset_index(drop=True)
            d_va = donors_dev.iloc[valid_idx].reset_index(drop=True)

            assert_no_donor_leakage(d_tr, d_va, fold_id)

            log.info(
                "%s | train donors=%d | valid donors=%d",
                fold_id,
                d_tr.nunique(),
                d_va.nunique(),
            )

            # Step 1: stability filter on the training fold only
            stable_genes, _ = run_bootstrap_rf_stability(
                X_train=X_tr,
                y_train=y_tr,
                donors_train=d_tr,
                fold_id=fold_id,
                logger=log,
                seed=repeat_seed,
            )

            X_tr_stable = X_tr[stable_genes.index.tolist()]
            X_va_stable = X_va[stable_genes.index.tolist()]

            # Step 2: try every feature selector
            for fs_name in fs_methods:
                selected, ranking = run_selector(
                    name=fs_name,
                    X_train=X_tr_stable,
                    y_train=y_tr,
                    donors_train=d_tr,
                    fold_id=f"{fold_id}-{fs_name}",
                    logger=log,
                )

                if len(ranking) < 2:
                    log.info(
                        "%s | %s ranking too short, skipping",
                        fold_id,
                        fs_name,
                    )
                    continue

                # Step 3: choose panel size
                panel_genes, auc_curve, elbow_k = elbow_optimal_subset(
                    X_train=X_tr_stable,
                    y_train=y_tr,
                    donors_train=d_tr,
                    ranked_genes=ranking,
                    fold_id=f"{fold_id}-{fs_name}",
                    logger=log,
                    n_folds=n_inner_folds,
                    seed=repeat_seed,
                )

                panels_per_fold[(repeat_idx, fold_idx, fs_name)] = panel_genes

                # Step 4: benchmark classifiers on this panel
                clf_rows = _benchmark_classifiers_one_panel(
                    X_tr=X_tr_stable[panel_genes],
                    y_tr=y_tr,
                    donors_tr=d_tr,
                    X_va=X_va_stable[panel_genes],
                    y_va=y_va,
                    fold_id=f"{fold_id}-{fs_name}",
                    n_inner_folds=n_inner_folds,
                    seed=repeat_seed,
                    logger=log,
                )

                for row in clf_rows:
                    all_rows.append(
                        {
                            "repeat": repeat_idx + 1,
                            "outer_fold": fold_idx,
                            "fs_method": fs_name,
                            "panel_size": elbow_k,
                            **row,
                        }
                    )

    results_df = pd.DataFrame(all_rows)

    log.info("Outer benchmark finished. Total result rows: %d", len(results_df))
    return results_df, panels_per_fold


def summarise_benchmark(
    results_df: pd.DataFrame,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """
    Summarise benchmark results across outer folds.

    The output table contains the mean, standard deviation,
    and number of AUC values for each FS + classifier combination.
    """
    log = logger or logging.getLogger()

    if results_df.empty:
        raise ValueError("results_df is empty, so there is nothing to summarise.")

    summary = (
        results_df
        .groupby(["fs_method", "classifier"])["auc"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .sort_values("mean", ascending=False)
        .reset_index(drop=True)
    )

    log.info("Top benchmark combinations:")
    for _, row in summary.head(5).iterrows():
        log.info(
            "  %s + %s | mean AUC = %.4f | sd = %.4f | n = %d",
            row["fs_method"],
            row["classifier"],
            row["mean"],
            row["std"],
            int(row["count"]),
        )

    return summary

# ------------------------------------------------------------------
# Final model selection
# ------------------------------------------------------------------

def pick_winner(
    results_df: pd.DataFrame,
    logger: logging.Logger | None = None,
) -> dict:
    """
    Pick the best feature-selector + classifier combination.

    I do this by taking the combination with the highest mean outer-fold AUC.
    """
    log = logger or logging.getLogger()

    if results_df.empty:
        raise ValueError("results_df is empty. Run the outer benchmark first.")

    summary = (
        results_df
        .groupby(["fs_method", "classifier"])["auc"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .sort_values("mean", ascending=False)
        .reset_index(drop=True)
    )

    best_row = summary.iloc[0]

    winner = {
        "fs_method": str(best_row["fs_method"]),
        "classifier": str(best_row["classifier"]),
        "mean_auc": float(best_row["mean"]),
        "std_auc": float(best_row["std"]) if not pd.isna(best_row["std"]) else 0.0,
        "n_folds": int(best_row["count"]),
    }

    log.info(
        "Winner: %s + %s | mean AUC = %.4f",
        winner["fs_method"],
        winner["classifier"],
        winner["mean_auc"],
    )

    return winner


def build_final_signature(
    X_dev: pd.DataFrame,
    y_dev: pd.Series,
    donors_dev: pd.Series,
    fs_name: str,
    logger: logging.Logger | None = None,
    *,
    max_k: int = MRMR_TOP_K,
    n_folds: int = N_INNER_FOLDS,
    seed: int = SEED,
) -> dict:
    """
    Refit the winning feature selector on the full development set.

    Then I run the elbow step again on the full development data to choose
    the final number of genes.

    Returns:
    - fs_method
    - signature (list of genes)
    - signature_size
    - ranking
    - auc_curve
    """
    log = logger or logging.getLogger()
    log.info("Building final signature with %s", fs_name)

    _, ranking = run_selector(
        name=fs_name,
        X_train=X_dev,
        y_train=y_dev,
        donors_train=donors_dev,
        fold_id=f"final-{fs_name}",
        logger=log,
    )

    signature, auc_curve, elbow_k = elbow_optimal_subset(
        X_train=X_dev,
        y_train=y_dev,
        donors_train=donors_dev,
        ranked_genes=ranking,
        fold_id=f"final-{fs_name}",
        logger=log,
        max_k=max_k,
        n_folds=n_folds,
        seed=seed,
    )

    final = {
        "fs_method": fs_name,
        "signature": signature,
        "signature_size": len(signature),
        "ranking": ranking,
        "auc_curve": auc_curve.tolist() if len(auc_curve) > 0 else [],
        "elbow_k": elbow_k,
    }

    log.info(
        "Final signature built with %d genes",
        len(signature),
    )

    return final


# ------------------------------------------------------------------
# Final model training
# ------------------------------------------------------------------

def make_classifier(
    classifier_name: str,
    best_params: dict | None = None,
) -> object:
    """
    Create one classifier object by name.

    If best_params is provided, they are passed into the model.
    """
    params = best_params or {}
    name = classifier_name.lower()

    if name == "svm":
        return SVC(
            probability=True,
            class_weight="balanced",
            random_state=SEED,
            **params,
        )

    if name == "logreg":
        return LogisticRegression(
            class_weight="balanced",
            max_iter=5000,
            random_state=SEED,
            **params,
        )

    if name == "rf":
        return RandomForestClassifier(
            class_weight="balanced",
            random_state=SEED,
            n_jobs=-1,
            **params,
        )

    if name == "xgboost":
        return XGBClassifier(
            objective="binary:logistic",
            eval_metric="logloss",
            use_label_encoder=False,
            random_state=SEED,
            n_jobs=-1,
            **params,
        )

    if name == "nb":
        return GaussianNB(**params)

    if name == "knn":
        from sklearn.neighbors import KNeighborsClassifier
        return KNeighborsClassifier(**params)

    raise ValueError(f"Unknown classifier: {classifier_name}")


def train_final_model(
    X_dev: pd.DataFrame,
    y_dev: pd.Series,
    donors_dev: pd.Series,
    final_signature: list[str],
    classifier_name: str,
    logger: logging.Logger | None = None,
    *,
    n_inner_folds: int = N_INNER_FOLDS,
    seed: int = SEED,
) -> dict:
    """
    Train the final model on the full development set.

    I use only the genes in the final signature.
    Hyperparameters are tuned again with donor-stratified inner CV on the
    development set only.

    Returns:
    - fitted scaler
    - fitted model
    - best parameters
    - inner CV AUC
    - signature genes
    """
    log = logger or logging.getLogger()
    log.info(
        "Training final %s model on %d signature genes",
        classifier_name,
        len(final_signature),
    )

    X_panel = X_dev[final_signature]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_panel)

    inner_cv = StratifiedGroupKFold(
        n_splits=n_inner_folds,
        shuffle=True,
        random_state=seed,
    )
    inner_splits = list(inner_cv.split(X_scaled, y_dev, groups=donors_dev))

    classifier_grid = _classifier_grid()
    if classifier_name not in classifier_grid:
        raise ValueError(f"Classifier {classifier_name!r} is not available")

    estimator, param_grid = classifier_grid[classifier_name]

    search = GridSearchCV(
        estimator=estimator,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=inner_splits,
        n_jobs=-1,
    )
    search.fit(X_scaled, y_dev)

    fitted = {
        "scaler": scaler,
        "model": search.best_estimator_,
        "best_params": search.best_params_,
        "inner_cv_auc": float(search.best_score_),
        "signature": list(final_signature),
        "classifier": classifier_name,
    }

    log.info(
        "Final model trained | best params = %s | inner CV AUC = %.4f",
        search.best_params_,
        float(search.best_score_),
    )

    return fitted


# ------------------------------------------------------------------
# Held-out evaluation
# ------------------------------------------------------------------

from sklearn.metrics import f1_score, roc_auc_score, accuracy_score, confusion_matrix
import numpy as np

def evaluate_on_holdout(
    fitted: dict,
    X_hold: pd.DataFrame,
    y_hold: pd.Series,
    logger: logging.Logger | None = None,
) -> dict:
    log = logger or logging.getLogger()

    genes = fitted["signature"]
    scaler = fitted["scaler"]
    model = fitted["model"]

    X_panel = X_hold[genes]
    X_scaled = scaler.transform(X_panel)

    probability = model.predict_proba(X_scaled)[:, 1]
    prediction = (probability >= 0.5).astype(int)

    auc = float(roc_auc_score(y_hold, probability))
    accuracy = float(accuracy_score(y_hold, prediction))
    confusion = confusion_matrix(y_hold, prediction)

    tn, fp, fn, tp = confusion.ravel()

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    f1 = float(f1_score(y_hold, prediction))

    log.info(
        "Held-out results | AUC = %.4f | ACC = %.3f | n = %d",
        auc,
        accuracy,
        len(y_hold),
    )

    return {
        "auc": auc,
        "accuracy": accuracy,
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "f1": f1,
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "tn": int(tn),
        "confusion_matrix": confusion.tolist(),
        "y_true": y_hold.tolist(),
        "y_pred": prediction.tolist(),
        "y_proba": probability.tolist(),
        "n_samples": int(len(y_hold)),
        "n_genes": int(len(genes)),
        "genes": list(genes),
    }

# ------------------------------------------------------------------
# Small save helper
# ------------------------------------------------------------------

def save_json(
    obj: dict,
    path: Path,
    logger: logging.Logger | None = None,
) -> None:
    """
    Save a dictionary as a JSON file.
    """
    log = logger or logging.getLogger()

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)

    log.info("Saved JSON file to %s", path)



#----------------------------------------------------
# TEST phase metrics
# ---------------------------------------------------



def bootstrap_auc_ci(
    ytrue: np.ndarray,                                              # Array of true binary labels from the held-out set.
    yproba: np.ndarray,                                             # Array of predicted probabilities for the positive class.
    nbootstrap: int = 1000,                                         # Number of bootstrap resamples to draw.
    seed: int = SEED,                                               # Random seed for reproducibility.
    logger: logging.Logger | None = None,                           # Optional logger for progress messages.
) -> dict:                                                          # Return a dictionary with CI results.
    """Compute a bootstrap confidence interval for held-out AUC."""  # Short description of the function.


    log = logger or logging.getLogger()                             # Use the provided logger, or fall back to the root logger.


    ytrue = np.asarray(ytrue)                                       # Convert labels to a NumPy array.
    yproba = np.asarray(yproba)                                     # Convert probabilities to a NumPy array.


    if len(ytrue) != len(yproba):                                   # Check that both arrays have the same length.
        raise ValueError("ytrue and yproba must have the same length.")  # Stop if lengths do not match.


    rng = np.random.RandomState(seed)                               # Create a reproducible random-number generator.


    observedauc = float(roc_auc_score(ytrue, yproba))               # Compute the observed AUC on the original held-out predictions.


    bootaucs = []                                                   # Create an empty list to store bootstrap AUC values.


    for i in range(nbootstrap):                                     # Repeat the bootstrap procedure nbootstrap times.
        idx = rng.choice(np.arange(len(ytrue)), size=len(ytrue), replace=True)  # Sample held-out rows with replacement.
        if len(np.unique(ytrue[idx])) < 2:                          # Skip bootstrap samples that contain only one class.
            continue                                                # Move to the next bootstrap iteration.
        auci = float(roc_auc_score(ytrue[idx], yproba[idx]))        # Compute AUC on the resampled labels and probabilities.
        bootaucs.append(auci)                                       # Store the bootstrap AUC.


    cilow = float(np.percentile(bootaucs, 2.5)) if bootaucs else None   # Compute the 2.5th percentile if bootstrap values exist.
    cihigh = float(np.percentile(bootaucs, 97.5)) if bootaucs else None # Compute the 97.5th percentile if bootstrap values exist.


    result = {                                                      # Build the output dictionary.
        "observed_auc": observedauc,                                # Store the observed held-out AUC.
        "ci_lo": cilow,                                             # Store the lower CI bound.
        "ci_hi": cihigh,                                            # Store the upper CI bound.
        "n_bootstraps": int(nbootstrap),                            # Store the requested number of bootstrap iterations.
        "n_valid_bootstrap": int(len(bootaucs)),                    # Store the number of valid bootstrap samples actually used.
        "bootstrap_aucs": [float(x) for x in bootaucs],             # Store all bootstrap AUC values.
    }                                                               # End dictionary.


    log.info(                                                       # Write a summary line to the log.
        "Bootstrap AUC CI | observed AUC %.4f | 95%% CI [%.4f, %.4f] | valid bootstrap %d",  # Log message template.
        observedauc,                                                # First value in the log message.
        cilow if cilow is not None else float("nan"),               # Lower CI value or NaN if missing.
        cihigh if cihigh is not None else float("nan"),             # Upper CI value or NaN if missing.
        len(bootaucs),                                              # Number of valid bootstrap samples.
    )                                                               # End logging call.


    return result                                                   # Return the bootstrap CI results.


def permutation(
    Xdev: pd.DataFrame,                                             # Development-set gene matrix.
    ydev: pd.Series,                                                # Development-set labels.
    donorsdev: pd.Series,                                           # Development-set donor IDs.
    finalsignature: list[str],                                      # Final gene panel used for the permutation test.
    classifiername: str,                                            # Name of the final classifier.
    bestparams: dict | None = None,                                 # Best hyperparameters for the final classifier.
    logger: logging.Logger | None = None,                           # Optional logger.
    npermutations: int = 1000,                                      # Number of donor-level permutations.
    nfolds: int = N_INNER_FOLDS,                                      # Number of donor-stratified CV folds.
    seed: int = SEED,                                               # Random seed.
) -> dict:                                                          # Return a dictionary with permutation-test results.
    """Run a donor-stratified permutation test on the development set."""  # Short description.

    log = logger or logging.getLogger()                             # Use provided logger or root logger.

    Xpanel = Xdev[finalsignature].copy()                            # Restrict development data to the final gene panel.

    cv = StratifiedGroupKFold(                                      # Build donor-stratified cross-validation splits.
        n_splits=nfolds,                                            # Number of folds.
        shuffle=True,                                               # Shuffle donors before splitting.
        random_state=seed,                                          # Seed for reproducibility.
    )                                                               # End CV object.

    def meancvauc(yvector: pd.Series) -> float:                     # Nested helper to compute mean CV AUC for one label vector.
        foldaucs = []                                               # Create a list to store fold AUCs.

        for trainidx, valididx in cv.split(Xpanel, yvector, groups=donorsdev):  # Loop over donor-stratified CV folds.
            Xtr = Xpanel.iloc[trainidx]                             # Extract training genes for this fold.
            Xva = Xpanel.iloc[valididx]                             # Extract validation genes for this fold.
            ytr = yvector.iloc[trainidx]                            # Extract training labels for this fold.
            yva = yvector.iloc[valididx]                            # Extract validation labels for this fold.

            scaler = StandardScaler()                               # Create a fresh scaler inside the fold.
            Xtrscaled = scaler.fit_transform(Xtr)                   # Fit the scaler on training data only.
            Xvascaled = scaler.transform(Xva)                       # Apply the same scaler to validation data.

            model = make_classifier(classifiername, bestparams)      # Recreate the classifier with the chosen hyperparameters.
            model.fit(Xtrscaled, ytr)                               # Fit the model on the training fold.

            proba = model.predict_proba(Xvascaled)[:, 1]            # Get predicted probabilities on the validation fold.
            foldauc = float(roc_auc_score(yva, proba))                # Compute AUC for this fold.
            foldaucs.append(foldauc)                                # Store the fold AUC.

        return float(np.mean(foldaucs))                             # Return the mean AUC across folds.

    observedauc = meancvauc(ydev)                                   # Compute the observed donor-stratified CV AUC on the true labels.

    donortable = pd.DataFrame({                                     # Build a one-row-per-donor table.
        "donor": donorsdev.astype(str).values,                      # Donor IDs as strings.
        "y": ydev.values,                                           # Original labels.
    }).drop_duplicates(subset=["donor"])                            # Keep one row per donor only.

    rng = np.random.RandomState(seed)                               # Create reproducible RNG for permutations.
    nullaucs = []                                                   # Create a list to store null AUC values.

    for i in range(npermutations):                                  # Repeat donor-level label shuffling.
        shuffledlabels = rng.permutation(donortable["y"].values)    # Shuffle donor labels.
        donortolabel = dict(zip(donortable["donor"], shuffledlabels))  # Map each donor to a permuted label.
        yperm = donorsdev.astype(str).map(donortolabel).astype(int) # Expand permuted donor labels back to sample level.

        try:                                                        # Try to evaluate this permutation.
            auci = meancvauc(yperm)                                 # Compute mean CV AUC under the permuted labels.
            nullaucs.append(float(auci))                            # Store the null AUC.
        except Exception as e:                                      # Catch rare fold failures.
            log.warning("Permutation %d skipped: %s", i + 1, e)     # Write a warning and skip this permutation.

    pvalue = (1 + sum(a >= observedauc for a in nullaucs)) / (1 + len(nullaucs))  # Compute permutation p-value with +1 correction.

    result = {                                                      # Build the output dictionary.
        "observed_auc": float(observedauc),                         # Store observed AUC.
        "null_aucs": [float(x) for x in nullaucs],                  # Store null distribution.
        "p_value": float(pvalue),                                   # Store corrected p-value.
        "n_permutations": int(len(nullaucs)),                       # Store number of successful permutations.
    }                                                               # End dictionary.

    log.info(                                                       # Log a summary line.
        "Permutation test | observed AUC %.4f | p %.4f | successful permutations %d",  # Log message template.
        result["observed_auc"],                                     # Observed AUC.
        result["p_value"],                                          # Permutation p-value.
        result["n_permutations"],                                   # Number of successful permutations.
    )                                                               # End logging call.

    return result                                                   # Return the permutation-test results.

def prune_correlated_genes(
    X: pd.DataFrame,
    ranked_genes: list[str],
    correlation_threshold: float,
    logger: logging.Logger | None = None,
    label: str = "final_panel",
) -> list[str]:
    """
    Keep genes in ranked order, dropping later genes if they are too correlated
    with any already-kept earlier gene.
    """
    log = logger or logging.getLogger()

    if len(ranked_genes) <= 1:
        return list(ranked_genes)

    corr_matrix = X[ranked_genes].corr(method="pearson")
    corr_abs = corr_matrix.abs()
    kept_genes: list[str] = []
    pruned_rows: list[dict] = []

    log.info(
        "%s Pearson pruning: %d input genes at |r| >= %.2f",
        label,
        len(ranked_genes),
        correlation_threshold,
    )

    for gene in ranked_genes:
        if not kept_genes:
            kept_genes.append(gene)
            continue

        # correlations of this candidate with already-kept genes
        corr_series = corr_abs.loc[gene, kept_genes]
        max_corr = corr_series.max()
        if max_corr >= correlation_threshold:
            # find the kept gene with which correlation is maximal
            partner = corr_series.idxmax()
            log.info(
                "%s | pruned %-12s |r|=%.3f with %s",
                label,
                gene,
                float(max_corr),
                partner,
            )
            pruned_rows.append(
                {
                    "pruned_gene": gene,
                    "kept_partner": partner,
                    "pearson_r": float(corr_matrix.loc[gene, partner]),
                    "abs_r": float(max_corr),
                }
            )
        else:
            kept_genes.append(gene)

    log.info(
        "%s Pearson pruning kept %d of %d genes at |r| < %.2f",
        label,
        len(kept_genes),
        len(ranked_genes),
        correlation_threshold,
    )

    # Optional: return pruned_rows as well if you want to save them outside.
    if pruned_rows:
        X.attrs[f"{label}_pearson_pruned"] = pruned_rows

    return kept_genes

def finalise_panel(
    results_df: pd.DataFrame,
    X_dev: pd.DataFrame,
    y_dev: pd.Series,
    donors_dev: pd.Series,
    logger: logging.Logger | None = None,
    stability_n_iter: int = FOLD_BOOTSTRAP_ITERATIONS,
    stability_n_trees: int = FOLD_BOOTSTRAP_TREES,
    stability_top_percentile: int = FOLD_BOOTSTRAP_TOP_PERCENTILE,
    stability_fallback_top_n: int = FOLD_STABLE_FALLBACK_TOP_N,
    stability_min_stable: int = FOLD_MIN_STABLE_GENES,
    seed: int = SEED,
) -> dict:
    """
    Choose the winning FS+classifier, rebuild the final signature on the full
    development set, then apply full-development bootstrap RF stability and
    final Pearson redundancy pruning to produce the final panel.
    """
    log = logger or logging.getLogger()

    # 1. Winner from outer benchmark
    winner = pick_winner(results_df, logger=log)

    # 2. Refit winning FS on full dev to get full ranking + elbow info
    final_signature_info = build_final_signature(
        X_dev=X_dev,
        y_dev=y_dev,
        donors_dev=donors_dev,
        fs_name=winner["fs_method"],
        logger=log,
        seed=seed,
    )

    full_ranking = list(final_signature_info["ranking"])

    # 3. Full-development bootstrap RF stability using FINAL_STABILITY_THRESHOLD
    stable_genes, stability_freq = run_bootstrap_rf_stability(
        X_train=X_dev,
        y_train=y_dev,
        donors_train=donors_dev,
        fold_id="final_panel",
        logger=log,
        n_iter=stability_n_iter,
        n_trees=stability_n_trees,
        top_percentile=stability_top_percentile,
        stability_threshold=FINAL_STABILITY_THRESHOLD,
        fallback_top_n=stability_fallback_top_n,
        min_stable=stability_min_stable,
        seed=seed,
    )

    # genes that passed the final stability rule, ordered by FS ranking
    high_freq_gene_set = set(stable_genes.index.tolist())
    high_freq_ranked_genes = [g for g in full_ranking if g in high_freq_gene_set]

    if not high_freq_ranked_genes:
        raise RuntimeError(
            "No genes survived the final stability threshold. "
            "Consider lowering FINAL_STABILITY_THRESHOLD."
        )

    # 4. Pearson redundancy pruning, preserving ranking
    final_panel = prune_correlated_genes(
        X=X_dev[high_freq_ranked_genes],
        ranked_genes=high_freq_ranked_genes,
        correlation_threshold=FINAL_PEARSON_THRESHOLD,
        logger=log,
        label="final_panel",
    )

    log.info(
        "Final panel intersection | ranked=%d stable=%d ranked_stable=%d",
        len(full_ranking),
        len(stable_genes),
        len(high_freq_ranked_genes),
    )

    if not final_panel:
        raise RuntimeError("Pearson pruning removed all genes from the final panel.")

    # 5. Panel table for export
    panel_table = pd.DataFrame(
        {
            "gene": final_panel,
            "rank_in_ranking": [full_ranking.index(gene) + 1 for gene in final_panel],
            "bootstrap_selection_freq": [
                float(stability_freq.get(gene, 0.0)) for gene in final_panel
            ],
        }
    ).sort_values("rank_in_ranking").reset_index(drop=True)

    log.info(
        "Final panel built FS %s classifier %s panel size %d",
        winner["fs_method"],
        winner["classifier"],
        len(final_panel),
    )

    return {
        "winner": winner,
        "panel": final_panel,
        "panel_table": panel_table,
        "high_freq_ranked_genes": high_freq_ranked_genes,
        "stability": stability_freq,
        "stable_genes": stable_genes,
        "ranking": full_ranking,
        "elbow_k": int(final_signature_info["elbow_k"]),
        "auc_curve": list(final_signature_info["auc_curve"]),
    }                                    

