"""
==============================================================================
helper_functions.py
-----------------------------------------------------------------------------
WHAT THIS FILE DOES:
    I use this file to keep the main pipeline script cleaner and easier to follow.
    It contains the functions for loading the data, making donor-stratified splits,
    running the feature-selection methods, benchmarking models, building the
    final gene signature, and testing the final model on the held-out donors.

    The overall aim is to support the main ML pipeline for
    classifying SLE versus healthy controls from gene-expression data.
    The main script (sle_pipeline.py) only orchestrates these steps.

WHERE THIS FILE MUST BE SAVED:
    In the same directory as sle_pipeline.py, otherwise the imports 
    from/to other scripts will not work.

BEFORE RUNNING:
    Open a terminal window within this dir, and do an import check 
    (to catch syntax and naming inconsistencies from the previous 
    scripts to this one), by running the following line: 

    python -m py_compile helper_functions.py  --> for syntax check
    python -c "import helper_functions"       --> imports check
=============================================================================
"""

# --------------------------
# IMPORT STEP
# --------------------------

# importing widely used packages

from collections import Counter
from pathlib import Path
from typing import Iterable, Tuple
import json
import logging

import numpy as np
import pandas as pd

# importing everything configured in my previous scripts

from requirements import *
from configurations import *

# ------------------------------------------------------------------
# Basic setup functions
# ------------------------------------------------------------------

def setup_logging(log_path: Path = LOG_PATH):
    """
    Set up a logger so every message goes to both a log file and the terminal.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True) # setting up parent directory to create log file

    logging.basicConfig(                                            # configure logging settings
        level=logging.INFO,                                         # set up lesser level of message : INFO (as opposed to error or warning) 
        format="%(asctime)s | %(levelname)-7s | %(message)s",       # set up datetime format + message
        datefmt="%Y-%m-%d %H:%M",
        handlers=[
            logging.FileHandler(log_path, mode="w"),                # set up handlers for file and terminal text writing
            logging.StreamHandler(),                                
        ],
    )

    log = logging.getLogger()
    log.info("Logger started. Writing log to %s", log_path)
    return log

# ---------------------------------------------------------------------------
# Donor-stratification CROSS VALIDATION functions
# ---------------------------------------------------------------------------

def donor_stratified_cv(n_splits= N_OUTER_FOLDS, seed = SEED):
    """
    This function creates donor-stratified cross-validation folds.
    Stratified = so that the Healthy vs SLE class ratios are preserved per fold.
    Samples from the same donor always stay in the same fold to retain train-test independence.
    This is to be used within all other functions that make stratified splits.
    """
    return StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=seed,
    )


def assert_no_donor_leakage(train_donors, test_donors, fold_label):
    """
    This function checks that no donor appears in both train and test, as a safety guard.
    If the same donor is found in both sets, the function stops the run with an error.

    1. train_donors & test_donors are iterable 
    2. Group overlap below is written with "&", meaning INTERSECTION, so it should always be empty.
    """
    overlap = set(train_donors) & set(test_donors)      
    
    if overlap:   # if not empty
        raise RuntimeError(
            f"[{fold_label}] Donor leakage found: {len(overlap)} overlapping donor(s). "
        )


def summarise_split(name, y, donors, logger):
    """
    Log a short summary of one dataset split:
    number of donors, samples, SLE samples, and healthy samples.

    name   = label for the split, e.g. "development" or "held-out"
    y      = 0/1 label for each sample (1 = SLE, 0 = Healthy)
    donors = donor ID for each sample
    """
    # count how many samples are SLE (1) and how many are Healthy (0)
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
    holdout_frac: float = HOLDOUT_FRACTION,
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

def run_bootstrap_boruta_stability(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    donors_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    n_iter: int = N_BOOTSTRAP,
    stability_threshold: float = FINAL_STABILITY_THRESHOLD,
    fallback_top_n: int = FOLD_STABLE_FALLBACK_TOP_N,
    min_stable: int = FOLD_MIN_STABLE_GENES,
    max_iter_boruta: int = BORUTA_MAX_ITER,
    seed: int = SEED,
) -> Tuple[pd.Series, pd.Series]:
    """
    Run a donor-level bootstrap stability check using Boruta's own
    confirmed/rejected decision, instead of a raw top-percentile RF
    importance cutoff.

    Kursa (2014) - the same paper this project's stability methodology is
    based on - found Boruta to be the most self-consistent of the
    RF-based gene selection approaches it compared. Since Boruta is also
    the winning method here, resampling and rerunning Boruta itself (not
    a generic RF importance cutoff) is the more literature-consistent
    way to check how stable its selection actually is.

    In each bootstrap iteration, I resample donors with replacement and
    rerun Boruta on that resampled data. At the end, I keep genes that
    Boruta confirms often enough across iterations.

    Calls run_boruta_selector(), defined later in this file, in the
    Feature selectors section.

    Returns:
    - stable_genes: genes that passed the stability rule
    - frequency: confirmation frequency for all genes
    """
    log = logger or logging.getLogger()
    log.info("%s | Running bootstrap Boruta stability filter (%d iterations)", fold_id, n_iter)

    # A quiet logger for the 200 inner Boruta calls, so their own
    # "Running Boruta..." / "Boruta selected N genes" lines don't flood
    # the real log - we only want our own progress lines below.
    quiet_log = logging.getLogger("bootstrap_boruta_inner")
    quiet_log.setLevel(logging.WARNING)

    gene_names = X_train.columns.tolist()
    confirm_counts = pd.Series(0.0, index=gene_names)

    unique_donors = donors_train.unique()
    rng = np.random.RandomState(seed)
    valid_iterations = 0

    for i in range(n_iter):
        # Resample donors, not rows
        sampled_donors = rng.choice(
            unique_donors,
            size=len(unique_donors),
            replace=True,
        )

        mask = donors_train.isin(sampled_donors)
        X_boot = X_train.loc[mask]
        y_boot = y_train.loc[mask]

        # Skip if one class disappears
        if y_boot.nunique() < 2:
            continue

        try:
            confirmed, _ = run_boruta_selector(
                X_train=X_boot,
                y_train=y_boot,
                fold_id=f"{fold_id}-boot{i + 1}",
                logger=quiet_log,
                max_iter=max_iter_boruta,
                seed=seed + i,
            )
        except Exception as e:
            log.warning("%s | Bootstrap iteration %d failed, skipping (%s)", fold_id, i + 1, e)
            continue

        confirm_counts.loc[confirmed] += 1
        valid_iterations += 1

        if (i + 1) % 10 == 0 or (i + 1) == n_iter:
            log.info(
                "%s | Bootstrap progress: %d/%d (valid=%d)",
                fold_id, i + 1, n_iter, valid_iterations,
            )

    if valid_iterations == 0:
        raise RuntimeError(f"{fold_id}: no valid bootstrap iterations completed")

    frequency = (confirm_counts / valid_iterations).sort_values(ascending=False)
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


def _knee_cutoff(
    scores: pd.Series,
    fold_id: str,
    method_label: str,
    logger: logging.Logger,
    fallback_n: int = 10,
) -> int:
    """
    Find a natural cutoff point on a descending score curve using knee detection.

    scores       = one numeric score per gene (e.g. RF importance, mRMR relevance)
    fallback_n   = how many genes to keep if no clear knee is found

    If KneeLocator cannot find a genuine bend in the curve (this happens on
    flat or very noisy score curves), we fall back to a small fixed number
    instead of crashing the run, and log exactly when this happens.
    """
    sorted_scores = scores.sort_values(ascending=False)
    n_available = len(sorted_scores)

    kl = KneeLocator(
        range(1, n_available + 1),
        sorted_scores.values,
        curve="convex",
        direction="decreasing",
    )

    if kl.knee is not None:
        k = int(kl.knee)
        logger.info(
            "%s | %s knee found at k=%d (out of %d candidate genes)",
            fold_id, method_label, k, n_available,
        )
    else:
        k = min(fallback_n, n_available)
        logger.info(
            "%s | %s: No elbow point in the curve, so fallback to %d genes took place (out of %d candidate genes).",
            fold_id, method_label, k, n_available,
        )

    return k



# ------------------------------------------------------------------
# Feature selectors
# ------------------------------------------------------------------

def run_mrmr_selector(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    fold_id: str = "fold?",
    logger: logging.Logger | None = None,
    *,
    fallback_n: int = 10,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run mRMR feature selection.

    First I try to use the external mrmr package. If that is not available,
    I use a simple fallback based on:
    - mutual information = relevance
    - average absolute correlation = redundancy

    Instead of a fixed number of genes, the panel size is chosen with knee
    detection on the relevance score curve, so mRMR decides its own cutoff
    independently of every other feature-selection method.

    Returns:
    - selected genes
    - ranking of genes from best to worst
    """
    log = logger or logging.getLogger()
    n_genes = X_train.shape[1]
    log.info("%s | Running mRMR (ranking all %d candidate genes)", fold_id, n_genes)

    try:
        from mrmr import mrmr_classif  # type: ignore

        # Rank ALL genes (K = n_genes) instead of stopping at a fixed number,
        # and ask for the relevance scores so we can knee-detect the cutoff.
        full_order, relevance, _ = mrmr_classif(
            X=X_train, y=y_train, K=n_genes, return_scores=True, show_progress=False
        )

        k = _knee_cutoff(relevance, fold_id, "mRMR (package, F-statistic relevance)", log, fallback_n)
        selected = list(full_order[:k])

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

    k = _knee_cutoff(relevance, fold_id, "mRMR (fallback, mutual information relevance)", log, fallback_n)

    corr_abs = X_scaled.corr().abs().fillna(0.0)

    selected: list[str] = []
    remaining = list(X_train.columns)

    # Start with the most relevant gene
    first_gene = relevance.sort_values(ascending=False).index[0]
    selected.append(first_gene)
    remaining.remove(first_gene)

    while len(selected) < k and remaining:
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
    fallback_n: int = 10,
    seed: int = SEED,
) -> Tuple[list[str], list[str]]:
    """
    Run Random Forest feature selection based on feature importance.

    The genes are ranked by RF importance, and a natural cutoff point is
    found with knee detection instead of a fixed top-k number, so RF
    importance decides its own panel size independently of the other
    feature-selection methods.
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

    k = _knee_cutoff(importance, fold_id, "RF importance", log, fallback_n)
    selected = ranking[:k]

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
# Classifier benchmark helpers
# ------------------------------------------------------------------

def _classifier_grid(y_train: pd.Series) -> dict:
    """
    Create the classifier grid used in the benchmark step.

    Each entry contains:
    - the classifier object
    - the parameter grid for GridSearchCV

    y_train is used to compute XGBoost's scale_pos_weight, since XGBoost
    has no class_weight parameter the way the other classifiers do.
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
        scale_pos_weight = (y_train == 0).sum() / (y_train == 1).sum()
        grid["xgboost"] = (
            XGBClassifier(
                objective="binary:logistic",
                eval_metric="logloss",
                random_state=SEED,
                n_jobs=-1,
                scale_pos_weight=scale_pos_weight,
            ),
            {
                "max_depth": list(XGB_DEPTH_VALUES),
                "learning_rate": list(XGB_LEARN_VALUES),
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

    for clf_name, (estimator, param_grid) in _classifier_grid(y_tr).items():
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
    1. run each feature selector directly on the fold's full gene set,
    2. each feature selector decides its own panel size,
    3. benchmark all classifiers on that panel,
    4. store the outer-fold results.

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

            # Step 1: try every feature selector directly on the fold's full
            # gene set. (A fold-level bootstrap stability pre-filter used to
            # run here first, narrowing the candidate pool before selection -
            # removed after reviewing real results on this 44-gene panel,
            # where the filter's 70% threshold rarely passed naturally and
            # the fallback path was doing most of the work instead of the
            # intended stability logic.)
            for fs_name in fs_methods:
                selected, ranking = run_selector(
                    name=fs_name,
                    X_train=X_tr,
                    y_train=y_tr,
                    donors_train=d_tr,
                    fold_id=f"{fold_id}-{fs_name}",
                    logger=log,
                )

                if len(selected) < 2:
                    log.info(
                        "%s | %s selected too few genes, skipping",
                        fold_id,
                        fs_name,
                    )
                    continue

                # Step 2: each method above already decided its own panel size -
                # no separate elbow-curve re-sizing step here anymore.
                panels_per_fold[(repeat_idx, fold_idx, fs_name)] = selected

                # Step 3: benchmark classifiers on this panel
                clf_rows = _benchmark_classifiers_one_panel(
                    X_tr=X_tr[selected],
                    y_tr=y_tr,
                    donors_tr=d_tr,
                    X_va=X_va[selected],
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
                            "panel_size": len(selected),
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
    seed: int = SEED,
) -> dict:
    """
    Refit the winning feature selector on the full development set.

    Unlike the outer-fold benchmark loop, this keeps the winning method's
    own actual selected genes directly - no separate elbow-curve re-sizing
    step, consistent with how Phase 2 now works.

    Returns:
    - fs_method
    - signature (the winning method's own selected genes)
    - signature_size
    - ranking (full ordered ranking of every candidate gene, kept for
      display/ordering purposes downstream, NOT as the selection itself)
    """
    log = logger or logging.getLogger()
    log.info("Building final signature with %s", fs_name)

    selected, ranking = run_selector(
        name=fs_name,
        X_train=X_dev,
        y_train=y_dev,
        donors_train=donors_dev,
        fold_id=f"final-{fs_name}",
        logger=log,
    )

    final = {
        "fs_method": fs_name,
        "signature": selected,
        "signature_size": len(selected),
        "ranking": ranking,
    }

    log.info("Final signature built with %d genes", len(selected))

    return final


# ------------------------------------------------------------------
# Final model training
# ------------------------------------------------------------------

def make_classifier(
    classifier_name: str,
    best_params: dict | None = None,
    y_train: pd.Series | None = None,
) -> object:
    """
    Create one classifier object by name.

    If best_params is provided, they are passed into the model.
    y_train is only needed for XGBoost, to compute scale_pos_weight
    for class-imbalance handling (XGBoost has no class_weight parameter).
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
        if y_train is None:
            raise ValueError("make_classifier('xgboost', ...) needs y_train to compute scale_pos_weight.")
        scale_pos_weight = (y_train == 0).sum() / (y_train == 1).sum()
        return XGBClassifier(
            objective="binary:logistic",
            eval_metric="logloss",
            use_label_encoder=False,
            random_state=SEED,
            n_jobs=-1,
            scale_pos_weight=scale_pos_weight,
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

    classifier_grid = _classifier_grid(y_dev)
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



# ------------------------------------------------------------------
# Held-out evaluation statistics (bootstrap CI, permutation test)
# ------------------------------------------------------------------

def bootstrap_auc_ci(
    y_true: np.ndarray,
    y_proba: np.ndarray,
    n_bootstrap: int = 1000,
    seed: int = SEED,
    logger: logging.Logger | None = None,
) -> dict:
    """Compute a bootstrap confidence interval for held-out AUC."""
    log = logger or logging.getLogger()

    y_true = np.asarray(y_true)
    y_proba = np.asarray(y_proba)

    if len(y_true) != len(y_proba):
        raise ValueError("y_true and y_proba must have the same length.")

    rng = np.random.RandomState(seed)
    observed_auc = float(roc_auc_score(y_true, y_proba))

    boot_aucs = []
    for i in range(n_bootstrap):
        idx = rng.choice(np.arange(len(y_true)), size=len(y_true), replace=True)
        if len(np.unique(y_true[idx])) < 2:
            # a resample with only one class has no AUC - skip it
            continue
        auc_i = float(roc_auc_score(y_true[idx], y_proba[idx]))
        boot_aucs.append(auc_i)

    ci_lo = float(np.percentile(boot_aucs, 2.5)) if boot_aucs else None
    ci_hi = float(np.percentile(boot_aucs, 97.5)) if boot_aucs else None

    result = {
        "observed_auc": observed_auc,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "n_bootstraps": int(n_bootstrap),
        "n_valid_bootstrap": int(len(boot_aucs)),
        "bootstrap_aucs": [float(x) for x in boot_aucs],
    }

    log.info(
        "Bootstrap AUC CI | observed AUC %.4f | 95%% CI [%.4f, %.4f] | valid bootstrap %d",
        observed_auc,
        ci_lo if ci_lo is not None else float("nan"),
        ci_hi if ci_hi is not None else float("nan"),
        len(boot_aucs),
    )

    return result


def permutation(
    X_dev: pd.DataFrame,
    y_dev: pd.Series,
    donors_dev: pd.Series,
    final_signature: list[str],
    classifier_name: str,
    best_params: dict | None = None,
    logger: logging.Logger | None = None,
    n_permutations: int = 1000,
    n_folds: int = N_INNER_FOLDS,
    seed: int = SEED,
) -> dict:
    """Run a donor-stratified permutation test on the development set."""
    log = logger or logging.getLogger()

    X_panel = X_dev[final_signature].copy()

    cv = StratifiedGroupKFold(
        n_splits=n_folds,
        shuffle=True,
        random_state=seed,
    )

    def mean_cv_auc(y_vector: pd.Series) -> float:
        # Run one donor-stratified CV pass and return the mean AUC across folds.
        fold_aucs = []

        for train_idx, valid_idx in cv.split(X_panel, y_vector, groups=donors_dev):
            X_tr = X_panel.iloc[train_idx]
            X_va = X_panel.iloc[valid_idx]
            y_tr = y_vector.iloc[train_idx]
            y_va = y_vector.iloc[valid_idx]

            scaler = StandardScaler()
            X_tr_scaled = scaler.fit_transform(X_tr)
            X_va_scaled = scaler.transform(X_va)

            model = make_classifier(classifier_name, best_params, y_tr)
            model.fit(X_tr_scaled, y_tr)

            proba = model.predict_proba(X_va_scaled)[:, 1]
            fold_auc = float(roc_auc_score(y_va, proba))
            fold_aucs.append(fold_auc)

        return float(np.mean(fold_aucs))

    # Observed AUC using the real (non-shuffled) labels
    observed_auc = mean_cv_auc(y_dev)

    # One row per donor, so permutation shuffles at the donor level,
    # not the sample level - a donor's repeated samples must stay together
    donor_table = pd.DataFrame({
        "donor": donors_dev.astype(str).values,
        "y": y_dev.values,
    }).drop_duplicates(subset=["donor"])

    rng = np.random.RandomState(seed)
    null_aucs = []

    for i in range(n_permutations):
        shuffled_labels = rng.permutation(donor_table["y"].values)
        donor_to_label = dict(zip(donor_table["donor"], shuffled_labels))
        y_perm = donors_dev.astype(str).map(donor_to_label).astype(int)

        try:
            auc_i = mean_cv_auc(y_perm)
            null_aucs.append(float(auc_i))
        except Exception as e:
            log.warning("Permutation %d skipped: %s", i + 1, e)

    # +1 correction avoids a p-value of exactly zero
    p_value = (1 + sum(a >= observed_auc for a in null_aucs)) / (1 + len(null_aucs))

    result = {
        "observed_auc": float(observed_auc),
        "null_aucs": [float(x) for x in null_aucs],
        "p_value": float(p_value),
        "n_permutations": int(len(null_aucs)),
    }

    log.info(
        "Permutation test | observed AUC %.4f | p %.4f | successful permutations %d",
        result["observed_auc"],
        result["p_value"],
        result["n_permutations"],
    )

    return result

def prune_correlated_genes(
    X: pd.DataFrame,
    ranked_genes: list[str],
    correlation_threshold: float,
    logger: logging.Logger | None = None,
    label: str = "final_panel",
) -> Tuple[list[str], pd.DataFrame]:
    """
    Keep genes in ranked order, dropping later genes if they are too correlated
    with any already-kept earlier gene.

    Returns:
    - kept genes
    - the full Pearson correlation matrix over ranked_genes (computed here
      anyway to decide the pruning - returned so callers don't need to
      recompute the same matrix a second time)
    """
    log = logger or logging.getLogger()

    if len(ranked_genes) <= 1:
        return list(ranked_genes), X[ranked_genes].corr(method="pearson")

    corr_matrix = X[ranked_genes].corr(method="pearson")
    corr_abs = corr_matrix.abs()
    kept_genes: list[str] = []

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
        else:
            kept_genes.append(gene)

    log.info(
        "%s Pearson pruning kept %d of %d genes at |r| < %.2f",
        label,
        len(kept_genes),
        len(ranked_genes),
        correlation_threshold,
    )

    return kept_genes, corr_matrix

def finalise_panel(
    results_df: pd.DataFrame,
    X_dev: pd.DataFrame,
    y_dev: pd.Series,
    donors_dev: pd.Series,
    logger: logging.Logger | None = None,
    stability_n_iter: int = N_BOOTSTRAP,
    stability_max_iter_boruta: int = BORUTA_MAX_ITER,
    stability_fallback_top_n: int = FOLD_STABLE_FALLBACK_TOP_N,
    stability_min_stable: int = FOLD_MIN_STABLE_GENES,
    seed: int = SEED,
) -> dict:
    """
    Choose the winning FS+classifier, rebuild the final signature on the full
    development set, then apply full-development bootstrap Boruta stability
    and final Pearson redundancy pruning to produce the final panel.
    """
    log = logger or logging.getLogger()

    # 1. Winner from outer benchmark
    winner = pick_winner(results_df, logger=log)

    # 2. Refit winning FS on full dev - keep its own real selected genes
    final_signature_info = build_final_signature(
        X_dev=X_dev,
        y_dev=y_dev,
        donors_dev=donors_dev,
        fs_name=winner["fs_method"],
        logger=log,
        seed=seed,
    )

    signature = list(final_signature_info["signature"])
    full_ranking = list(final_signature_info["ranking"])

    # 3. Full-development bootstrap Boruta stability, restricted to the FS
    #    method's own selected genes only - not the entire development
    #    gene pool. Bootstrap and Pearson pruning below now only ever
    #    consider genes the winning method actually chose. Uses Boruta's
    #    own confirmed/rejected decision each resample, not a raw RF
    #    importance percentile cutoff - see run_bootstrap_boruta_stability.
    stable_genes, stability_freq = run_bootstrap_boruta_stability(
        X_train=X_dev[signature],
        y_train=y_dev,
        donors_train=donors_dev,
        fold_id="final_panel",
        logger=log,
        n_iter=stability_n_iter,
        stability_threshold=FINAL_STABILITY_THRESHOLD,
        fallback_top_n=stability_fallback_top_n,
        min_stable=stability_min_stable,
        max_iter_boruta=stability_max_iter_boruta,
        seed=seed,
    )

    # genes that passed the final stability rule, ordered by the FS
    # method's own full ranking (full_ranking here is used ONLY to decide
    # DISPLAY ORDER among the already-correctly-scoped stable genes -
    # every gene in stable_genes is already guaranteed to be inside
    # `signature`, since that's the only thing bootstrap was run on)
    high_freq_gene_set = set(stable_genes.index.tolist())
    high_freq_ranked_genes = [g for g in full_ranking if g in high_freq_gene_set]

    if not high_freq_ranked_genes:
        raise RuntimeError(
            "No genes survived the final stability threshold. "
            "Consider lowering FINAL_STABILITY_THRESHOLD."
        )

    # 4. Pearson redundancy pruning, preserving ranking
    final_panel, full_correlation_matrix = prune_correlated_genes(
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
        "correlation_matrix": full_correlation_matrix,
    }

