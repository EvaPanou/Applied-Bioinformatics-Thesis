#-----------------------------------------------------------------------------------------
# This script contains all the commands for configured parameters, paths, or column names
# that you will need to have imported for the main pipeline script to run.

# IT WILL NEED TO BE LOCATED IN THE SAME DIRECTORY AS THE MAIN PIPELINE SCRIPT "sle_pipeline.py"
# in order to be direclty called into action

# Defaults are explained in "METHODS.md". Folder Layout is explain in "folder_layout.md". 

# Resolve paths are relative to wherever the script lives, so the same configurtion

# Configurations
# ------------------------------------------------------------------------------------------

# (stdlib Path is for filesystem paths; everything else lives in requirements.py)
from pathlib import Path
from typing import Tuple

# Paths 
# --------
# Path() handles the spaces and brackets in the folder names natively.
CODE_DIR    = Path(__file__).resolve().parent       # CODE_DIR is the folder this file lives in
PROJECT_DIR = CODE_DIR.parent                       # PROJECT_DIR is one level up 
INPUT_DIR   = PROJECT_DIR / "[01] Input"
OUTPUT_DIR  = PROJECT_DIR / "[03] Output"

# Colour palette definition
# --------------------------
CB_BLUE   = "#0072B2"             # for the SLE class / primary curve
CB_ORANGE = "#E69F00"             # for healthy / comparator
CB_GREY   = "#999999"             # Grey is for null distributions / context


# File names/paths 
# ------------------
DATA_FILENAME    = "ML_FeatureMatrix_AllDEGs_with_Metadata.tsv"
LOG_FILENAME     = "pipeline_log.txt"

DATA_PATH        = INPUT_DIR  / DATA_FILENAME
TABLES_DIR       = OUTPUT_DIR / "tables"
FIGURES_DIR      = OUTPUT_DIR / "figures"
LOG_PATH         = OUTPUT_DIR / LOG_FILENAME


# Metadata columns
# ------------------
SAMPLE_ID_COL        = "SampleID"
DONOR_ID_COL         = "DonorID"
LABEL_COL            = "Condition"        
LABEL_TEXT_COL       = "Condition_label"
TIMEPOINT_COL        = "Timepoint"
TIMEPOINT_LABEL_COL  = "Time_label"


METADATA_COLUMNS: Tuple[str, ...] = (                       # for helper_functions.py
    SAMPLE_ID_COL,
    DONOR_ID_COL,
    LABEL_COL,
    LABEL_TEXT_COL,
    TIMEPOINT_COL,
    TIMEPOINT_LABEL_COL,
)


# Reproducibility 
# -----------------------
SEED = 42


# Splits 
# ------------------------
HOLDOUT_FRAC    = 0.20
N_OUTER_FOLDS   = 5
N_INNER_FOLDS   = 5
N_REPEATS       = 1


# Feature selectors 
# --------------------
FS_METHODS = ("mrmr", "lasso", "elasticnet", "svm_rfe", "rf", "boruta")   # to feed data into the benchmarking

MRMR_TOP_K = 30      # also the per-method top-k cap and bootstrap top-k
LASSO_C_VALUES = (0.001, 0.01, 0.1, 1.0, 10.0)
ELASTICNET_C_VALUES = (0.001, 0.01, 0.1, 1.0, 10.0)
ELASTICNET_L1_RATIOS = (0.2, 0.5, 0.8)
RF_IMPORTANCE_TREES = 300
SVM_RFE_STEP = 1
SVM_RFE_MIN_FEATURES = 2
SVM_RFE_INNER_FOLDS = 3
BORUTA_MAX_ITER = 100       


# Classifier benchmark 
# --------------------------------
# Hyperparameter grids

CLASSIFIERS = ("svm", "logreg", "rf", "xgboost", "nb", "knn")

LR_C_VALUES         = (0.001, 0.01, 0.1, 1.0, 10.0)
SVM_C_VALUES        = (0.1, 1.0, 10.0)
SVM_GAMMA_VALUES    = ("scale", "auto")
RF_DEPTH_VALUES     = (None, 5, 10)
RF_LEAF_VALUES      = (1, 2, 5)
XGB_DEPTH_VALUES    = (3, 5, 7)
XGB_LR_VALUES       = (0.01, 0.05, 0.1)
XGB_N_ESTIMATORS    = (100, 200, 300)
NB_VAR_SMOOTHING    = (1e-9, 1e-8, 1e-7)
KNN_NEIGHBORS       = (3, 5, 7, 9)

# Per-fold bootstrap RF stability filter 
# ---------------------------------------------
FOLD_BOOTSTRAP_ITERATIONS = 400
FOLD_BOOTSTRAP_TREES = 200    # meaning RF trees
FOLD_BOOTSTRAP_TOP_PERCENTILE = 20    # keep top 80% by RF importance per iter
FOLD_STABILITY_THRESHOLD = 0.70   # gene must be in top-80% in >=70% of iterations
FOLD_STABLE_FALLBACK_TOP_N = 20     # if too few pass, take top-20 by frequency
FOLD_MIN_STABLE_GENES = 5      # below this, use the fallback


# Final-panel bootstrap stability column 
# -------------------------------------------
N_BOOTSTRAP = 200
FINAL_STABILITY_THRESHOLD = 0.80
FINAL_PEARSON_THRESHOLD = 0.90


# Held-out evaluation 
# -------------------------------------
N_PERMUTATIONS         = 200
N_HOLDOUT_BOOTSTRAP    = 200