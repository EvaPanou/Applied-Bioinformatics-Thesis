
"""
=============================================================================
requirements.py
-----------------------------------------------------------------------------
WHAT THIS FILE DOES:
  It imports every package and tool the rest of the pipeline needs. 
  The rest of the scripts then write  `from requirements import *` so they
  do not have to repeat all of these import lines.

BEFORE YOU RUN THE PIPELINE (one time only):
  Install the packages by opening a terminal in this folder and running:

      pip install -r requirements.txt 
      # txt file within this directory, for mass installation of required packages

WHERE THIS FILE MUST BE SAVED:
  In the same directory as sle_pipeline.py, otherwise the imports will not work.
=============================================================================

"""
# -----------------------------------------------------------------------------
# Standard-library tools (built into Python, nothing to install)
# -----------------------------------------------------------------------------
import json                          # read and write the .json output files
import logging                       # write progress messages to the log file
import os                            # talk to the operating system (e.g. set the seed)
import random                        # Python's built-in random-number generator
from collections import Counter      # quick way to count how many of each label
from dataclasses import dataclass    # tidy way to bundle related values into one object
from pathlib import Path             # handle file paths safely on any operating system
from typing import Iterable, Tuple   # type hints describing what a function takes/returns


# -----------------------------------------------------------------------------
# Plotting 
# -----------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Scientific and ML packages
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
from boruta import BorutaPy
from xgboost import XGBClassifier


# scikit-learn: for ML
from sklearn.ensemble import RandomForestClassifier   # Random Forest model
from sklearn.feature_selection import RFECV            # used by the SVM-RFE selector
from sklearn.linear_model import (
    LogisticRegression,      # plain logistic regression
    LogisticRegressionCV,    # logistic regression that tunes itself with cross-validation
)
from sklearn.metrics import (
    accuracy_score,     # fraction of predictions that were correct
    confusion_matrix,   # table of correct vs incorrect predictions per class
    f1_score,           # the f1 score, which balances precision and recall into one number
    roc_auc_score,      # the AUC score (the main metric used in this pipeline)
    roc_curve,          # the coordinates needed to draw an ROC curve
)
from sklearn.model_selection import (
    GridSearchCV,            # tries every combination of settings to find the best
    StratifiedGroupKFold,    # cross-validation that keeps each donor in ONE fold only
    StratifiedShuffleSplit,  # used for the donor-level train/test split
)
from sklearn.naive_bayes import GaussianNB         # Naive Bayes classifier
from sklearn.preprocessing import StandardScaler   # rescales features onto a common scale
from sklearn.svm import SVC                        # Support Vector Machine classifier


# -----------------------------------------------------------------------------
# SHAP Plotting
# -----------------------------------------------------------------------------
# If this raises an error and it is not installed, the pipeline simply skips the SHAP plots.
#
# The flag is HAS_SHAP with NO leading underscore on purpose: `import *` skips
# names starting with "_", so an underscore would hide it from the other scripts.
try:
    import shap
    HAS_SHAP = True
    logging.getLogger("shap").setLevel(logging.WARNING)  # stop SHAP flooding the log
except Exception:
    shap = None
    HAS_SHAP = False