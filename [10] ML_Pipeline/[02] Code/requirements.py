
# -------------------------------------------------------------------------------------------------------------------
# This script contains all the commands for python packages that you will need to have installed in your computer
# and then you will need to have them imported for the main pipeline script to run.


# IT WILL NEED TO BE LOCATED IN THE SAME DIRECTORY AS THE MAIN PIPELINE SCRIPT "sle_pipeline.py"
# in order to be directly called into action
# ------------------------------------------------------------------------------------------------------------------

# INSTALLATIONS
# -------------------------------------
# We first try to import each package; if that fails, we pip-install it into the same
# Python interpreter that is running this script (sys.executable guarantees the right
# environment) and then continue.
# We exclude packages that ship with Python (json, os, sys, subprocess, importlib, ...).

# Format: (pip_name, import_name)
#   - pip_name    -> the name on PyPI, e.g. "scikit-learn"
#   - import_name -> the name you `import` in Python, e.g. "sklearn"

import importlib
import subprocess
import sys


_REQUIRED = [
    ("numpy",           "numpy"),
    ("pandas",          "pandas"),
    ("scikit-learn",    "sklearn"),
    ("matplotlib",      "matplotlib"),
    ("shap",            "shap"),
    ("mrmr-selection",  "mrmr"),
    ("kneed",           "kneed"),
    ("xgboost",         "xgboost"),
    ("Boruta",          "boruta"),
]


for pip_name, import_name in _REQUIRED:
    try:
        importlib.import_module(import_name)
    except ImportError:
        print(f"[requirements] '{import_name}' not found -- installing '{pip_name}' ...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])
        except subprocess.CalledProcessError:
            if pip_name == "shap":
                print(f"[requirements] '{pip_name}' install failed -- continuing without it")
            else:
                raise



# IMPORTS
# -------------------------------------
# standard library
import json
import logging
import os
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple



# matplotlib
import matplotlib
matplotlib.use("Agg")                 # safe on any laptop, no display required
import matplotlib.pyplot as plt



# others
import numpy as np
import pandas as pd
from boruta import BorutaPy
from xgboost import XGBClassifier


# scikit-learn
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV         # used by SVM-RFE
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import (
    accuracy_score, confusion_matrix, roc_auc_score, roc_curve,
)
from sklearn.model_selection import (
    GridSearchCV, StratifiedGroupKFold, StratifiedShuffleSplit,
)
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


# optional: SHAP 
# Renamed _HAS_SHAP -> HAS_SHAP so it gets exported by `from requirements import *`
# (Python's `import *` skips names starting with an underscore).
try:
    import shap
    HAS_SHAP = True
    # SHAP's kernel explainer prints per-sample arithmetic at INFO,
    # which floods our log file -- silence its noise.
    logging.getLogger("shap").setLevel(logging.WARNING)
except Exception:
    shap = None                       # type: ignore
    HAS_SHAP = False