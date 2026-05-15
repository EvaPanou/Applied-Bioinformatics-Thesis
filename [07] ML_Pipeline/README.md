# Donor-Stratified SLE Diagnostic ML Pipeline

This directory provides a visual-methodology overview of the donor-aware machine-learning pipeline developed for the identification and evaluation of a compact diagnostic gene signature for systemic lupus erythematosus (SLE).

The purpose of these flowcharts is to summarize the main analytical logic of the pipeline without replacing the detailed methods and results documentation. The figures show how the workflow progresses from the prepared gene-expression input matrix, through donor-aware model development and feature-selection benchmarking, to final gene-panel definition and sealed held-out evaluation.

A more detailed methodology/results README, including references, full output interpretation, performance discussion, and selected result plots, will be provided separately within the "[03] Output" folder.

---

## Folder Structure

```text
.
├── README.md
├── [01] Input/ = includes the singluar matrix file that was used as input in this pipeline
├── [02] Code/ = includes five Python script files that built & configured this ML pipeline
├── [03] Ouput/ = includes all output tables and figures, as well as result metrics, to be used for the discussion
└── [04] Flowchart/ = contains overall & sectional flowcharts for a better presentation of the methodology
```

---

## Starting Data

The pipeline starts from a prepared machine-learning input dataset containing differentially expressed gene-derived features and sample metadata.

The input data used in the visualized workflow include:

- **507 samples**
- **158 donors**
- repeated samples across multiple timepoints (max: 5)
- gene-expression features used as candidate predictors
- metadata including donor identity and disease condition (SLE vs Healthy)

Because multiple samples may originate from the same donor, all splitting and validation steps are performed using a donor-aware strategy. This prevents samples from the same biological individual from appearing in both training and testing subsets.

---

## Overall Pipeline Flowchart

The complete pipeline is summarized in the PDF file below:

[Open Overall Pipeline Flowchart](<./[04] Flowchart/Overall_Pipeline_Flowchart.pdf>)

The overall flowchart presents the full methodology as four connected sections:

1. **Input Data and Donor Split**
2. **Algorithm Benchmarking on the Development Set**
3. **Final Gene-Signature Panel Definition**
4. **Held-Out Test Set Evaluation**

The main flow of the pipeline is:

```text
Prepared DEG + metadata input file
        ↓
Donor-stratified development (train) / held-out (test) split
        ↓
Development-set feature-selector/classifier benchmarking
        ↓
Final compact gene-signature panel derivation
        ↓
Sealed held-out test-set evaluation
        ↓
Validated SLE diagnostic gene signature & Feature Selection / Classifier Combination

```
## How the Sections Connect

The four section-specific diagrams should be read as consecutive stages of a single workflow.

Section 1 defines the two major data branches: the development cohort and the sealed held-out cohort. Section 2 uses only the development cohort to compare feature-selection and classification strategies. Section 3 uses the selected winning strategy and the full development cohort to derive the final compact gene panel. Section 4 reconnects the sealed held-out cohort for the first time and uses it for final evaluation only.

This structure ensures that model selection and final evaluation remain separated.

---

## Section 1 — Input Data and Donor Split

Section 1 explains how the prepared input matrix is loaded and divided into development and held-out subsets at donor level.

![Section 1 — Input Data and Donor Split](<./[04] Flowchart/Flowchart_Section1.png>)

The pipeline begins with the prepared DEG and metadata input file. The dataset contains **507 samples from 158 donors**. A donor-level stratified 80/20 split is then performed, producing:

- **Development cohort:** 126 donors / 404 samples
- **Held-out cohort:** 32 donors / 103 samples

The development cohort is used for feature selection, model benchmarking, and final panel construction. The held-out cohort is sealed after the initial split and is not used during model selection.

The main principle of this section is that the same donor must stay in only one subset.

**What this figure shows:**

- the starting input dataset;
- package/configuration setup;
- donor-level stratified splitting;
- separation of development and held-out cohorts;
- the leakage-prevention rule that each donor remains in only one subset;
- generation of the split summary file.

---

## Section 2 — Algorithm Benchmarking on the Development Set

Section 2 describes how the development cohort is used to compare feature-selection and classification strategies.

![Section 2 — Algorithm Benchmarking on the Development Set](<./[04] Flowchart/Flowchart_Section2.png>)

Only the development cohort enters this stage. The workflow uses donor-aware outer 5-fold cross-validation, where each outer fold separates training and testing donors without donor overlap.

Within each outer training split, bootstrap Random Forest stability filtering is applied to identify genes that are repeatedly prioritized across resampling iterations. The stabilized feature space is then used for feature-selection and classifier benchmarking.

The feature-selection methods include:

- mRMR
- LASSO
- Elastic Net
- SVM-RFE
- Random Forest importance
- Boruta

The classifiers include:

- Logistic Regression
- Naive Bayes
- Support Vector Machine
- k-nearest neighbours
- Random Forest
- XGBoost

Inner donor-aware 5-fold cross-validation is used for hyperparameter tuning. After tuning, each model is refitted on the outer training donors and evaluated on the outer testing donors. The outer-fold ROC AUC values are collected for every feature-selector/classifier pair, and the best-performing combination is selected based on mean outer-fold AUC.

**What this figure shows:**

- development-set-only benchmarking;
- donor-aware outer cross-validation;
- bootstrap RF stability filtering;
- multiple feature-selection methods;
- multiple classifier families;
- inner cross-validation for hyperparameter tuning;
- outer-fold AUC collection;
- selection of the best feature-selector/classifier combination.

---

## Section 3 — Final Gene-Signature Panel Definition

Section 3 explains how the final compact gene panel is derived after the best-performing modelling strategy has been selected.

![Section 3 — Final Gene-Signature Panel Definition](<./[04] Flowchart/Flowchart_Section3.png>)

The winning feature selector is rerun on the full development cohort to generate a global ranking of candidate genes. A full-development bootstrap Random Forest stability analysis is then performed to estimate per-gene selection frequency.

Genes with high bootstrap stability are retained as candidate panel genes. In the visualized workflow, genes with selection frequency greater than or equal to **0.80** are kept for further consideration.

An elbow curve is inspected for candidate panel size, but it is not treated as the final decision criterion. The final panel is instead refined through stability selection and redundancy reduction.

Pearson correlation pruning is used to remove highly correlated genes. Genes are pruned when correlation exceeds the redundancy threshold:

```text
|r| ≥ 0.90
```

This produces a compact, non-redundant final gene panel. In the current flowchart, the final panel contains **6 genes**.

**What this figure shows:**

- selection of the winning feature-selector/classifier strategy;
- rerunning the winning selector on the full development cohort;
- global gene ranking;
- full-development bootstrap RF stability analysis;
- high-frequency gene retention;
- elbow-curve inspection;
- Pearson correlation-based redundancy pruning;
- final compact gene-signature panel definition.

---

## Section 4 — Held-Out Test Set Evaluation

Section 4 describes the final evaluation stage using the sealed held-out donor cohort.

![Section 4 — Held-Out Test Set Evaluation](<./[04] Flowchart/Flowchart_Section4.png>)

After the final panel is defined, the winning classifier is trained on the full development cohort using only the selected final-panel genes. The trained model is then applied once to the sealed held-out donors.

The held-out evaluation includes:

- ROC AUC
- sensitivity
- specificity
- accuracy
- F1-score
- confusion matrix

To assess uncertainty, the workflow includes bootstrap confidence interval estimation for held-out AUC using **200 resamples**. A donor-stratified permutation test with **200 permutations** is also performed on the development set to generate a null-performance reference.

A logistic-regression comparator is trained and evaluated on the same final gene panel as a transparent baseline model. SHAP-based interpretation is then applied to the held-out donors to summarize the contribution of the selected genes to model predictions.

**What this figure shows:**

- final training on the full development cohort;
- one-time evaluation on sealed held-out donors;
- held-out performance metric calculation;
- bootstrap confidence interval estimation;
- donor-stratified permutation testing;
- logistic-regression comparator evaluation;
- SHAP-based model interpretation;
- validated SLE diagnostic gene signature.
