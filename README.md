# Applied-Bioinformatics-Thesis

**Identifying biomarkers related to auto-immune diseases using Machine Learning methods**

Master's Thesis — Applied Bioinformatics, Aristotle University of Thessaloniki
Eva Panou

> This repository includes the full pipeline for identifying a useful blood-transcriptomic gene signature for **Systemic Lupus Erythematosus (SLE)**

---

## 1. Research Question & Hypothesis

**Main hypothesis:** a small, stable set of differentially expressed genes, identified through a donor-aware and stability-based machine-learning approach applied to blood transcriptomic data, is sufficient to discriminate SLE patients from healthy controls with high, generalisable diagnostic accuracy on unseen donors.

SLE is a chronic, multisystem autoimmune disease marked by substantial heterogeneity in presentation. At the transcriptomic level, its most consistent hallmark is overexpression of interferon-inducible genes in peripheral blood (the "interferon signature"), reported across many independent cohorts. This thesis asks whether that signal — and complementary non-interferon signal — can be distilled into a small, reproducible gene panel that:

1. is derived through a methodologically rigorous, leakage-free pipeline (donor-aware splitting throughout, since the dataset is longitudinal with multiple samples per donor);
2. generalises to a sealed, never-touched-during-development held-out cohort;
3. is not simply tracking pregnancy status or pregnancy complications rather than SLE itself, given the dataset's pregnancy-cohort design (see §6).

---

## 2. Data Source

The dataset is **GSE108497**, from Hong et al., *"Longitudinal profiling of human blood transcriptome in healthy and lupus pregnancy,"* *J Exp Med* 2019 ([doi.org/10.1084/jem.20190185](https://doi.org/10.1084/jem.20190185)). It is a longitudinal, multicenter microarray study from the PROMISSE cohort: whole-blood transcriptomes from SLE-pregnant, healthy-pregnant, SLE-non-pregnant, and healthy-non-pregnant women, sampled at up to five timepoints per donor (four pregnancy windows plus postpartum), profiled on Illumina HT-12 V4 beadchips.

Rather than reprocessing the raw microarray data from scratch, this project sources it through **[ADEx (Autoimmune Diseases Explorer)](https://adex.genyo.es/)** — a curated database that reprocesses public autoimmune-disease omics datasets through a single, homogeneous pipeline (Martorell-Marugán et al., *BMC Bioinformatics* 2021, [doi.org/10.1186/s12859-021-04268-4](https://doi.org/10.1186/s12859-021-04268-4)).

Working from the ADEx-processed version (rather than raw GEO data) means normalization and probe-to-gene mapping are already handled upstream in a standardized way — but the ADEx metadata still required substantial manual correction for this thesis.

### Useful links

- Source paper: [Hong et al. 2019, *J Exp Med*](https://doi.org/10.1084/jem.20190185)
- GEO accession: [GSE108497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108497)
- ADEx database: [adex.genyo.es](https://adex.genyo.es/)
- ADEx paper: [Martorell-Marugán et al. 2021, *BMC Bioinformatics*](https://doi.org/10.1186/s12859-021-04268-4)

---

## 3. Repository Structure

| Folder | Description | README(s) |
|---|---|---|
| `00_Metadata/` | Metadata assembly & exploratory data analysis | [README](./00_Metadata/README.md) · [04_Analysis_Results](./00_Metadata/04_Analysis_Results/README.md) |
| `01_RawData_&_PCA/` | Raw expression QC/outlier removal + ALASCA longitudinal modeling | [README](./01_RawData_%26_PCA/README.md) · [02_Analysis_Results](./01_RawData_%26_PCA/02_Analysis_Results/README.md) · [04_ALASCA_Output](./01_RawData_%26_PCA/04_ALASCA_Output/README.md) |
| `02_DGE_Analysis/` | Differential gene expression (limma pipeline) | [README](./02_DGE_Analysis/README.md) · [02_Results](./02_DGE_Analysis/02_Results/README.md) · [A_SLE_vs_Healthy](./02_DGE_Analysis/02_Results/A_SLE_vs_Healthy/README.md) · [B_Pooled_SLE_vs_Healthy](./02_DGE_Analysis/02_Results/B_Pooled_SLE_vs_Healthy/README.md) · [C_Dream_SLE_vs_Healthy](./02_DGE_Analysis/02_Results/C_Dream_SLE_vs_Healthy/README.md) · [Gene_Clusters](./02_DGE_Analysis/02_Results/Gene_Clusters/README.md) · [DE_Genes](./02_DGE_Analysis/02_Results/DE_Genes/README.md) |
| `03_Feature_Validation/` | Validation of the DEG feature set (union/intersection methods) | [README](./03_Feature_Validation/README.md) · [02_Results](./03_Feature_Validation/02_Results/README.md) · [MA_Plots](./03_Feature_Validation/02_Results/MA_Plots/README.md) · [PCA_Plots](./03_Feature_Validation/02_Results/PCA_Plots/README.md) · [PERMANOVA](./03_Feature_Validation/02_Results/PERMANOVA/README.md) · [Volcano_Plots](./03_Feature_Validation/02_Results/Volcano_Plots/README.md) |
| `04_ML_Prerequisites/` | ML-ready input table construction | [README](./04_ML_Prerequisites/README.md) |
| `05_ML_Pipeline/` | Donor-aware benchmarking, panel derivation, held-out evaluation | [README](./05_ML_Pipeline/README.md) · [01_Code](./05_ML_Pipeline/01_Code/README.md) · [02_Output](./05_ML_Pipeline/02_Output/README.md) |
| `06_Final_Panel_Validation/` | Confounder & specificity validation of the final 3-gene panel | [README](./06_Final_Panel_Validation/README.md) · [02_Validation_Output](./06_Final_Panel_Validation/02_Validation_Output/README.md) |
| `.gitattributes` / `.gitignore` | Git LFS tracking rules and ignore patterns | — |
| `LICENSE` | Repository license | — |
| `README.md` | This file | — |

Each stage's subdirectory has files and folders that follow a shared convention:

- **File-numbering convention:** within each folder, files are numbered `00_`, `01_`, `02_`, … in the order the pipeline touches them. 
    -- `00_` files are **inputs** carried in from the previous stage; 
    -- odd numbers from `01_` onwards are **code** (notebooks or scripts); 
    -- even numbers from `02_` onwards are **outputs** (results, tables, figures) produced by that code.


## 4. Methodology Overview

The pipeline runs as seven sequential stages, each stage's output becoming the next stage's input:

1. **`00_Metadata`** — GEO/ADEx metadata is assembled, cleaned, and explored: correcting never-pregnant (NP) donor mislabeling, deduplicating at the donor level, fixing unit inconsistencies in gestational-age fields, and resolving within-donor inconsistencies by majority vote. NP donors are held out throughout as a specificity check rather than discarded.
2. **`01_RawData_&_PCA`** — Outlier detection on the raw expression matrix (IQR-based and PCA + Mahalanobis/chi-squared methods), plus ALASCA longitudinal modeling of expression trajectories over pregnancy timepoints (NP donors excluded, since they lack a genuine repeated time series).
3. **`02_DGE_Analysis`** — Differential expression analysis via `limma`, using a pooled random-effects meta-analysis approach (`metafor::rma()`) rather than naive per-timepoint pooling, to identify genes significantly different between SLE and healthy donors.
4. **`03_Feature_Validation`** — Validation of the DEG feature set produced in stage 2, comparing selection methods (union vs. intersection across DE approaches) via PCA, PERMANOVA, and volcano/MA plots.
5. **`04_ML_Prerequisites`** — Merges the validated 44-gene `Method_All_Union` expression panel with corrected sample metadata into a single ML-ready input table (pure formatting — no modeling).
6. **`05_ML_Pipeline`** — The core modeling stage: donor-stratified development/held-out splitting, nested cross-validation benchmarking of feature-selection methods against classifiers, bootstrap-stability-based final panel derivation, and sealed held-out evaluation with SHAP interpretation.
7. **`06_Final_Panel_Validation`** — Confirms the final 3-gene panel tracks SLE status specifically — not pregnancy complication severity, pregnancy status alone, or batch — using likelihood-ratio tests against the held-out never-pregnant donor cohort.

A visual flow diagram tying all seven repository stages together is planned but not yet built; this section will link to it once available. (A separate, narrower flowchart set already exists inside `05_ML_Pipeline` — a 4-section diagram covering just that stage's donor-split → benchmarking → panel-derivation → held-out-evaluation workflow — see that stage's own README.)

---

## 5. Methodological Safeguards

Several design decisions recur across stages, aimed at preventing leakage and overfitting given the dataset's longitudinal (repeated-measures) structure:

- **Donor-aware splitting** — all samples from a given donor stay in a single subset across every split and cross-validation fold, everywhere in the pipeline.
- **Sealed held-out evaluation** — the held-out cohort is untouched during feature selection, model tuning, or panel construction, and is used exactly once.
- **Nested cross-validation** — separate inner (hyperparameter tuning) and outer (performance estimation) loops.
- **Bootstrap stability filtering** — genes are prioritized by how consistently they're selected under donor-level resampling, not by a single run's ranking.
- **Permutation testing** — a null AUC distribution from label permutation is used to assess whether observed performance exceeds chance.
- **Specificity validation** — the final panel is explicitly tested against alternative explanations (pregnancy, complications, batch) rather than assumed to be SLE-specific by construction.

---

## 6. Main Results

**Current final panel:** **HERC5, BATF2, SPATS2L**, derived via Boruta feature selection + SVM classification on the 44-gene `Method_All_Union` candidate space (bootstrap selection frequency ≥0.80, Pearson pruning at |r|≥0.90), and confirmed in `06_Final_Panel_Validation` to track SLE status specifically rather than pregnancy status or complication severity.

**Sealed held-out performance** (27 donors / 96 samples, never touched during model or panel selection):

| Metric | Value |
|---|---|
| AUC | 0.861 (95% CI 0.774–0.921, 200 bootstrap resamples) |
| Sensitivity | 0.810 |
| Specificity | 0.788 |
| Accuracy | 0.802 |
| F1 | 0.843 |
| Permutation test | observed AUC 0.923 (development-set CV) vs. null mean 0.497, *p* = 0.005 (200 donor-stratified permutations) |

A logistic-regression comparator on the same 3-gene panel scored slightly *higher* on the held-out set (AUC 0.884, 95% CI 0.803–0.939) than the winning SVM — reported transparently as a baseline check rather than smoothed over; with only 3 features, a linear model is a reasonable competitor to SVM here. SHAP analysis on the held-out set ranks the panel's contribution as HERC5 > BATF2 > SPATS2L, consistent across the beeswarm and mean-|SHAP| plots.

---

## 8. Key References

### Data Source Paper
- Hong S, Banchereau R, Maslow B-SL, et al. Longitudinal profiling of human blood transcriptome in healthy and lupus pregnancy. *J Exp Med*. 2019;216(5):1154–1169. https://doi.org/10.1084/jem.20190185

### Other References
- Martorell-Marugán J, López-Domínguez R, García-Moreno A, et al. A comprehensive database for integrated analysis of omics data in autoimmune diseases. *BMC Bioinformatics*. 2021;22:343. https://doi.org/10.1186/s12859-021-04268-4
- Kursa MB, Rudnicki WR. Boruta – a system for feature selection. *Fundamenta Informaticae*. 2010;101(4):271–285.
- Lundberg SM, Lee S-I. A unified approach to interpreting model predictions. *NeurIPS*. 2017;30:4765–4774.
- Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. *Patterns*. 2023;4(9):100804.
