# Applied-Bioinformatics-Thesis

> *This repository will contain all the relevant steps and commands to follow in order to identify biomarkers related to auto-immune diseases, using Machine Learning methods*

## Info about the dataset

### The data used in this repository is from the research paper *"Longitudinal profiling of human blood transcriptome in healthy and lupus pregnancy"*.

### The raw microarray dataset described in the manuscript is deposited in the NCBI GEO (accession no. GSE108497), along with supplementary files like series matrices. It can be also be found processed in the Autoimmune Diseases Explorer (ADEx), with the same GEO Accession number, along with the processed metadata file. Here the processed ADEx version was used to push the analysis forward. However the processed ADEx metadata file was further modified to include more info, during the first steps of the analysis.

### Useful links

- Research Paper Source: [Longitudinal profiling of human blood transcriptome in healthy and lupus pregnancy](https://doi.org/10.1084/jem.20190185)
- GEO Accession Display: [Series GSE108497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108497)
- Autoimmune Diseases Explorer: [ADEx](https://adex.genyo.es/)

## Goal of this project
The aim of this project is to identify transcriptomic biomarkers of the autoimmune disease systemic lupus erythematosus (SLE) in pregnant women, by applying machine-learning methods to longitudinal microarray expression data from healthy and SLE individuals. The ultimate goal is to enable early, personalized monitoring of disease activity in individuals with suspected SLE. By preprocessing and further analysis of the data, high-dimensional gene-expression data will be transformed into optimized feature sets, before splitting into training and test sets. Then, training and comparison of multiple Machine Learning classifiers will take place, to predict disease status from gene-expression patterns. The final step, is to carry out comparisons across algorithms and feature sets to identify the best-performing model, and run a Negative Control for Accuracy to confirm its robustness and applicability.

## What have I done so far

### Completed Tasks

1. Metadata Assembly
2. Exploratory Metadata Analysis
3. Exploratory Raw Data Analysis
4. Raw Data Preprocessing
5. Initial PCA
6. Post-PCA Analysis

### Upcoming Tasks
7. DE & Further Processing
8. Proof-of-Concept - ML using DE-selected genes as features
9. Train-Test Split
10. Feature Engineering for ML
11. ML Experiments with Top-10 Genes via:
- DE
- PCA
- Univariate Analysis
12. Model Evaluation & Visualization
13. External Test Validation
14. Negative Control Analysis

