# Purpose of this folder

The LitReview folder is the place to present the goal of my work, expand on what I have done and compare with how similar research has tackled analogous questions in the past

## Goal of my work
The specific objective of this project is to identify transcriptomic biomarkers of SLE in pregnant women through the application of machine-learning approaches to longitudinal microarray expression data from healthy and SLE individuals. The ultimate goal is to facilitate early and personalized monitoring of disease activity in individuals at risk of developing SLE, and attempting to utilise one of the many microarray datasets accumulated over past decades that remain an invaluable resource. 

## What I have done
To achieve the above objectives, the analytical pipeline includes:
* Data Preprocessing and Quality Control
  - Exploration of both metadata and gene-expression data.
  - Detection and removal of outliers.
  - Filtering of low-expressed and low-variance genes to reduce noise.
* Differential Expression Analysis (DEA)
  - Identification of differentially expressed genes (DEGs) across both condition (SLE vs. healthy) and time points.
  - Application of appropriate multiple-testing correction to control false discovery rates.
* Feature Optimization and Dimensionality Reduction
  - Integration of all identified DEGs into a unified feature set.
  - Implementation of Principal Component Analysis (PCA) to visualize data structure and assess sample separation.
* Machine Learning Model Development
  - Splitting the dataset into training and testing subsets.
  - Training and comparison of multiple ML classifiers to predict disease status from gene-expression patterns.
* Model Evaluation and Validation
  - Assessment of model performance using quantitative metrics (e.g., accuracy, precision, recall, F1-score).
  - Execution of a Negative Control for Accuracy to verify robustness and ensure the generalizability of results.

## Literature Review of 3 Similar Papers

The 3 Papers chosen for this are the: 
1.	**Chowdhury et al, 2023:** *Drought‑responsive genes in tomato: meta‑analysis of gene expression using machine learning*
2.	**Mishra et al, 2022:** *Longitudinal multi‑omics analysis identifies early blood‑based predictors of anti‑TNF therapy response in inflammatory bowel disease*
3.	**Lefol et al, 2023:** *TiSA: TimeSeriesAnalysis––a pipeline for the analysisof longitudinal transcriptomics data*
