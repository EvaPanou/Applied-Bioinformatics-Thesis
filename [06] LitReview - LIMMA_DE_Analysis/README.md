# Purpose of this folder

The LitReview folder is the place to present the goal of my work, expand on what I have done and compare with how similar research has tackled analogous questions in the past

## Goal of my work
The specific objective of this project is to identify transcriptomic biomarkers of SLE in pregnant women through the application of machine-learning approaches to longitudinal microarray expression data from healthy and SLE individuals. The ultimate goal is to facilitate early and personalized monitoring of disease activity in individuals at risk of developing SLE, and attempting to utilise one of the many microarray datasets accumulated over past decades that remain an invaluable resource. 

## What I have done
To achieve the above objectives, the analytical pipeline includes:
* Data Preprocessing and Quality Control
  - Exploration of both metadata and gene-expression data.
  - Detection and removal of outliers.
  - Filtering of low-expressed to reduce noise.
* Differential Expression Analysis (DEA)
  - Application of the LIMMA linear modeling framework for the identification of differentially expressed genes (DEGs) across both condition (SLE vs. healthy) and time points.
      - Condition-based comparisons: SLE vs. Healthy at each gestational time point.
      - Temporal comparisons: sequential time-point contrasts within each condition. 
  - Application of appropriate multiple-testing correction to control false discovery rates.
  - Construction of union DEG heatmaps with hierarchical clustering and clear annotation of condition and time point.
* Feature Optimization and Dimensionality Reduction
  - Integration of all identified DEGs into a unified feature set.
  - Filtering of top 5% most variable genes to refine machine learning input features.
  - Implementation of Principal Component Analysis (PCA) to visualize data structure and assess sample separation (both before and after DEG anaylisis)

## Selection of 3 Similar Papers

The 3 Papers chosen for this are the following, which are also included, in the same numerical order, in pdf format within this folder: 
1.	**Chowdhury et al, 2023:** *Drought‑responsive genes in tomato: meta‑analysis of gene expression using machine learning* (tomato paper)
2.	**Mishra et al, 2022:** *Longitudinal multi‑omics analysis identifies early blood‑based predictors of anti‑TNF therapy response in inflammatory bowel disease*
3.	**Lefol et al, 2023:** *TiSA: TimeSeriesAnalysis––a pipeline for the analysisof longitudinal transcriptomics data*

## Literature Overview
Each paper applies a multi-step workflow for transcriptomic longitudinal data analysis, but they differ in key aspects such as the application of machine learning. 
Following data acquisition, the **processing and normalization steps** are universally performed; Mishra et al. and TiSA both use DESeq2 for normalization and differential expression analysis, whereas Chowdhury et al. employs edgeR after read alignment, reflecting some minor variations in analytical software selection tailored to each data type and study scale. A common methodological thread among the three is the time-series or longitudinal characterization of transcriptomes (mostly RNA, but TiSA is applicable to Microarrays as well). **Differential gene expression** and **functional enrichment analyses** are prominent features in all workflows, though their implementations vary: TiSA and Chowdhury provide detailed temporal clustering and integrate gene ontology studies, while Mishra et al. combines transcriptomic and methylomic data, leveraging co-expression networks and disease phenotype stratification for a multi-omics perspective. Distinctively, the application of **machine learning** is foregrounded in Mishra et al., that builds predictive models for anti-TNF therapy response; while Chowdhury et al. applies XGBoost classification to refine candidate drought-responsive genes and validate their predictive accuracy using biomarker networks and QTL overlap. 
### Focus on DEG
Firstly, Mishra et al. combine traditional pairwise differential expression analysis against baseline (using DESeq2) with a comple,emtary longitudinal modeling (ImpulseDE2), by rigorously defining true DEGs as those significant and present in both models (FDR < 0.05); this is followed by GO enrichment on DEGs and co-expression modules, contrasting treatment arms (anti-TNF vs. anti-integrin) and culminating in machine learning–based feature selection for biomarker discovery through random forests. 
The TiSA pipeline complements this by performing differential expression via DESeq2 or limma, but explicitly contrasts both experimental conditions at every time point (conditional DEG analysis) and evaluates temporal changes within each group (temporal DEG analysis); it sets significance thresholds and leverages advanced visualization (heatmaps) and clustering (PART method) to distill both consistent and dynamic gene responses within the time course, thus offering a powerful template for multi-axis DE analysis and trajectory grouping. 
Lastly, Chowdhury et al. does something similar to the above, by emphasizing robust timepoint-specific DEG identification with edgeR (tight FDR and fold change cutoffs), and they also built a heatmap to interpret the longitudinal change of DEGs per timepoint. They assembled a global, non-redundant set of DEGs from across all timepoints, and filtered by expression variance (>95%) in order to produce a high-confidence feature matrix suitable for machine learning, thereby rigorously capturing both immediate and delayed responses but without explicit continuous modeling. 
### CONCLUSION
Based on these approaches, a strong course of action for a longitudinal transcriptomic project would be to implement a dual DE strategy — combining pairwise/conditional and longitudinal/temporal modeling (Mishra et al., TiSA), visualizing {and clustering expression patterns}, and finally constructing and variance-filtering a comprehensive set of DEGs (Chowdhury et al.) for downstream machine learning and functional analysis; this integrated strategy maximizes sensitivity to dynamic biological signals while supporting robust, interpretable biomarker selection.​

## NEXT STEPS
* Machine Learning Model Development
  - Splitting the dataset into training and testing subsets.
  - Training and comparison of multiple ML classifiers to predict disease status from gene-expression patterns.
* Model Evaluation and Validation
  - Assessment of model performance using quantitative metrics (e.g., accuracy, precision, recall, F1-score).
  - Execution of a Negative Control for Accuracy to verify robustness and ensure the generalizability of results.
