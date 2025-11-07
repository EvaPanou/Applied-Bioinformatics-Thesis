# DE Genes lists - ML Feature Extraction

This section describes the feature extraction phase that bridges the differential expression (DE) analysis with the forthcoming machine learning (ML) modeling. All tab-delimited files in this folder were derived from the **limma_DEG_pipeline.R**, applied to the SLE pregnancy dataset to identify biologically meaningful gene expression signatures distinguishing **SLE vs. Healthy** conditions across gestation. 

## 1. Methodological Context

Following differential gene expression analysis, gene-level statistics (logFC, adjusted p-values) were used to identify **significant and biologically informative gene sets**, as a unified gene set across conditions and timepoints. Genes were filtered using thresholds of *FDR < 0.05* and *|logFC| ≥ 1*.  
Two complementary selection strategies were then applied:

1. **Biological relevance:** retention of all DEGs representing key SLE-associated pathways (interferon response, neutrophil activation, complement regulation).  
2. **Statistical informativeness:** selection of genes within the **top 5% of variance** across samples, using `matrixStats::rowVars`.  

The resulting subsets were combined into hierarchical feature tiers for subsequent machine learning tasks, balancing interpretability, biological signal strength, and numerical stability.

---

## 2. File Index

| File Name | Description | Type | Intended Use |
|------------|--------------|------|---------------|
| **AllDEGs.tsv** | Unique list of all DEGs (union of pairwise and temporal analyses). Represents the biologically grounded baseline feature set. | Gene list | Input for broad biological modeling and feature overlap analysis. |
| **HighVariance_genes_top5pct.tsv** | Top 5% of genes by variance across samples. Captures expression variability independent of condition. | Gene list | Used for variance-based feature filtering to enhance discriminative power. |
| **DEG_HighVariance.tsv** | Intersection of DEGs and high-variance genes. Defines the compact, biologically enriched core set. | Gene list | Serves as the principal gene set for feature matrix construction and interpretation. |
| **ML_FeatureMatrix_AllDEGs_log2.tsv** | Expression matrix of all DEGs (log2-transformed). Rows = samples; columns = genes. | Feature matrix | Suitable for exploratory unsupervised analysis (PCA, clustering). |
| **ML_FeatureMatrix_DEG_HighVariance_log2.tsv** | Expression matrix of DEG∩HighVariance features (log2-transformed). | Feature matrix | Intermediate feature matrix for initial supervised modeling. |
| **ML_FeatureMatrix_DEG_HighVariance_Zscored.tsv** | Z-score standardized version of the previous matrix. | Scaled feature matrix | Ideal for ML algorithms requiring centered and scaled input (e.g., linear models, PCA, regularized regression). |

---

## 3.  Interpretation of Feature Construction

The design of these feature sets reflects a **two-layered rationale**: preserving biological interpretability while ensuring numerical robustness. By starting from DEGs, the pipeline retains **genes with mechanistic relevance**—chiefly those associated with interferon signaling (*IFI44L*, *OAS1*, *MX1*), neutrophil activity (*CEACAM8*, *ELANE*), and immunoregulation (*STAT1*, *IRF7*). 

The hierarchical organization (AllDEGs → HighVariance → Intersection) provides flexibility for subsequent modeling and evaluation of model generalizability versus interpretability. In particular, the **intersection feature matrix** represents a *parsimonious yet biologically faithful embedding* of the transcriptomic signal, suitable for diagnostic or prognostic model development.


In summary, the "DE_Genes/` folder represents a critical interface between **biological discovery** and **predictive modeling**. These feature matrices and gene lists capture the molecular essence of SLE pregnancy at the expression level, providing the foundation for downstream algorithms to learn, differentiate, and ultimately interpret the immunological divergence between healthy and lupus-affected pregnancies.

