# ML PIPELINE INPUT DATA FILE

`ML_FeatureMatrix_AllDEGs_with_Metadata.tsv`

It is the final data file produced by the ML preprocessing script of step [[06] Machine_Learning_Preparation](../[06] Machine_Learning_Preparation). 

It basically combines:

- metadata columns from the `GSE108497_updated_metadata_II.csv`. This was the updated metadata file created from the merging of data from ADEx and NCBI ([step [02] Metadata](../[02] Metadata/01_Metadata_Assembly/)
- DEG expression features from `ML_FeatureMatrix_AllDEGs_log2.tsv`. This was an expression matrix containing log₂-transformed     expression values for all 212 genes selected as Differentially Expressed Genes (DEGs) in the LIMMA analysis [step [04] LIMMA_DE_analysis](../[04] LIMMA_DE_analysis/DE_Genes/).

This format allows ML models to access both sample labels and gene expression features while preserving donor information for proper validation strategies.
