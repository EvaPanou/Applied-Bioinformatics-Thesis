## ML PIPELINE INPUT DATA FILE

**ML_FeatureMatrix_AllDEGs_with_Metadata.tsv**

This is the final data file produced by the ML preprocessing script of step  
[\[06\] Machine_Learning_Preparation](<../../[06] Machine_Learning_Preparation/>).

It combines:

- Metadata columns from `GSE108497_updated_metadata_II.csv`.  
  This updated metadata file was created through the merging of data from ADEx and NCBI  
  ([step \[02\] Metadata](<../../[02] Metadata/01_Metadata_Assembly/>)).

- DEG expression features from `ML_FeatureMatrix_AllDEGs_log2.tsv`.  
  This expression matrix contains log₂-transformed expression values for all 212 genes identified as Differentially Expressed Genes (DEGs) during the LIMMA analysis  
  ([step \[04\] LIMMA_DE_analysis](<../../[04] LIMMA_DE_analysis/DE_Genes/>)).

This format allows machine learning models to access both sample labels and gene expression features while preserving donor information required for proper validation strategies.
