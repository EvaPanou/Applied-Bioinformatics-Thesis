# 02_DGE_Analysis / 02_Results / DE_Genes

## Folder Overview

The combined, cross-method results — **start here for any downstream work.** Every gene-selection combination across Methods A, B, and C (union, intersection, and named subsets), each paired with a ready-to-use ML feature matrix. Produced by Step 10 of `../../01_limma_DEG_pipeline.R`.

## Folder Structure & File Reference

This folder doesn't follow the 00/odd/even convention — it's pure output, all from one script step, organized instead by a naming convention: `Method_<A/B/C/All/B_C>_<Union/Intersection/Final>[_Metrics]`.

| File | Genes/Rows | Description |
|---|---|---|
| `Method_A_Union.tsv` | 65 | Program A, significant at ≥1 of 5 timepoints (before pooling) |
| `Method_A_Union_Metrics.tsv` | 65 | Same 65 genes, with meta-analysis statistics attached |
| `ML_FeatureMatrix_Method_A_Union_log2.tsv` | 475 samples × 65 genes | Expression matrix for the `Method_A_Union` panel |
| `Method_A_Final.tsv` | 42 | Program A, properly pooled via meta-analysis |
| `ML_FeatureMatrix_Method_A_Final_log2.tsv` | 475 × 42 | Expression matrix for `Method_A_Final` |
| `Method_B_Final.tsv` | 44 | Method B's final list |
| `ML_FeatureMatrix_Method_B_Final_log2.tsv` | 475 × 44 | Expression matrix for `Method_B_Final` |
| `Method_C_Final.tsv` | 44 | Method C's final list (identical set to B) |
| `ML_FeatureMatrix_Method_C_Final_log2.tsv` | 475 × 44 | Expression matrix for `Method_C_Final` |
| `Method_B_C_Intersection.tsv` | 44 | B ∩ C (same 44, since B and C fully agree) |
| `ML_FeatureMatrix_Method_B_C_Intersection_log2.tsv` | 475 × 44 | Expression matrix for `Method_B_C_Intersection` |
| `Method_All_Intersection.tsv` | 42 | A ∩ B ∩ C (**the primary validated panel**) |
| `ML_FeatureMatrix_Method_All_Intersection_log2.tsv` | 475 × 42 | Expression matrix for `Method_All_Intersection` |
| `Method_All_Union.tsv` | 44 | A ∪ B ∪ C |
| `Method_All_Union_Metrics.tsv` | 44 | Same 44, with per-method membership flags (`In_A`/`In_B`/`In_C`) and Method A's heterogeneity statistics attached |
| `ML_FeatureMatrix_Method_All_Union_log2.tsv` | 475 × 44 | Expression matrix for `Method_All_Union` — **the file `04_ML_Prerequisites` actually builds on** |
| `All_Genes_Metrics.tsv` | 9,055 | Complete, unfiltered meta-analysis reference (every gene tested, significant or not) — the one file that doesn't follow the naming convention above, since it isn't a combination of anything |
| `ML_FeatureMatrix_Method_All_Intersection_with_Metadata.tsv` | 475 | The 42-gene panel plus `SampleID`/`DonorID`/`Condition`/`Condition_label`/`Timepoint`/`Time_label` — built specifically for `sle_pipeline.py`'s loader, which expects those 6 metadata columns alongside expression data in one file |

---

## Results

### The core cross-method comparison (Step 10's own console message)

```text
Final A (meta-analysis): 42 | Final B (dupCor pooled): 44 | Final C (dream pooled): 44 | All three agree: 42 | B and C agree: 44
```

**Method A's pooled result (42) is a strict subset of Method B/C's shared 44** — every gene Method A's meta-analysis confirms is also confirmed by both B and C, but B and C additionally agree on 2 more genes that A's per-timepoint-then-pooled approach doesn't independently reach. This is why `Method_All_Intersection` (A ∩ B ∩ C = 42) is narrower than `Method_All_Union` (A ∪ B ∪ C = 44) by exactly those same 2 genes.

**`Method_All_Union_Metrics.tsv`** — header and first 2 rows:

```text
Gene,In_A,In_B,In_C,N_methods,pooled_logFC,adj.P.Val,I2,QEp,Low_heterogeneity
BATF2,TRUE,TRUE,TRUE,3,1.344,3.8e-42,44.97,0.121,TRUE
CMPK2,TRUE,TRUE,TRUE,3,1.487,1.6e-40,47.73,0.108,TRUE
```

This is the single most useful file for understanding *why* a gene made it into the final panel — `N_methods` shows how many of the 3 methods agree on it, and `Low_heterogeneity` (QEp > 0.05) flags whether Method A's per-timepoint estimates were consistent with each other for that gene, i.e. whether it's a genuinely time-invariant signal rather than one propped up by averaging across inconsistent timepoints.

**`Method_A_Union_Metrics.tsv`** — the 65-gene union with the same meta-analysis stats:

```text
Gene,pooled_logFC,pooled_pval,I2,tau2,QEp,adj.P.Val,Low_heterogeneity
HES4,1.680,1.0e-69,6.28,0.00287,0.375,7.7e-67,TRUE
HERC6,1.363,1.4e-50,48.30,0.01994,0.100,4.1e-48,TRUE
```

### A worth-knowing observation: two separate ML-export paths from this stage

The script builds ML-ready output two different ways, for two different consumers, and it's easy to conflate them:

1. **The 7 plain `ML_FeatureMatrix_<name>_log2.tsv` files** — one per candidate gene list, expression-only, no metadata. These are generic and panel-agnostic.
2. **`ML_FeatureMatrix_Method_All_Intersection_with_Metadata.tsv`** — a single, separately-built file specifically formatted for `sle_pipeline.py`'s loader (which requires 6 exact metadata columns alongside expression). Its panel is controlled by one variable in the script, `ML_PANEL_FOR_PYTHON <- "Method_All_Intersection"` — the **42-gene** intersection panel, by default.

**What actually got used downstream doesn't match either of these directly**: `04_ML_Prerequisites`'s own README states its input is `ML_FeatureMatrix_Method_All_Union_log2.tsv` — the plain, metadata-free **44-gene** union matrix (option 1 above), not the pre-built "with_Metadata" convenience file (which defaults to the 42-gene intersection). `04_ML_Prerequisites` then does its own metadata merge from scratch rather than using the ready-made file this script already built for that purpose. Not an error — both files are valid, real outputs of this script — just worth knowing if you're tracing exactly which panel size (42 vs. 44) shows up at which later stage. The 44-gene `Method_All_Union` panel is the one that actually propagates forward through `04_ML_Prerequisites` and `05_ML_Pipeline`.

**Used downstream by:** `04_ML_Prerequisites` merges `ML_FeatureMatrix_Method_All_Union_log2.tsv` (44 genes, 475 samples) with corrected metadata to build the ML pipeline's input table — this is the file that determines every gene `05_ML_Pipeline`'s benchmarking and final-panel derivation can possibly select from.
