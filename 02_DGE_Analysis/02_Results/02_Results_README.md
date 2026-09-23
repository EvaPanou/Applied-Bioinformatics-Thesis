# 02_Results

All output from `01_limma_DEG_pipeline_v3.R`. Five subfolders (one per method,
plus two combining folders) and six top-level files.

## Subfolders

### `A_SLE_vs_Healthy/` — Method A, per-timepoint results
For each of the 5 gestational timepoints (`<16 weeks`, `16-23 weeks`,
`24-31 weeks`, `32-40 weeks`, `PP`):
- `A_FULL_SLE_vs_Healthy_at_<timepoint>.tsv` — every gene tested, unfiltered, ranked by p-value
- `A_UP_FDR0.05_logFC1_..._at_<timepoint>.tsv` — up-regulated genes passing the threshold
- `A_DN_FDR0.05_logFC-1_..._at_<timepoint>.tsv` — down-regulated genes passing the threshold

15 files total (3 per timepoint x 5 timepoints).

### `B_Pooled_SLE_vs_Healthy/` — Method B, pooled result
- `B_FULL_SLE_vs_Healthy_Pooled.tsv` — every gene, unfiltered
- `B_UP_FDR0.05_logFC1_SLE_vs_Healthy_Pooled.tsv` / `B_DN_...` — filtered up/down lists

### `C_Dream_SLE_vs_Healthy/` — Method C, pooled result
- `C_FULL_SLE_vs_Healthy_Dream.tsv` — every gene, unfiltered
- `C_UP_FDR0.05_logFC1_SLE_vs_Healthy_Dream.tsv` / `C_DN_...` — filtered up/down lists
- `MethodC_dream_convergence_failures.tsv` — only appears if any gene failed to fit; absent means every gene converged

### `Gene_Clusters/`
Method A's DEG union, split into up to 3 groups by expression pattern
(k-means clustering on the heatmap data): `Heatmap_Cluster_Up_regulated_in_SLE_genes.tsv`,
`..._Mixed_time_dependent_genes_genes.tsv`, `..._Down_regulated_in_SLE_genes.tsv`.
Categorized raw gene lists for manual review, not an automated enrichment result.

### `DE_Genes/` — the combined, cross-method results (start here for downstream work)

**Naming convention:** `Method_<A/B/C/All/B_C>_<Union/Intersection/Final>[_Metrics]`

| File | Genes | What it contains |
|---|---|---|
| `Method_A_Union.tsv` | 65 | Program A, significant at >=1 of 5 timepoints (before pooling) |
| `Method_A_Union_Metrics.tsv` | 65 | Same 65 genes, with meta-analysis statistics attached |
| `Method_A_Final.tsv` | 42 | Program A, properly pooled via meta-analysis |
| `Method_B_Final.tsv` | 44 | Method B's final list |
| `Method_C_Final.tsv` | 44 | Method C's final list |
| `Method_All_Intersection.tsv` | 42 | A and B and C (**the primary validated panel**) |
| `Method_B_C_Intersection.tsv` | 44 | B and C |
| `Method_All_Union.tsv` | 44 | A or B or C |
| `Method_All_Union_Metrics.tsv` | 44 | Same as above, with per-method membership flags (`In_A`/`In_B`/`In_C`) and Method A's heterogeneity statistics attached |
| `All_Genes_Metrics.tsv` | ~9,055 | Complete, unfiltered meta-analysis reference (every gene tested, significant or not) — the one exception to the naming convention above, since it isn't a combination of anything |

Each of the 7 gene-list files (all except the two `_Metrics` files and
`All_Genes_Metrics`) has a matching `ML_FeatureMatrix_<name>_log2.tsv` —
samples x genes, log2 expression, ready to use as ML input.

One additional file, `ML_FeatureMatrix_Method_All_Intersection_with_Metadata.tsv`,
is formatted specifically for the downstream Python ML pipeline (adds
`SampleID`/`DonorID`/`Condition`/`Condition_label`/`Timepoint`/`Time_label`
columns alongside the 42-gene expression data).

**For downstream work:** use `Method_All_Intersection` (42 genes) as the
primary panel — every gene here is independently confirmed by all three
statistical methods, the strongest validated result this analysis produces.

## Top-level files

| File | What it shows |
|---|---|
| `Heatmap_A_DEG_Union_FDR0.05.png` | Method A's 65-gene union, z-scored, split by Condition and timepoint |
| `Heatmap_FinalPanel_ConditionOnly.png` | The 42-gene final panel, Condition-only (no timepoint split) |
| `UpSet_A_UP_FDR0.05.png` / `UpSet_A_DN_FDR0.05.png` | Overlap of significant genes across Method A's 5 timepoints |
| `ForestPlot_MethodB_vs_MethodC.png` | Method B vs. Method C pooled effect sizes (logFC with 95% CI), top genes |
| `Global_Time_by_Condition_Ftest.tsv` | Per-gene test for a Condition x Time interaction (all ~9,055 genes); 0 genes reach significance at FDR<0.05 |
