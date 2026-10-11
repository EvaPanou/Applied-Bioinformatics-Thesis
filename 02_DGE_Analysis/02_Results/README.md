# 02_DGE_Analysis / 02_Results

## Folder Overview

All output from `../01_limma_DEG_pipeline.R`. Five subfolders (one per method, plus two combining/summary folders) and six top-level files. This folder itself holds only the six top-level cross-method summary outputs — the bulk of the actual per-method results live in the subfolders, each documented in its own README.

## Folder Structure & File Reference

| File / folder | Description |
|---|---|
| [`A_SLE_vs_Healthy/`](./A_SLE_vs_Healthy/README.md) | Method A — per-timepoint results (15 files) |
| [`B_Pooled_SLE_vs_Healthy/`](./B_Pooled_SLE_vs_Healthy/README.md) | Method B — pooled duplicateCorrelation result (3 files) |
| [`C_Dream_SLE_vs_Healthy/`](./C_Dream_SLE_vs_Healthy/README.md) | Method C — pooled dream() result (3 files) |
| [`Gene_Clusters/`](./Gene_Clusters/README.md) | Method A's DEG union, split by expression pattern (1 file — see its own README for why) |
| [`DE_Genes/`](./DE_Genes/README.md) | Combined, cross-method results — start here for downstream work (18 files) |
| `Heatmap_A_DEG_Union_FDR0.05.png` | Method A's 65-gene union, z-scored, split by Condition and timepoint |
| `Heatmap_FinalPanel_ConditionOnly.png` | The 42-gene final panel, Condition-only (no timepoint split) |
| `UpSet_A_UP_FDR0.05.png` | Overlap of significant up-regulated genes across Method A's 5 timepoints |
| `UpSet_A_DN_FDR0.05.png` | Overlap of significant down-regulated genes across Method A's 5 timepoints |
| `ForestPlot_MethodB_vs_MethodC.png` | Method B vs. Method C pooled effect sizes (logFC with 95% CI), top genes |
| `Global_Time_by_Condition_Ftest.tsv` | Per-gene test for a Condition × Time interaction (all ~9,055 genes); 0 genes reach significance at FDR<0.05 |

---

## Results — top-level files

### `Global_Time_by_Condition_Ftest.tsv`

The formal joint F-test (Step 12 of the parent script) for whether the SLE-vs-Healthy effect changes shape across gestation. 9,055 rows, one per gene, ranked by F-statistic. Top row:

```text
ConditionSLE.Time16.23.weeks,ConditionSLE.Time24.31.weeks,ConditionSLE.Time32.40.weeks,ConditionSLE.TimePP,AveExpr,F,P.Value,adj.P.Val,Gene
-0.007,0.049,0.055,0.470,13.549,7.229,1.18e-05,0.107,S100A8
```

**Even the single most significant gene by this test, `S100A8`, has adj.P.Val = 0.107 — above the 0.05 threshold.** No gene in the entire 9,055-gene set reaches significance, confirming quantitatively what the per-gene meta-analysis heterogeneity stats in `DE_Genes/Method_A_Union_Metrics.tsv` suggest qualitatively: the Condition effect is time-invariant across gestation. **This is the direct evidence behind pooling over time in Methods B and C**, and behind treating `01_RawData_&_PCA/04_ALASCA_Output`'s multivariate interaction-effect finding as confirmatory rather than a separate discovery.

### `Heatmap_A_DEG_Union_FDR0.05.png`

Method A's 65-gene union (all genes significant at ≥1 of the 5 timepoints), z-scored per gene, columns split by Condition then by timepoint within each:

![Heatmap of Method A DEG union](Heatmap_A_DEG_Union_FDR0.05.png)

A clean, consistent red (up-regulated)/blue (down-regulated) split between the SLE and Healthy column blocks across every timepoint — visually, the signal doesn't appear to shift in strength or direction across gestation, consistent with the F-test result above. The row dendrogram shows finer clustering by magnitude within the broadly up-regulated set, but — as established in `Gene_Clusters/README.md` — every one of these 65 genes' k-means cluster ultimately gets the same "Up-regulated in SLE" label; none of the visible row substructure crosses into genuinely down-regulated territory except `ORM1`.

### `Heatmap_FinalPanel_ConditionOnly.png`

The 42-gene final panel (`Method_All_Intersection`), Condition-only, no timepoint split:

![Heatmap of the final 42-gene panel](Heatmap_FinalPanel_ConditionOnly.png)

A sharp, near-binary split — SLE columns solid red, Healthy columns solid blue, for essentially every one of the 42 genes. This is the cleanest visual summary of the panel's separating power prior to any machine-learning step. The 42 gene labels visible in the plot (`IFITM3`, `DDX58`, `IRF7`, `OASL`, `ZBP1`, `HELZ2`, `SAMD9L`, `PARP14`, `EIF2AK2`, `BATF2`, `HES4`, `IFI27`, `IFI6`, `ISG15`, `IFI44L`, `XAF1`, `PARP12`, `LY6E`, `HERC6`, `IFI44`, `RSAD2`, `HERC5`, `IFIT1`, `EPSTI1`, `OAS1`, `OAS2`, `OAS3`, `MX1`, `TIMM10`, `USP18`, `SPATS2L`, `RTP4`, `OTOF`, `IFIT5`, `TRIM6`, `IFIT2`, `IFIH1`, `DDX60`, `TRIM22`, `CMPK2`, `IFIT3`, plus `DHX58`) match `DE_Genes/Method_All_Intersection.tsv` exactly.

### `UpSet_A_UP_FDR0.05.png`

Overlap of up-regulated genes across Method A's 5 timepoint contrasts:

![UpSet plot of up-regulated genes per timepoint](UpSet_A_UP_FDR0.05.png)

**31 genes are up-regulated at all 5 timepoints simultaneously** — the largest single bar by a wide margin — with the remaining bars (8, 7, 4, 3, and a run of smaller intersections down to 1) representing genes significant at some but not all timepoints. This is the visual, per-gene version of what the meta-analysis and F-test establish statistically: a large, stable core signal, plus a smaller shoulder of timepoint-specific or borderline genes.

### `UpSet_A_DN_FDR0.05.png`

The equivalent view for down-regulated genes:

![UpSet plot of down-regulated genes per timepoint](UpSet_A_DN_FDR0.05.png)

Only two timepoints (`<16 weeks`, `16-23 weeks`) have any down-regulated gene at all, and both show exactly one gene — the same gene in both cases, `ORM1` (confirmed directly from `A_SLE_vs_Healthy/A_DN_FDR0.05_logFC-1_SLE_vs_Healthy_at_<16.weeks.tsv` and the `16-23 weeks` equivalent). The other three timepoints (`24-31 weeks`, `32-40 weeks`, `PP`) have zero down-regulated genes and don't appear in this plot at all. This is the direct, per-timepoint evidence behind the dataset's overwhelmingly up-regulated signature, referenced throughout this project.

### `ForestPlot_MethodB_vs_MethodC.png`

Method B (duplicateCorrelation) vs. Method C (dream) pooled effect sizes, top genes by either method, with 95% confidence intervals:

![Forest plot comparing Method B and Method C](ForestPlot_MethodB_vs_MethodC.png)

Green (Method B) and orange (Method C) point estimates and intervals overlap closely for essentially every gene shown, from `IFI27` (the largest effect, logFC ≈3.3–3.4) down to `HELZ2` (the smallest shown, logFC ≈1.0–1.1) — visual confirmation that the two structurally different correlation-estimation approaches (one shared genome-wide value vs. one value per gene) converge on the same effect-size estimates in practice for this dataset, consistent with `Method_B_Final` and `Method_C_Final` both landing on the same 44 genes (see `DE_Genes/README.md`).
