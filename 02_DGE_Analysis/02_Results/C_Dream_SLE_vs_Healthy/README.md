# 02_DGE_Analysis / 02_Results / C_Dream_SLE_vs_Healthy

## Folder Overview

Method C's result: the same pooled question as Method B (SLE vs. Healthy, all timepoints combined), but fit via `variancePartition::dream()` — a true linear mixed model per gene, estimating the within-donor correlation separately for each gene rather than one shared genome-wide value. Produced by Step 6 of `../../01_limma_DEG_pipeline.R`.

## Folder Structure & File Reference

| File | Genes | Description |
|---|---|---|
| `C_FULL_SLE_vs_Healthy_Dream.tsv` | 9,055 | Every gene tested, unfiltered, ranked by p-value. Same columns as Method B's `FULL` table, plus one extra: `z.std` (the mixed-model's standardized z-statistic, specific to `dream()`'s output). |
| `C_UP_FDR0.05_logFC1_SLE_vs_Healthy_Dream.tsv` | 43 | Up-regulated at FDR<0.05, logFC>1. |
| `C_DN_FDR0.05_logFC-1_SLE_vs_Healthy_Dream.tsv` | 1 | Down-regulated — `ORM1` again, the same gene flagged throughout this stage. |

**43 + 1 = 44**, matching `Method_C_Final.tsv` in `../DE_Genes/` — and, checked directly (not just by count), **the exact same 44 genes as `Method_B_Final.tsv`.** Methods B and C fully converge on this dataset despite their structurally different correlation-estimation approach.

No `MethodC_dream_convergence_failures.tsv` file exists in this folder — the parent script only writes one if `dream()` fails to converge for any gene, and checks explicitly (`if (!is.null(dream_gene_errors) && length(dream_gene_errors) > 0)`). **Its absence here confirms `dream()` converged successfully for all 9,055 genes.**

---

## Results

**Top of `C_FULL_SLE_vs_Healthy_Dream.tsv`** — note the ranking differs from Methods A and B, which were both topped by `IFI44L`:

```text
logFC,CI.L,CI.R,AveExpr,t,P.Value,adj.P.Val,B,z.std,Gene
1.22,1.01,1.44,9.47,11.29,3.0e-21,1.6e-17,46.9,9.46,SAMD9L
1.98,1.64,2.33,10.45,11.25,3.5e-21,1.6e-17,46.6,9.45,EPSTI1
```

`SAMD9L` and `EPSTI1` top Method C specifically — a genuine ranking difference from Methods A/B (both topped by `IFI44L`), reflecting that per-gene correlation estimation can shift which gene comes out "most significant" even when the overall panel converges. This kind of reordering-without-panel-divergence is exactly what a well-behaved analysis should show: the *significant set* is robust to the choice of correlation model, even if the exact top rank isn't identical gene-for-gene.

**`ORM1`, down-regulated:**

```text
logFC,CI.L,CI.R,AveExpr,t,P.Value,adj.P.Val,B,z.std,Gene
-1.06,-1.44,-0.69,8.42,-5.60,1.1e-07,3.9e-06,7.83,-5.30,ORM1
```

**Used downstream by:**
- `C_UP` + `C_DN` genes become `Method_C_Final` (44 genes) in `../DE_Genes/`.
- Compared directly against Method B in `../ForestPlot_MethodB_vs_MethodC.png`.
- `Method_B_C_Intersection` in `../DE_Genes/` is — given the exact identity confirmed above — the same 44 genes as either method's `Final` list individually; the intersection doesn't remove anything here.
