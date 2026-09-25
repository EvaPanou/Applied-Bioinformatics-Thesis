# 02_DGE_Analysis / 02_Results / B_Pooled_SLE_vs_Healthy

## Folder Overview

Method B's result: a single pooled SLE-vs-Healthy comparison across all timepoints (`~0 + Condition`), via limma's `duplicateCorrelation()` with donor blocking — one consensus within-donor correlation estimated across the whole genome and applied to every gene. Produced by Step 5 of `../../01_limma_DEG_pipeline.R`.

## Folder Structure & File Reference

| File | Genes | Description |
|---|---|---|
| `B_FULL_SLE_vs_Healthy_Pooled.tsv` | 9,055 | Every gene tested, unfiltered, ranked by p-value. Columns: `logFC`, `CI.L`, `CI.R` (95% confidence interval — present here and in Method C, unlike Method A's per-timepoint tables, since `topTable(..., confint = TRUE)` was requested for the pooled models specifically), `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `Gene`. |
| `B_UP_FDR0.05_logFC1_SLE_vs_Healthy_Pooled.tsv` | 43 | Up-regulated at FDR<0.05, logFC>1. |
| `B_DN_FDR0.05_logFC-1_SLE_vs_Healthy_Pooled.tsv` | 1 | Down-regulated — just `ORM1`, the same single down-regulated gene found throughout this stage (Method A's `<16 weeks`/`16-23 weeks`, and — as it turns out — Method C below too). |

**43 + 1 = 44**, matching `Method_B_Final.tsv` (44 genes) in `../DE_Genes/` exactly.

---

## Results

**Top of `B_FULL_SLE_vs_Healthy_Pooled.tsv`:**

```text
logFC,CI.L,CI.R,AveExpr,t,P.Value,adj.P.Val,B,Gene
2.65,2.36,2.95,8.92,17.63,5.4e-54,4.8e-50,112.0,IFI44L
2.35,2.08,2.62,11.66,17.34,1.2e-52,5.5e-49,108.9,ISG15
1.83,1.62,2.04,8.40,17.14,1.0e-51,3.0e-48,106.8,OASL
```

`IFI44L` tops Method B's pooled result too, at a comparable effect size (logFC 2.65) to what it showed at most individual timepoints in Method A (2.98–4.04 depending on timepoint) — pooling brings the estimate roughly into the middle of that per-timepoint range, as expected for an inverse-variance-style combination.

**`ORM1`** is the sole down-regulated gene:

```text
logFC,CI.L,CI.R,AveExpr,t,P.Value,adj.P.Val,B,Gene
-1.05,-1.35,-0.74,8.42,-6.75,4.2e-11,1.9e-9,14.5,ORM1
```

**Used downstream by:**
- `B_UP` + `B_DN` genes become `Method_B_Final` (44 genes) in `../DE_Genes/`.
- Compared directly against Method C's equivalent result in `../ForestPlot_MethodB_vs_MethodC.png` (documented in the parent `02_Results/README.md`) — the two methods' pooled effect sizes overlap closely for every gene shown there.
- `Method_B_Final` and Method C's `Method_C_Final` turn out to be identical (44 genes each) — see `../DE_Genes/README.md` for the full cross-method comparison.
