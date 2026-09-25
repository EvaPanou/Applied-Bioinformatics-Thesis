# 02_DGE_Analysis / 02_Results / A_SLE_vs_Healthy

## Folder Overview

Method A's results: SLE vs. Healthy tested separately at each of the 5 gestational timepoints, via limma with `duplicateCorrelation()` donor blocking. Produced by Step 2–4 of `../../01_limma_DEG_pipeline.R`. For each timepoint, three files: the full ranked gene table, the FDR<0.05/logFC>1 up-regulated subset, and the equivalent down-regulated subset.

## Folder Structure & File Reference

Pure output, 15 files — 3 per timepoint × 5 timepoints, all flat (no further subfolders). Every `FULL` table has the same 7 columns: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `Gene` (limma's standard `topTable()` output, with `Gene` added as an explicit column since row names aren't preserved by `fwrite()`).

| File | Genes | Description |
|---|---|---|
| `A_FULL_SLE_vs_Healthy_at_<16.weeks.tsv` | 9,055 | Every gene tested at `<16 weeks`, unfiltered, ranked by p-value |
| `A_UP_FDR0.05_logFC1_SLE_vs_Healthy_at_<16.weeks.tsv` | 49 | Up-regulated at `<16 weeks` |
| `A_DN_FDR0.05_logFC-1_SLE_vs_Healthy_at_<16.weeks.tsv` | 1 | Down-regulated at `<16 weeks` (`ORM1`) |
| `A_FULL_SLE_vs_Healthy_at_16-23.weeks.tsv` | 9,055 | Every gene tested at `16-23 weeks`, unfiltered |
| `A_UP_FDR0.05_logFC1_SLE_vs_Healthy_at_16-23.weeks.tsv` | 49 | Up-regulated at `16-23 weeks` |
| `A_DN_FDR0.05_logFC-1_SLE_vs_Healthy_at_16-23.weeks.tsv` | 1 | Down-regulated at `16-23 weeks` (`ORM1`) |
| `A_FULL_SLE_vs_Healthy_at_24-31.weeks.tsv` | 9,055 | Every gene tested at `24-31 weeks`, unfiltered |
| `A_UP_FDR0.05_logFC1_SLE_vs_Healthy_at_24-31.weeks.tsv` | 45 | Up-regulated at `24-31 weeks` |
| `A_DN_FDR0.05_logFC-1_SLE_vs_Healthy_at_24-31.weeks.tsv` | 0 | No down-regulated genes at this timepoint |
| `A_FULL_SLE_vs_Healthy_at_32-40.weeks.tsv` | 9,055 | Every gene tested at `32-40 weeks`, unfiltered |
| `A_UP_FDR0.05_logFC1_SLE_vs_Healthy_at_32-40.weeks.tsv` | 32 | Up-regulated at `32-40 weeks` (the smallest UP set of the 5) |
| `A_DN_FDR0.05_logFC-1_SLE_vs_Healthy_at_32-40.weeks.tsv` | 0 | No down-regulated genes at this timepoint |
| `A_FULL_SLE_vs_Healthy_at_PP.tsv` | 9,055 | Every gene tested postpartum, unfiltered |
| `A_UP_FDR0.05_logFC1_SLE_vs_Healthy_at_PP.tsv` | 47 | Up-regulated postpartum |
| `A_DN_FDR0.05_logFC-1_SLE_vs_Healthy_at_PP.tsv` | 0 | No down-regulated genes postpartum |

---

## Results

**`<16 weeks`** — top 2 genes by significance:

```text
logFC,AveExpr,t,P.Value,adj.P.Val,B,Gene
2.98,8.92,13.64,5.9e-36,5.4e-32,70.7,IFI44L
4.04,10.86,13.25,2.8e-34,1.3e-30,66.9,IFI27
```

**`16-23 weeks`**:

```text
logFC,AveExpr,t,P.Value,adj.P.Val,B,Gene
3.06,8.92,13.51,2.1e-35,1.9e-31,69.5,IFI44L
2.35,10.45,13.09,1.3e-33,6.0e-30,65.5,EPSTI1
```

**`32-40 weeks`** (the smallest UP set, 32 genes) — top genes here are different from the earliest timepoints:

```text
logFC,AveExpr,t,P.Value,adj.P.Val,B,Gene
1.19,10.53,9.36,3.2e-19,2.1e-15,32.8,ZBP1
2.02,11.66,9.32,4.6e-19,2.1e-15,32.5,ISG15
```

**`PP`** (postpartum) — signal persists after delivery, though somewhat attenuated (note the smaller B-statistics relative to the pregnancy timepoints):

```text
logFC,AveExpr,t,P.Value,adj.P.Val,B,Gene
3.45,10.86,8.09,5.2e-15,4.7e-11,23.4,IFI27
2.34,8.92,7.64,1.3e-13,5.7e-10,20.3,IFI44L
```

`IFI44L` places in the top 2 genes at 4 of the 5 timepoints (all but `32-40 weeks`) — the single most consistently top-ranked gene across this whole folder, and the same gene topping Method B's pooled result (see `../B_Pooled_SLE_vs_Healthy/README.md`).

**The down-regulated side is essentially empty**: only `<16 weeks` and `16-23 weeks` have any down-regulated gene at all, and both contain exactly one — `ORM1` — with `24-31 weeks`, `32-40 weeks`, and `PP` having zero. See `../UpSet_A_DN_FDR0.05.png` (documented in the parent `02_Results/README.md`) for the visual confirmation.

**Used downstream by:**
- The **UP/DOWN sets across all 5 timepoints** feed `../Heatmap_A_DEG_Union_FDR0.05.png` and both UpSet plots directly.
- The **union** of every UP+DOWN gene across all 5 timepoints becomes `Method_A_Union` (65 genes) in `../DE_Genes/`.
- The **`FULL` tables' logFC and standard error** (derived from `logFC / t`) feed Step 9's random-effects meta-analysis, producing the pooled `Method_A_Final` panel (42 genes) — the per-timepoint results here are never used as final significance calls on their own; they're pooled first.
