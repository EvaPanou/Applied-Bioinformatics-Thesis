# 01_RawData_&_PCA / 04_ALASCA_Output

## Folder Overview

Every output of `../03_ALASCA_Analysis.R`: a repeated-measures multivariate decomposition (RM-ASCA+, via the ALASCA package) of the corrected non-NP expression matrix, run alongside `02_DGE_Analysis` and `03_Feature_Validation` as a third, independently-built line of evidence for the SLE signature — not built by first selecting genes that already differ by Condition, but by decomposing the whole transcriptome's variance into named model terms.

This folder previously had two versions of its own README uploaded (`04_ALASCA_Output_README.md` and `04_ALASCA_Output_README (1).md`); the `(1)` version is the complete one (333 lines, including a literature section and a methodological discussion missing from the shorter file), and this document is built from that fuller version. Worth deleting the shorter duplicate from the actual repo.

## Folder Structure

Pure output — numbered by generation order rather than the 00/odd/even convention:

```text
04_ALASCA_Output/
├── 01_time_point_effect_plot.png
├── 01_condition_effect_plot.png
├── 01_time_by_condition_effect_plot.png
├── 02_condition_only_pca_by_condition.png
├── 03_interaction_pca_by_condition.png
└── 03_alasca_model_no_NP_3effects.rds
```

## File Reference

| File | What it is |
|---|---|
| `01_time_point_effect_plot.png` | Pure gestational-time effect (shared across both groups) — unaffected by the Condition/interaction split, effect 1 |
| `01_condition_effect_plot.png` | Pure Condition main effect, interaction removed — PC1 explains 100% of this effect's variance, top loadings are the same ISG genes the DEG panel found |
| `01_time_by_condition_effect_plot.png` | The interaction alone — top loadings are a completely different, non-ISG gene set, meaning the interaction isn't being driven by the SLE signature |
| `02_condition_only_pca_by_condition.png` | 2-point score plot for the Condition-only effect (see mechanism below for why it's exactly 2, not more) |
| `03_interaction_pca_by_condition.png` | 10-row (5 timepoint × 2 condition) score plot for the interaction alone — see worked example below |
| `03_alasca_model_no_NP_3effects.rds` | The full saved model object — R-only, see below |

---

## What ALASCA actually computes — the mechanism behind every plot here

Three steps happen for every effect (`time_point`, `Condition`, `time_point:Condition`), always in this order:

1. **Per-gene model fit** — `value ~ time_point * Condition + (1 | Donor_id)` is fit separately for every one of the 9,055 genes, giving each gene its own full set of regression coefficients.
2. **Build a small "effect matrix" for one term at a time.** For the requested effect, ALASCA uses only the coefficients belonging to that term to compute a *predicted* value for every combination of levels that term involves, with everything else held at a shared reference point. For `time_point:Condition`, that's 5 timepoints × 2 conditions = 10 rows — not 475 rows for individual samples. This is a table of predicted *group cells*, not of samples.
3. **PCA on that small matrix, separately from anything else in the project.** The 10×9,055 (or 2×9,055, for `Condition` alone) effect matrix goes into `prcomp()`. `get_scores(effect, component)` returns each row's coordinate on the requested PC.

## Why cells collapse onto identical scores — the mechanism, not an approximation

Under standard reference coding, `Condition = Healthy` and `time_point = <16 weeks` are both baseline levels. The interaction coefficients only exist to describe *SLE's extra deviation at each non-baseline timepoint* — there is no such parameter as "Healthy's interaction at 24-31 weeks" or "SLE's interaction at baseline," because those cells are defined to have zero interaction contribution by the model's own parameterization. It isn't that they were estimated and happened to come out at zero — there was never a free parameter there to estimate. So step 2 computes the *identical* predicted value for every Healthy row and for the SLE/`<16 weeks` row, since all six of those rows are "no interaction present" cells by construction.

### Worked example: `03_interaction_pca_by_condition.png`

![Group scores — Time x Condition interaction](03_interaction_pca_by_condition.png)

The console output from `str(scores_pc1)` for this effect (as recorded when the script was run — R's `str()` truncates a 10-element vector to the first 5 by default):

| time_point | Condition | PC1 score |
|---|---|---|
| <16 weeks | Healthy | 9.34 |
| <16 weeks | SLE | 9.34 |
| 16-23 weeks | Healthy | 9.34 |
| 16-23 weeks | SLE | -7.27 |
| 24-31 weeks | Healthy | 9.34 |

This explains the plot exactly: **6 of the 10 rows share the value 9.34** — all 5 Healthy timepoints, plus SLE's own `<16 weeks` row, since that one is also a "no interaction yet" cell. Plotted with `alpha = 0.8`, 6 overlapping layers at one point render as solid, opaque blue; the single orange layer underneath (SLE at `<16 weeks`) is fully there in the data, just visually swamped by the 5 blue layers on top of it — visible in the plot above as the single blue dot around (9, 1) with four separate orange dots scattered below and to the left of it (16-23wk, 24-31wk, 32-40wk, PP) — the only rows where the interaction term actually has something to estimate.

### The other four plots

**`01_time_point_effect_plot.png`** — pure gestational-time effect, shared across both groups by construction (hence only one line/color in the legend):

![Pure time effect](01_time_point_effect_plot.png)

PC1 explains 74.24% of this effect's variance. Top loadings are myeloid/neutrophil-associated genes (`MCEMP1`, `SLC22A4`, `S100A8`, among others) — a plausible pregnancy-associated innate-immune signature, structurally distinct from the interferon-stimulated genes driving the Condition effect below.

**`01_condition_effect_plot.png`** — the pure SLE-vs-Healthy main effect, interaction removed:

![Pure Condition effect](01_condition_effect_plot.png)

PC1 explains 100% of this effect's variance (expected — a 2-level factor with no interaction has exactly one degree of freedom to distribute variance across). Top loadings are the same interferon-stimulated genes (`IFI44L`, `IFI27`, `EPSTI1`, `RSAD2`, `IFI44`, and others) found throughout `02_DGE_Analysis`'s DEG panel — this is the ALASCA-side confirmation of the same signal, arrived at independently, without first selecting genes by differential expression.

**`01_time_by_condition_effect_plot.png`** — the interaction alone, across all 5 timepoints:

![Interaction effect across timepoints](01_time_by_condition_effect_plot.png)

PC1 explains 44.3% of this effect's variance. Top loadings (`PRMT9`, `SEC61G`, `NDC1`, `MPC1`, and others) are a completely different, non-interferon gene set from the Condition effect's loadings above — visual confirmation that whatever the interaction term is capturing, it isn't acting on the SLE signature itself.

**`02_condition_only_pca_by_condition.png`** — the same Condition effect as a 2-point score plot:

![Condition-only group scores](02_condition_only_pca_by_condition.png)

Exactly 2 points (Healthy, SLE) by construction — a factor with 2 levels and no interaction produces one score per level. The y-axis sits at ~1e-13 — numerically zero, confirming PC1 alone captures essentially all of this effect's variance, consistent with the 100% figure above.

## On the `.rds` file — why it can't be opened outside R

`03_alasca_model_no_NP_3effects.rds` is confirmed to be a standard, uncorrupted gzip-compressed R serialization file (112 MB) — but attempting to read it with Python's `pyreadr` fails with *"Invalid file, or file has unsupported features."* That's expected, not a sign of a bad save: `saveRDS()` on an ALASCA model saves an **R6 reference-class object** — a live environment holding methods (`get_scores()`, `get_loadings()`, `plot()`, etc.), active bindings, and closures — not a plain data structure like a data.frame or list. Tools built to read RDS files outside R (`pyreadr`, `rio`, similar) only support the portable, "pure data" RDS structures; reference-class objects with embedded functions and environments can only be reconstructed by an actual R session with the `ALASCA` package loaded, since reconstructing the object means re-attaching its methods, not just its data.

**Getting the exact full table yourself, without re-running the model:** the object stored in `interaction_scores` at the end of `../03_ALASCA_Analysis.R` already *is* the full 10-row merged table (`plot_group_scores()`'s return value). If your R session is still open from the run, just type `interaction_scores` at the console to print all 10 rows with exact values. If not, `readRDS("03_alasca_model_no_NP_3effects.rds")` reloads the full model, and `get_scores(alasca_model, effect = 3, component = 1)` / `component = 2` regenerate the PC1/PC2 tables exactly as they were the first time — no re-fitting needed, since the model itself is what's saved.

## Literature backing this methodology

- **RM-ASCA+ itself**: Madssen, T.S., Giskeødegård, G.F., Smilde, A.K., & Westerhuis, J.A. (2021). Repeated measures ASCA+ for analysis of longitudinal intervention studies with multivariate outcome data. *PLOS Computational Biology*, 17(11), e1009585. This is the original method paper — the three-step decomposition described above (per-gene fit → effect matrix → PCA on the effect matrix) is exactly the procedure this paper introduces, not something specific to the R implementation.
- **The ALASCA R package** (what actually ran here): Jarmund, A.H., Madssen, T.S., & Giskeødegård, G.F. (2022). ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods. *Frontiers in Molecular Biosciences*, 9, 962431.
- **Foundational ASCA** (the non-repeated-measures ancestor both papers above build on): Smilde, A.K., Jansen, J.J., Hoefsloot, H.C., Lamers, R.J., Van Der Greef, J., & Timmerman, M.E. (2005). ANOVA-simultaneous component analysis (ASCA): a new tool for analyzing designed metabolomics data. *Bioinformatics*, 21(13), 3043-3048.
- **Reference/dummy coding and why interaction cells collapse the way they do**: standard linear-model theory, e.g. Faraway, J.J. (2014). *Linear Models with R* (2nd ed.). CRC Press — not specific to ALASCA, this is the same coding scheme underlying `time_point * Condition` in any R model formula, including the limma models in `02_DGE_Analysis`.

## Why this justifies pooling over time in the DEG comparisons

Two structurally independent methods now agree that Condition's effect on gene expression does not meaningfully change shape across gestation:

1. `02_DGE_Analysis`'s `Global_Time_by_Condition_Ftest.tsv` — a per-gene univariate test — found 0 of 9,055 genes with a significant Condition × Time interaction.
2. This ALASCA run — a multivariate decomposition of the whole transcriptome at once, built without selecting any genes in advance — shows the interaction effect's top-loaded genes are a completely different, non-ISG set from the Condition effect's, meaning the interaction isn't acting on the SLE signature at all.

The general statistical principle this supports — that a non-significant interaction term should be dropped in favor of the simpler, pooled model, and that *failing* to drop it can itself distort the main-effect estimate — is argued directly in Engqvist, L. (2005). The mistreatment of covariate interaction terms in linear model analyses of behavioural and evolutionary ecology studies. *Animal Behaviour*, 70(4), 967-971.

**One honesty caveat worth keeping attached to this, since Engqvist's own paper is partly a warning against exactly this kind of reasoning done carelessly**: the concern that paper raises is about *post-hoc* interaction removal — dropping a term only after seeing it's inconvenient, then presenting the simplified model as if it were the plan all along, which inflates false confidence. That critique doesn't apply here in the way it usually would: Method B and Method C's pooled, no-interaction models (see `02_DGE_Analysis`) were never a simplification chosen *because* the interaction turned out non-significant — they were two of three co-equal, independently specified methods from the start, run in parallel with Method A's per-timepoint approach, not after it. This ALASCA result and the earlier F-test are better read as *confirmatory* evidence that the pooled design was reasonable all along, not as the justification that produced it after the fact.

**Practical upshot for further Condition-only comparisons**: pooling over time (Method B/C's `~0 + Condition` approach, no interaction term) is the statistically supported default going forward, not something that needs timepoint stratification re-added "just in case" — two independent lines of evidence now back the assumption that would be needed to justify it.
