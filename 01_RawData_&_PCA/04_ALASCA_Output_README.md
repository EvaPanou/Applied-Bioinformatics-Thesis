# 04_ALASCA_Output

Repeated-measures multivariate decomposition (RM-ASCA+, via the ALASCA
package) of the corrected non-NP expression matrix, run alongside
`02_DGE_Analysis` and `03_Feature_Validation` as a third, independently-built
line of evidence for the SLE signature -- not built by first selecting genes
that already differ by Condition, but by decomposing the whole transcriptome's
variance into named model terms.

## What ALASCA actually computes -- the mechanism behind every plot here

Three steps happen for every effect (`time_point`, `Condition`,
`time_point:Condition`), always in this order:

1. **Per-gene model fit** (this already happened before any plotting code
   runs -- it's the "Calculating LMM coefficients" line in the console log).
   `value ~ time_point * Condition + (1 | Donor_id)` is fit separately for
   every one of the 9,055 genes, giving each gene its own full set of
   regression coefficients.

2. **Build a small "effect matrix" for one term at a time.** For the
   requested effect, ALASCA uses only the coefficients belonging to that
   term to compute a *predicted* value for every combination of levels that
   term involves, with everything else held at a shared reference point.
   For `time_point:Condition`, that's 5 timepoints x 2 conditions = 10 rows
   -- not 475 rows for individual samples. This is a table of predicted
   *group cells*, not of samples.

3. **PCA on that small matrix, separately from anything else in the
   project.** The 10x9,055 (or 2x9,055, for `Condition` alone) effect matrix
   goes into `prcomp()`. `get_scores(effect, component)` returns each row's
   coordinate on the requested PC -- that's the `score` column, one number
   per group cell.

## Why cells collapse onto identical scores -- the mechanism, not an approximation

Under standard reference coding, `Condition = Healthy` and
`time_point = <16 weeks` are both baseline levels. The interaction
coefficients only exist to describe *SLE's extra deviation at each
non-baseline timepoint* -- there is no such parameter as "Healthy's
interaction at 24-31 weeks" or "SLE's interaction at baseline," because
those cells are defined to have zero interaction contribution by the
model's own parameterization. It isn't that they were estimated and
happened to come out at zero -- there was never a free parameter there to
estimate. So step 2 computes the *identical* predicted value for every
Healthy row and for the SLE/`<16 weeks` row, since all six of those rows
are "no interaction present" cells by construction.

### Worked example: `03_interaction_pca_by_condition.png`

The console output from `str(scores_pc1)` for this effect, verified
directly (not read off the plot):

| time_point | Condition | PC1 score |
|---|---|---|
| <16 weeks | Healthy | 9.34 |
| <16 weeks | SLE | 9.34 |
| 16-23 weeks | Healthy | 9.34 |
| 16-23 weeks | SLE | -7.27 |
| 24-31 weeks | Healthy | 9.34 |

(R's `str()` truncates a 10-element vector to the first 5 by default --
the remaining 5 rows, `24-31wk/SLE` through `PP/SLE`, exist in the full
`interaction_scores` object but weren't printed to console. See "Getting
the exact full table" below for how to pull them without re-running
anything.)

This table explains the plot exactly: **6 of the 10 rows share the value
9.34** -- all 5 Healthy timepoints, plus SLE's own `<16 weeks` row, since
that one is also a "no interaction yet" cell. Plotted with `alpha = 0.8`,
6 overlapping layers at one point render as solid, opaque blue; the single
orange layer underneath (SLE at `<16 weeks`) is fully there in the data,
just visually swamped by the 5 blue layers on top of it. The 4 visibly
separate orange dots are SLE's remaining 4 timepoints (16-23wk, 24-31wk,
32-40wk, PP) -- the only rows where the interaction term actually has
something to estimate.

## On the `.rds` file -- why it can't be opened outside R

`03_alasca_model_no_NP_3effects.rds` is confirmed to be a standard,
uncorrupted gzip-compressed R serialization file (112 MB) -- but attempting
to read it with Python's `pyreadr` fails with *"Invalid file, or file has
unsupported features."* That's expected, not a sign of a bad save: `saveRDS()`
on an ALASCA model saves an **R6 reference-class object** -- a live
environment holding methods (`get_scores()`, `get_loadings()`, `plot()`,
etc.), active bindings, and closures -- not a plain data structure like a
data.frame or list. Tools built to read RDS files outside R (`pyreadr`,
`rio`, similar) only support the portable, "pure data" RDS structures;
reference-class objects with embedded functions and environments can only
be reconstructed by an actual R session with the `ALASCA` package loaded,
since reconstructing the object means re-attaching its methods, not just
its data.

**Getting the exact full table yourself, without re-running the model:**
the object stored in `interaction_scores` at the end of `03_ALASCA_Analysis.R`
already *is* the full 10-row merged table (`plot_group_scores()`'s return
value). If your R session is still open from the run, just type
`interaction_scores` at the console to print all 10 rows with exact
values. If not, `readRDS("03_alasca_model_no_NP_3effects.rds")` reloads the
full model, and `get_scores(alasca_model, effect = 3, component = 1)` /
`component = 2` regenerate the PC1/PC2 tables exactly as they were the
first time -- no re-fitting needed, since the model itself is what's saved.

## Files in this folder

| File | What it is |
|---|---|
| `01_time_point_effect_plot.png` | Pure gestational-time effect (shared across both groups) -- unaffected by the Condition/interaction split, effect 1 |
| `01_condition_effect_plot.png` | Pure Condition main effect, interaction removed -- PC1 explains 100% of this effect's variance, top loadings are the same ISG genes the DEG panel found |
| `01_time_by_condition_effect_plot.png` | The interaction alone -- top loadings are a completely different, non-ISG gene set, meaning the interaction isn't being driven by the SLE signature |
| `02_condition_only_pca_by_condition.png` | 2-point score plot for the Condition-only effect (see mechanism above for why it's exactly 2, not more) |
| `03_interaction_pca_by_condition.png` | 10-row (5 timepoint x 2 condition) score plot for the interaction alone -- see worked example above |
| `03_alasca_model_no_NP_3effects.rds` | The full saved model object -- R-only, see above |
