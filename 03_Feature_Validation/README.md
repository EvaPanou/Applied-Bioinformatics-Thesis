# 03_Feature_Validation

## Folder Overview

Validates `02_DGE_Analysis`'s gene selection by directly comparing "before" (all 9,055 genes, blind to Condition) against "after" (the DEG panels) on the same sample set, using two complementary tools: PCA (visual) and PERMANOVA (a formal effect-size test). Also produces volcano and MA plots restricted to each "after" panel, to show the spread of effect size/significance within an already-validated gene set rather than to discover anything new. This stage doesn't select or filter genes itself — it's purely a sanity check on stage `02`'s output.

## Folder Structure & File Reference

Per the repo-wide convention, `00_` files are inputs, odd numbers are code, even numbers are outputs:

| File / folder | Type | Description |
|---|---|---|
| `00_final_filtered_expression_non_NP_donors.tsv` | Input | The same 9,055 genes × 475 samples expression matrix used throughout `02_DGE_Analysis` — copied locally rather than referenced across folders. |
| `00_non_NP_metadata.csv` | Input | The same 489-sample metadata used in `02_DGE_Analysis`. |
| `00_Method_All_Union.tsv` | Input | The "lenient" after-panel — 44 genes (A or B or C), copied from `02_DGE_Analysis/02_Results/DE_Genes/`. |
| `00_Method_All_Intersection.tsv` | Input | The "strict" after-panel — 42 genes (A and B and C), copied from the same source. |
| `00_Method_All_Union_Metrics.tsv` | Input | Per-gene statistics spanning all 3 methods — the source for the volcano/MA plots' logFC and significance values, deliberately not tied to any single method's own output table. |
| `01_Feature_Validation.R` | Code | Loads the inputs, then runs 4 steps: volcano/MA plots restricted to each after-panel, PCA on all 3 panels (before + 2 after), and PERMANOVA on all 3 — wrapped in `run_feature_validation()`, called once at the end of the file. |
| `02_Results/` | Output (folder) | 4 subfolders: `MA_Plots/`, `PCA_Plots/`, `PERMANOVA/`, `Volcano_Plots/` — see [`02_Results/README.md`](./02_Results/README.md). |

---

## Methodology

### The "before vs. after" framing, and what it actually proves

Three gene panels are compared throughout this stage:

- **`All_Genes`** (before) — every gene, unfiltered, entirely blind to Condition.
- **`DEG_Union`** (after, lenient) — the 44-gene `Method_All_Union` panel.
- **`Final_Panel`** (after, strict) — the 42-gene `Method_All_Intersection` panel.

The script's own header comment is explicit about calibrating expectations here: with ~44 genes nearly all moving in one coherent direction (the ISG signature), a clean PC1 split in the "after" panels is close to expected *by construction* — genes individually selected because they differ by Condition, in a shared direction, mechanically create Condition-aligned covariance once you restrict to just them. A split appearing isn't surprising on its own. What's genuinely informative instead:

1. Whether the "before" plot reveals a **data-quality problem** independent of Condition (an outlier donor, batch clustering) — real QC value the "after" panels can't give you.
2. The **size of the jump** in PERMANOVA R² from before to after — a quantifiable measure of how much a validated panel concentrates the signal, not just a yes/no on separation.
3. That the "after" R² for Condition is computed **at the donor level**, not a naive sample-level number that repeated timepoints could inflate.

### PCA (`run_pca()`)

Standard `prcomp()` on each panel's expression matrix (centered and scaled, zero-variance genes dropped first), plotted as PC1 vs. PC2, colored by Condition, shaped by timepoint, with a 95% confidence ellipse per Condition (not per Condition×Timepoint combination — the plotting code explicitly groups the ellipse by Condition only, since `shape` is also a discrete aesthetic that would otherwise force 10 separate ellipses).

### PERMANOVA (`run_permanova_condition()` / `run_permanova_time()`)

Two separate tests per panel, each using a permutation scheme that matches what actually varies at what level — the same repeated-measures logic `duplicateCorrelation()`/`dream()` apply in `02_DGE_Analysis`, here applied to a multivariate distance-based test instead of a per-gene linear model:

- **Condition** is a between-donor factor — every sample from one donor shares the same label. Testing at the sample level would let a donor with 5 timepoints cast 5 "votes" for the same label (the same pseudoreplication problem the whole DEG pipeline exists to avoid). Fixed by collapsing each donor to one point (mean expression across that donor's own timepoints) before running `adonis2()` — so every donor contributes exactly one independent observation.
- **Time** is a within-donor factor — collapsing to one point per donor would destroy the thing being tested. Time is tested at the full sample level instead, with `strata = Donor_id` restricting permutations to within each donor, following `vegan::adonis2()`'s own documented approach for a factor nested within a blocking factor.

Both use Euclidean distance and 999 permutations.

### Volcano & MA plots

Both read from `Method_All_Union_Metrics.tsv` — a table spanning all 3 DGE methods, not any single method's own `topTable()` output — restricted to whichever panel's genes. `AveExpr` isn't in that metrics file (it only ever comes from a single method's own output, exactly what this script avoids), so it's computed directly from the raw expression matrix instead, tying it to the data rather than to any one method's fit. Every gene plotted has already passed the DEG panel's own significance filter, so every point is expected to sit past the dashed threshold lines — the plots show the *spread* of effect size/significance within an already-validated panel, not a search for new hits.

**How to run it:** Open in RStudio. **Before running**, edit the hardcoded absolute path at the top (`setwd("C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/03_Feature_Validation")`) — same local-machine-specific issue as the R scripts in `01_RawData_&_PCA` and `02_DGE_Analysis`. Source the whole file; the final line (`run_feature_validation()`) runs all 4 steps automatically.

**Environment:** R, via RStudio. Packages: `tidyverse`, `data.table`, `ggplot2`, `ggrepel` (non-overlapping gene labels on the volcano plot), `vegan` (for `adonis2()`, PERMANOVA). No R or package version numbers pinned in the script.

---

## Key Result

PERMANOVA's Condition R² jumps from **6.1%** (`All_Genes`, before) to **41.1%** (`DEG_Union`) to **41.6%** (`Final_Panel`) — a roughly 7× concentration of signal from restricting to the validated gene panels, at donor level (n=135), p=0.001 for every panel. Time's R² moves the opposite direction (6.6% → 2.1% → 1.6%), confirming the DEG panels concentrate Condition-related variance specifically, not variance in general. See [`02_Results/README.md`](./02_Results/README.md) for the full table and every plot.

## Next Stage

`04_ML_Prerequisites` merges the 44-gene `Method_All_Union` expression matrix with corrected metadata into the ML-ready input table.
