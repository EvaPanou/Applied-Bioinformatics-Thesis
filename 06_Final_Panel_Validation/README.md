# 06_Final_Panel_Validation

## Overview

Confirms that the final 3-gene panel — **HERC5, BATF2, SPATS2L** — tracks SLE status, and not pregnancy complication severity, pregnancy status alone, or technical batch. Two focused statistical tests, covering only the 3 genes that made the final panel, with batch folded in as a covariate in both rather than tested separately. This notebook deliberately does **not** re-test SLE vs. Healthy overall — that's already established on the full transcriptome by the primary DE/ML pipeline — and does not re-run `01_RawData_&_PCA`'s whole-transcriptome batch analysis.

## Pipeline context

```text
05_ML_Pipeline
      |
      v
06_Final_Panel_Validation   <- this stage (last)
```

## Folder Structure & File Reference

| File / folder | Type | Description |
|---|---|---|
| `00_expression_NP_only_locked.tsv` | Input | Expression for the 23 never-pregnant donors — 13,416 genes × 23 samples, unfiltered. |
| `00_final_filtered_expression_non_NP_donors.tsv` | Input | Expression for the 475 pregnant-cohort samples — 9,055 genes × 475 samples, post-QC. Same file used throughout `02_DGE_Analysis` onward. |
| `00_GSE108497_processed_metadata.csv` | Input | Full corrected metadata, all 512 samples (including the 23 NP donors), 29 columns. |
| `01_Panel_Validation.ipynb` | Code | Loads and combines the two expression files, fits two mixed-model specificity tests, and produces 3 diagnostic plots — 8 steps, no custom functions (`def`), matching this thesis's established code style. |
| `02_Validation_Output/` | Output (folder) | 3 CSVs (one per model, one combined) + 3 plots. See [`02_Validation_Output/README.md`](./02_Validation_Output/README.md). |

**A note on this stage's documentation:** two files were uploaded for this stage — the notebook itself, and a separate, very thorough README (`README_Panel_Specificity_Validation.md`) that the thesis author had already written, covering design decisions, rationale, and the real run's results in detail. This page is built directly from that document rather than reconstructing the same ground independently, since it's already the authoritative account of this stage.

---

## Why two expression files instead of one

There's no single unzipped file that covers all 498 samples together — the combined raw matrix only exists inside a zip archive, which this notebook avoids using. Step 1 loads the two flat files above and Step 2 restricts each to the 3 final-panel genes, transposes both from genes-as-rows to samples-as-rows, and concatenates them into one combined expression table before merging with metadata. This is a storage/upload constraint, not a methodological choice — the two files are two halves of one original matrix (their GEO sample IDs run on sequentially from one file into the other), split apart because the variance/low-expression filter was only ever computed on the non-NP training set.

## Why batch is adjusted for directly, in both models

**Complication status is confirmed to exist only for SLE donors.** Checked directly against the source paper (Hong et al., 2019, *J Exp Med*): complications (NC=46, PE=24, OC=22, summing to the 92 SLE-pregnant donors) are only ever assigned within "pregnant SLE patients." Healthy participants were pre-screened specifically to be low-risk, so complication status was never a meaningful thing to assign them. This means Model 1 (SLE-only) doesn't confound `complication_group` with `Condition` inside the model itself — `Condition` is fixed at SLE for every row there.

**What that doesn't rule out is batch.** `01_RawData_&_PCA`'s whole-transcriptome analysis already found batch comparable in size to Condition on one of its principal components — batch and biology are not fully independent in this dataset (batch 3 is 97.8% Healthy, as documented in the root README's methodological-safeguards note). If a batch happens to overlap disproportionately with, say, the PE donors specifically, an unadjusted complication test could mistake that processing difference for a complication effect — or the reverse, mask a real one.

**Folding batch into the same model as the predictor of interest** — rather than testing it separately — follows the standard approach for exactly this situation, per Leek et al. (2010, *Nat Rev Genet*) on batch effects correlating with the outcome under study, and Johnson et al. (2007, *Biostatistics*, the ComBat paper) on keeping the biological variable of interest in the same model used to estimate the batch term.

## What is a likelihood-ratio test (LRT)?

An LRT compares two *nested* models — a "full" model with the predictor of interest, and a "reduced" model without it — by asking how much better the full model explains the data: `2 × (full model's log-likelihood − reduced model's log-likelihood)` follows a chi-squared distribution under "the predictor doesn't matter," with degrees of freedom equal to how many extra parameters the full model has. Used consistently for both models, rather than mixing in a simpler Wald test for one — Wald tests are less reliable in small/moderate samples and don't extend as cleanly to a multi-level predictor (`complication_group`, 3 levels).

---

## Methodology, step by step

1. **Data input.** Loads both expression files and the metadata file, asserts all 3 panel genes are present in both expression files, prints shapes.
2. **Combine and merge.** Transposes and concatenates the two expression files, merges with metadata on sample ID, prints the donor-level 92/43/23 sanity check (SLE-pregnant / Healthy-pregnant / Healthy-never-pregnant). `pregnancy_status` is derived from `time_point == "NP"`, since the metadata carries no standalone pregnancy-status column.
3. **Complication-group derivation** (SLE donors only). Same PE → OC → NC logic used in earlier stages: PE if `pe == 1`; otherwise OC if any of `fd`, `nnd`, `pl_insuff`, `iugr`, `sga` is 1; otherwise NC. Prints the 46/24/22 sanity check.
4. **Model 1 — Complication specificity** (SLE donors only). `gene ~ complication_group + batch`, donor random intercept, LRT for the `complication_group` term, BH-FDR corrected across the 3 genes. Both models fit with `reml=False` (ML, not REML) — required for a valid LRT on fixed effects.
5. **Model 2 — Pregnancy specificity** (Healthy donors only). `gene ~ pregnancy_status + batch`, same LRT + BH-FDR approach, testing Pregnant vs. Never-Pregnant.
6. **Combined results.** Stacks both models (2 independent families, BH-FDR corrected separately across the 3 genes each — 6 tests total, no separate batch family since batch is a covariate, not a tested hypothesis), adds a per-gene pass/flag summary.
7. **Three visual counterparts**, one per grouping variable — not statistical tests, the visual complement to the LRT + FDR results above.
8. **Interpret results.** Same limitations-and-exclusions closing pattern as earlier stages' notebooks.

**How to run it:** Three plain filenames expected in the same folder as the notebook — no hardcoded absolute paths, unlike the R scripts in earlier stages. Run top to bottom in Jupyter.

**Environment:** Python. Packages: `pandas`, `numpy`, `matplotlib`, `scipy` (for `chi2`), `statsmodels` (`mixedlm`, `multipletests`). No `scikit-learn` or `umap-learn` needed here, since there's no genome-wide PCA/UMAP step in this notebook.

---

## Key Result

**All 6 FDR-corrected p-values sit above 0.05** — none of HERC5, BATF2, or SPATS2L show evidence of tracking complication severity or pregnancy status once SLE status and batch are accounted for. The panel's SLE-vs-Healthy signal isn't attributable to either confound in this cohort. See [`02_Validation_Output/README.md`](./02_Validation_Output/README.md) for the full numbers, including the one borderline result (SPATS2L, complication specificity, raw p=0.066) worth reading carefully rather than glossing over.

## Limitations

The donor random-intercept covariance structure (vs. Hong et al.'s spatial-power structure) and the LRT's chi-squared approximation are both slightly anti-conservative in general — this doesn't threaten the conclusions here specifically, since both models are built to detect a null result, and an anti-conservative test makes a true null *harder* to obtain cleanly, not easier. Batch adjustment guards against complication/pregnancy status being confounded with processing batch; it doesn't establish that `batch` fully captures every source of technical variation in the dataset.

## What's deliberately not in here

SLE vs. Healthy overall (already established in the primary DE/ML pipeline); the other 12 bootstrap-stable genes from `05_ML_Pipeline` that didn't make the final panel (out of scope here); a whole-transcriptome PCA/UMAP battery (batch is addressed inline as a covariate above instead).

## This is the final stage of the pipeline.

See the [root README](../README.md) for the overall repository structure and results summary.
