# 03_Feature_Validation / 02_Results / PCA_Plots

## Folder Overview

PCA (PC1 vs. PC2) for all 3 gene panels — the `All_Genes` before-baseline plus both after-panels — colored by Condition, shaped by timepoint, with a 95% confidence ellipse per Condition. Produced by Step 3 of `../../01_Feature_Validation.R`.

## Folder Structure & File Reference

| File | Genes | Stage | Description |
|---|---|---|---|
| `PCA_All_Genes.png` | 9,055 | Before | Every gene, unfiltered, blind to Condition |
| `PCA_DEG_Union.png` | 44 | After (lenient) | Restricted to the `DEG_Union` panel |
| `PCA_Final_Panel.png` | 42 | After (strict) | Restricted to the `Final_Panel` |

---

## Results

**`PCA_All_Genes.png`**:

![PCA, All Genes, before selection](PCA_All_Genes.png)

The two Condition ellipses (Healthy blue, SLE orange) are large and almost fully overlapping, centered close together — the whole-transcriptome view shows no meaningful Condition-driven separation, as expected: at this scale, donor-to-donor variation, technical noise, and biology unrelated to SLE all swamp the disease signal. No obvious single-donor or batch outlier is visible either, which is a useful negative QC result in its own right — nothing here suggests a data-quality problem independent of Condition.

**`PCA_DEG_Union.png`**:

![PCA, DEG Union panel, after selection](PCA_DEG_Union.png)

A visibly different picture: the Healthy ellipse sits mostly to the left (roughly PC1 -10 to -2), the SLE ellipse to the right (roughly PC1 -2 to 15) — overlapping in the middle, but clearly offset along PC1 in a way the before plot didn't show. This is the visual counterpart of the PERMANOVA jump documented in `../PERMANOVA/README.md` — expected once the plot is restricted to genes individually selected for differing by Condition, per the caveat in the parent stage's own README.

**`PCA_Final_Panel.png`**:

![PCA, Final Panel, after selection](PCA_Final_Panel.png)

Essentially the same separation pattern as `DEG_Union`, consistent with `Final_Panel`'s slightly higher PERMANOVA R² (41.6% vs. 41.1%) despite having 2 fewer genes — removing the 2 genes that don't reach cross-method consensus doesn't cost any visible separation.

**Used downstream by:** nothing further in the pipeline — visual QC artifacts for the thesis discussion, complementing the formal PERMANOVA test in `../PERMANOVA/`.
