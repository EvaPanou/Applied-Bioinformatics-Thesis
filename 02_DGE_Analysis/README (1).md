# 02_DGE_Analysis

Differential gene expression (DGE) analysis stage of the SLE pregnancy biomarker
discovery project. Takes the filtered, NP-donor-corrected expression data and
identifies genes that reliably distinguish SLE from Healthy, using three
independent statistical methods.

## Folder contents

| File / folder | What it is |
|---|---|
| `00_final_filtered_expression_non_NP_donors.tsv` | Input: expression matrix (genes x samples), log2 microarray values, NP donors already excluded, low-expression genes already filtered (9,055 genes, 475 samples) |
| `00_non_NP_metadata.csv` | Input: per-sample metadata (Sample ID, Donor_id, Condition, time_point) |
| `01_limma_DEG_pipeline_v3.R` | The analysis script. Source it and call `run_deg_pipeline()` to run the full analysis end to end (see script header for details) |
| `02_Results/` | All script output — see `02_Results/README.md` for the full breakdown |

## Methodology, in brief

Three independent statistical methods are used, each testing SLE vs. Healthy
with the same significance rule (**FDR < 0.05, |logFC| > 1**, ~2-fold change),
differing only in how the effect size and repeated-measures structure
(multiple timepoints per donor) are estimated:

- **Method A** — tests each of the 5 gestational timepoints separately, then
  pools the 5 estimates via random-effects meta-analysis (`metafor::rma`)
- **Method B** — a single pooled model across all timepoints, using limma's
  `duplicateCorrelation` (one shared within-donor correlation for all genes)
- **Method C** — the same pooled question as Method B, but via a mixed model
  (`variancePartition::dream`) that estimates the within-donor correlation
  separately for each gene

A gene's inclusion in the final panel is validated by agreement across some
or all of these three independently-built methods — not by any single test
alone.

## Key result

**42 genes** are found significant by all three methods simultaneously
(`Method_All_Intersection` in `02_Results/DE_Genes/`) — a dominant,
overwhelmingly up-regulated interferon-stimulated gene (ISG) signature
(`IFI44L`, `MX1`, `OAS1/2/3`, `ISG15`, `RSAD2`, and others), confirmed stable
across gestation by a genome-wide interaction test (0 of 9,055 genes show a
significant Condition x Time interaction).

## Next stage

QA plots, followed by ML feature selection on one of the gene panels in
`02_Results/DE_Genes/` (see that folder's section in `02_Results/README.md`
for which file to use and why).
