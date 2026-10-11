# RawData_&_PCA Processing

## Folder Overview

This is the second stage of the pipeline. It takes the corrected, non-NP-split metadata from `00_Metadata` and the raw microarray expression matrix, and does two structurally different things with them: a Python notebook (`01_RawData_Processing.ipynb`) performs sample-level outlier detection and gene-level filtering on the raw data, and an R script (`03_ALASCA_Analysis.R`) runs a repeated-measures multivariate decomposition (RM-ASCA+, via the ALASCA package) on the notebook's cleaned output, to characterize how gestational time and SLE status each independently shape the transcriptome. The two tools are used together specifically because outlier/filtering logic is native to Python's data-science stack, while the donor-aware repeated-measures modeling ALASCA performs has no equivalent tool in Python.

## Folder Structure

Per the repo-wide convention, `00_` files are inputs, odd numbers are code, even numbers are outputs:

```text
01_RawData_&_PCA/
├── 00_GSE108497_processed_metadata.csv     [input]  corrected metadata, 512 samples (from 00_Metadata)
├── 00_GSE108497_raw.zip                    [input]  raw expression matrix, zipped (13,416 genes × 512 samples)
├── 00_non_NP_metadata.csv                  [input]  the non-NP analytic cohort, 489 samples
├── 01_RawData_Processing.ipynb             [code]   outlier detection + gene filtering (Python)
├── 02_Analysis_Results/                    [output] every table/figure from 01 — see its own README
├── 03_ALASCA_Analysis.R                    [code]   RM-ASCA+ repeated-measures decomposition (R)
└── 04_ALASCA_Output/                       [output] every table/figure from 03 — see its own README
```

## File Reference

| File | Type | Description |
|---|---|---|
| `00_GSE108497_processed_metadata.csv` | Input | The full corrected metadata from `00_Metadata` — 512 samples × 29 columns, NP donors still included. |
| `00_GSE108497_raw.zip` | Input | Zipped raw expression matrix (`GSE108497_raw.tsv` inside, ~130 MB uncompressed) — 13,416 genes × 512 samples. Read directly from the zip without extracting to disk. |
| `00_non_NP_metadata.csv` | Input | The non-NP analytic cohort's metadata — 489 samples × 29 columns. (Note: despite the `00_` prefix suggesting it's untouched input, this file is actually produced mid-notebook — see Methodology below — and is really a bridge between the metadata stage and this one.) |
| `01_RawData_Processing.ipynb` | Code | Loads the raw matrix and metadata, checks basic data quality, locks NP donors aside, detects and removes outlier samples by two independent methods, filters out low-variance and low-expression genes, and prepares the ALASCA-ready files. |
| `02_Analysis_Results/` | Output (folder) | Every table and figure the notebook produces — 19 files. See [`02_Analysis_Results/README.md`](./02_Analysis_Results/README.md). |
| `03_ALASCA_Analysis.R` | Code | Reads the notebook's cleaned, filtered expression matrix and fits a repeated-measures ASCA+ model (`time_point * Condition + (1|Donor_id)`), decomposing variance into three separate effects: pure time, pure Condition, and their interaction. |
| `04_ALASCA_Output/` | Output (folder) | Every plot and the saved model object from the R script — 8 files (6 real outputs plus 2 README variants). See [`04_ALASCA_Output/README.md`](./04_ALASCA_Output/README.md). |

---

## Methodology

### `01_RawData_Processing.ipynb`

**What it does**, in the order the notebook runs:

1. **Loads** the raw expression matrix straight out of `00_GSE108497_raw.zip` (no extraction to disk) and the metadata pointed to by a hand-set variable at the top of the notebook (`metadata_filename`, `run_label`) — the notebook is written to support re-running against different donor subsets, and doesn't enforce which one is loaded.
2. **Basic structure checks**: confirms every value is numeric, counts genes (13,416) and samples (512), confirms zero missing values, and checks the overall value range (~4–14) to confirm the data has already been log2-transformed and normalized upstream by ADEx — no additional normalization is applied here.
3. **Boxplot of every sample** (raw, before any cleaning) — a first visual pass across all 512 samples.
4. **NP/non-NP split**: donors with `time_point == "NP"` are locked aside — saved to their own files and never touched again in this notebook. From this point on, `run_label` switches to `"non_NP_donors"` and every subsequent output filename carries that suffix.
5. **Sample-level outlier detection, by two independent methods** — removing confirmed technical outliers measurably increases power to detect true DE genes downstream, following Kauffmann & Huber (2010):
   - **Method 1 (IQR-based):** each sample's median expression is compared to the median-of-medians; samples beyond 1.5× the IQR are flagged as outliers, and beyond 3× the IQR as *extreme* outliers.
   - **Method 2 (PCA + Mahalanobis/χ²):** every sample is standardized, PCA-reduced to 2 components, and its squared Mahalanobis distance from the sample cloud's center is computed by hand (not via a shortcut function); samples outside the 95% confidence ellipse (χ² threshold, 2 df) are flagged.
   - **Intersection:** only samples flagged as *extreme* by Method 1 **and** flagged by Method 2 are actually removed — a deliberately conservative rule requiring agreement from both independent methods.
6. **Boxplot after outlier removal** — the same plot as step 3, for a visual before/after comparison.
7. **Gene-level variance filtering**: genes in the bottom 10th percentile of variance (computed without using Condition labels, to avoid biasing the downstream DE analysis — see Bourgon, Gentleman & Huber, 2010) are removed.
8. **All-genes backup saved** at this point — variance-filtered but *not yet* expression-filtered — kept as a separate file specifically so later stages can compare filtered vs. unfiltered results side by side.
9. **Gene-level low-expression filtering**: genes with mean expression at or below the 25th percentile are removed, producing the final filtered matrix.
10. **Distribution checks** on the final matrix (histogram, skewness).
11. **Exploratory sample-level PCA**, run twice — once on the fully filtered matrix, once on the all-genes backup — colored by Condition and by Timepoint, to visually compare whether the expression filter changes the sample-level structure.
12. **Saves the ALASCA-ready metadata** (`Donor_id`, `time_point`, `Condition` only, for the surviving samples) — the expression matrix itself doesn't need re-saving, since step 9's output already serves that purpose.

**How to run it:** Requires `00_GSE108497_processed_metadata.csv` and `00_GSE108497_raw.zip` in the same directory as the notebook. Run top to bottom in Jupyter; the first code cell installs/upgrades its own dependencies via `%pip install` if they're missing. Creates the `02_Analysis_Results/` folder itself if it doesn't exist.

**Environment:** Python (Anaconda-managed local kernel, consistent with the rest of the thesis — Python 3.11.5 per the project's stated environment). Packages imported: `pandas`, `scipy` (`stats`, `chi2`), `scikit-learn` (`StandardScaler`, `PCA`), `numpy`, `matplotlib.pyplot`, `seaborn`, `zipfile`, `os`. As with the metadata-stage notebooks, exact package version numbers aren't pinned in the notebook itself — a `pip freeze` from the actual environment would be needed to record them precisely here.

**One naming note:** a markdown cell in the notebook refers to the next script as `02_ALASCA_Analysis.R`, but the actual R script in this folder is `03_ALASCA_Analysis.R` — a leftover from renumbering at some point. Not a functional problem, just worth knowing if you're cross-referencing the notebook's own comments against the folder.

### `03_ALASCA_Analysis.R`

**What it does:** Loads the notebook's final filtered expression matrix and ALASCA-ready metadata, reshapes the wide expression matrix into the long format ALASCA requires (one row per sample × gene combination), then fits:

```r
ALASCA(df = expression_long_format_reshaped,
       formula = value ~ time_point * Condition + (1 | Donor_id),
       effects = c("time_point", "Condition", "time_point:Condition"),
       participant_column = "Donor_id")
```

This explicitly requests three separate effect matrices — pure gestational time, pure Condition (SLE vs. Healthy), and their interaction — rather than relying on ALASCA's `separate_effects = TRUE` default, which (per the package's own documentation) only ever splits a two-way interaction model into two bundled matrices, not three cleanly separable ones. Treating Condition as fully separable from the interaction is specifically justified here by `02_DGE_Analysis`'s own finding that 0 of 9,055 genes show a significant Condition × Time interaction — there's essentially no real interaction effect being obscured by asking for a "pure" Condition effect.

All three effects are then plotted individually (`01_*_effect_plot.png`), and group-level PC1/PC2 scores are plotted for the Condition effect (`02_condition_only_pca_by_condition.png`) and the interaction effect (`03_interaction_pca_by_condition.png`) — see `04_ALASCA_Output/README.md` for the full worked-example explanation of what these score plots actually show and why several points overlap by construction. The fitted model object is saved as an `.rds` file for later reuse without re-fitting.

**How to run it:** Open in RStudio. **Before running**, you'll need to edit the two hardcoded absolute paths near the top of the script (`input_folder`, `output_folder`) — they currently point to a specific local machine's folder layout (`C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/01_RawData/...`) rather than a path relative to the script's own location, and that local path uses an older folder name (`01_RawData`) that predates this repo's current `01_RawData_&_PCA` naming. Update both to point at `02_Analysis_Results` (input) and `04_ALASCA_Output` (output) relative to wherever you've placed this repo. Source the whole script; it creates the output folder itself if needed (every `ggsave()` call passes `create.dir = TRUE`).

**Environment:** R, run via RStudio. Packages: `ALASCA` (installed from GitHub via `devtools::install_github("andjar/ALASCA")` — not available on CRAN), `dplyr`, `ggplot2` (**≥3.5.0 required** — the script's own comment notes that `create.dir=` inside `ggsave()` needs that version or newer; an older `ggplot2` will fail on the save step). No specific R version is pinned in the script itself.

**One limitation worth flagging:** the script's `ALASCA()` call has `validate = TRUE` and `n_validation_runs = 20` present but commented out — meaning no bootstrap/permutation-based significance testing of the ALASCA effects themselves is actually run in this version. The script's own header comment is explicit that this is deliberate and unverified: the exact behavior of the three-way `effects` argument hadn't been tested against a live R install at the time of writing, and the comments instruct running `str(alasca_model)` immediately after fitting to confirm effect indices 1/2/3 actually correspond to time/Condition/interaction as assumed, before trusting anything downstream.

---

## Results

The actual outputs of both code artifacts — all tables, all figures, and the saved ALASCA model — are documented in the two subfolder READMEs, since this top-level folder itself contains no output files directly:

- **[`02_Analysis_Results/README.md`](./02_Analysis_Results/README.md)** — outlier detection results, filtering results, and exploratory PCA from the Python notebook.
- **[`04_ALASCA_Output/README.md`](./04_ALASCA_Output/README.md)** — the three decomposed effects (time, Condition, interaction) and what they show, including a worked-example explanation of the interaction score plot.

**Feeds downstream to:** `02_Analysis_Results/15_final_filtered_expression_non_NP_donors.tsv` is the expression matrix that `02_DGE_Analysis` builds on directly.
