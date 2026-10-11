# ==========================================
# FEATURE VALIDATION -- BEFORE vs AFTER PCA + PERMANOVA
# Built on 02_DGE_Analysis's outputs
# ==========================================
#
# Volcano/MA show ONE Condition comparison (SLE vs Healthy), using
# Method_All_Union_Metrics.tsv -- a table that spans all 3 methods
# (A/B/C), not one single method's own output -- restricted to each gene
# panel. Not a separate plot per method or per timepoint. This script's
# job is to validate the DEG stage's gene selection by comparing gene
# panels, not methods.
#
# BEFORE vs. AFTER validation:
# GENE_PANELS below is structured as one "before" baseline (All_Genes --
# every gene, unfiltered, blind to Condition) and two "after" panels
# (DEG_Union, Final_Panel -- the sets 02_DGE_Analysis actually validated).
# Running PCA/PERMANOVA on both sides is standard practice for
# sanity-checking a DEG result: the full transcriptome should show
# weak/messy Condition structure (dominated by donor variation, technical
# noise, and biology unrelated to SLE), and the DEG panels should show it
# much more cleanly. That comparison IS useful, but calibrate what it
# actually proves: with ~44 genes nearly all moving in ONE coherent
# direction (the ISG signature), a clean PC1 split in the "after" panels
# is close to expected by construction -- genes individually selected
# because they differ by Condition, in a shared direction, mechanically
# create Condition-aligned covariance once you restrict to just them. So
# a split appearing isn't surprising on its own. What IS genuinely
# informative:
#   (1) whether the "before" plot reveals a DATA QUALITY problem
#       (outlier donor, batch clustering) independent of Condition --
#       that's real QC value the "after" panels can't give you;
#   (2) the *size* of the jump in PERMANOVA R2 from "before" to "after"
#       -- a quantifiable measure of how much a validated panel
#       concentrates the signal, not just a yes/no on separation;
#   (3) that the "after" R2 for Condition is computed at the DONOR
#       level (see run_permanova_condition below) -- a properly
#       controlled effect size, not a naive sample-level number that
#       repeated timepoints could inflate.

# ---- INSTALLATIONS ------------------
# =====================================
# turn into comment once you have installed each of these once

# install.packages(c("tidyverse","data.table","vegan"), dependencies = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)   # general data wrangling
  library(data.table)  # fast data import/export
  library(ggplot2)     # plotting
  library(ggrepel)     # non-overlapping gene labels on the volcano plot
  library(vegan)       # PERMANOVA (adonis2)
})


# ---- FILE PATHS --------------------------
# ==========================================
# Everything below is a plain file name, not a path. All four input files
# sit in the same folder as this script (03_Feature_Validation), so R
# finds them automatically -- no folder-building needed.

setwd("C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/03_Feature_Validation")

EXPR_FILE <- "00_final_filtered_expression_non_NP_donors.tsv" # identical to the previous step's inputs, they are from the pre-processing raw data step
META_FILE <- "00_non_NP_metadata.csv"

DEG_UNION_FILE   <- "00_Method_All_Union.tsv"         # AFTER (lenient) panel -- 44 genes, A or B or C
FINAL_PANEL_FILE <- "00_Method_All_Intersection.tsv"  # AFTER (strict)  panel -- ~42-44 genes, cross-method-validated

# Per-gene stats for Volcano/MA -- spans all 3 methods (A/B/C), not tied
# to one method's own output table. Copied locally the same way the two
# gene-list files above were.
ALL_METHODS_METRICS_FILE <- "00_Method_All_Union_Metrics.tsv"
FDR_CUT <- 0.05
LFC_CUT <- 1

OUTPUT_DIRECTORY <- "02_Results/"
PCA_DIR       <- file.path(OUTPUT_DIRECTORY, "PCA_Plots")
PERMANOVA_DIR <- file.path(OUTPUT_DIRECTORY, "PERMANOVA")
V_DIR         <- file.path(OUTPUT_DIRECTORY, "Volcano_Plots")
MA_DIR        <- file.path(OUTPUT_DIRECTORY, "MA_Plots")
dir.create(PCA_DIR,       showWarnings = FALSE, recursive = TRUE)
dir.create(PERMANOVA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(V_DIR,         showWarnings = FALSE, recursive = TRUE)
dir.create(MA_DIR,        showWarnings = FALSE, recursive = TRUE)


# ---- VARIABLE DEFINITIONS -----------------
# ===========================================
# Identical constants/levels to 02_DGE_Analysis's limma script, so both
# stages always refer to Condition/Time/Donor the same way.
SAMPLE    <- "Sample"
CONDITION <- "Condition"
TIME      <- "time_point"
SUBJECT   <- "Donor_id"

TIME_POINT_LEVELS <- c("<16 weeks", "16-23 weeks", "24-31 weeks", "32-40 weeks", "PP")

# Colorblind-safe Condition palette (orange/blue), applied to every plot
# below that colors by Condition.
CONDITION_COLORS <- c(SLE = "#E69F00", Healthy = "#0072B2")

# One "before" baseline + two "after" panels (see header comment above
# for why it's framed this way). NA = no subsetting at all (every gene).
GENE_PANELS <- list(
  All_Genes   = NA_character_,   # BEFORE -- every gene, unfiltered, blind to Condition
  DEG_Union   = DEG_UNION_FILE,  # AFTER (lenient)
  Final_Panel = FINAL_PANEL_FILE # AFTER (strict)
)

# Used to label PCA titles and the PERMANOVA summary table so the
# before/after comparison is legible without cross-referencing GENE_PANELS.
stage_for_panel <- function(panel_name) {
  if (panel_name == "All_Genes") "Before (unfiltered)" else "After (DEG-selected)"
}


# ---- LOAD DATA HELPERS ---------------------
# =============================================

load_expression <- function(path) {
  x <- data.table::fread(path) |> as.data.frame()
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  x
}

load_metadata <- function(path) {
  m <- data.table::fread(path) |> as.data.frame()
  m[[SAMPLE]] <- trimws(m[[SAMPLE]])
  m
}

# Reading a gene list TSV (the Method_*.tsv files 02_DGE_Analysis writes all
# use a "Gene" column, but this also accepts gene/SYMBOL/symbol so older
# lists still work).
get_deg_genes_from_list <- function(path, expr, subset_label = "gene panel") {
  if (!file.exists(path)) stop(paste0("File not found: ", path))
  dg <- data.table::fread(path)
  col_gene <- intersect(names(dg), c("Gene", "gene", "SYMBOL", "symbol"))
  if (length(col_gene) == 0) stop(paste0("Cannot find gene column in ", subset_label))
  genes <- unique(dg[[col_gene[1]]])
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < 2) stop(paste0("Too few genes for PCA/PERMANOVA in ", subset_label))
  genes
}


# ---- VOLCANO & MA (one Condition comparison, restricted to a gene panel) ----
# ================================================================================
# Both take a small table (just the panel's genes) with logFC, AveExpr,
# adj.P.Val columns already in it -- all genes plotted have already passed
# the DEG panel's own significance filter, so every point here is
# expected to sit past the dashed threshold lines. The plot's job is to
# show the SPREAD of effect size/significance within the panel, not to
# find new hits.

plot_volcano <- function(gene_table, outfile, plot_title) {
  gene_table$neglog10FDR <- -log10(gene_table$adj.P.Val)

  p <- ggplot(gene_table, aes(x = logFC, y = neglog10FDR)) +
    geom_point(size = 2, color = "black") +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 100) +
    labs(title = plot_title, x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_minimal(base_size = 12)

  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
}

plot_MA <- function(gene_table, outfile, plot_title) {
  p <- ggplot(gene_table, aes(x = AveExpr, y = logFC)) +
    geom_point(size = 2, color = "black") +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 100) +
    labs(title = plot_title, x = "Average Expression (log2)", y = "log2 Fold Change") +
    theme_minimal(base_size = 12)

  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
}


# ---- PCA (generic, one gene panel at a time) ----
# ===================================================

run_pca <- function(mat, meta, title, outfile) {
  mat <- as.matrix(mat)
  keep <- apply(mat, 1, function(v) sd(v, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]

  pcs <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  df <- as.data.frame(pcs$x[, 1:2, drop = FALSE])
  df[[SAMPLE]] <- rownames(df)
  df <- df |> dplyr::left_join(meta, by = SAMPLE)

  p <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[CONDITION]], shape = .data[[TIME]])) +
    geom_point(size = 2, alpha = 0.9) +
    # group= here matters: without it, ggplot draws one ellipse per
    # Condition+Timepoint combination (10 ellipses) instead of one per
    # Condition (2 ellipses), since shape is also a discrete aesthetic.
    stat_ellipse(aes(group = .data[[CONDITION]]), level = 0.95, linewidth = 0.8) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(title = title, color = CONDITION, shape = "Time point") +
    theme_minimal(base_size = 12)

  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
  invisible(list(pca = pcs, scores = df))
}


# ---- PERMANOVA (generic, one gene panel at a time) ----
# =========================================================
#
# Two separate tests, each using the permutation scheme that actually
# matches what varies at what level -- this is the same repeated-measures
# reasoning duplicateCorrelation()/dream() apply in 02_DGE_Analysis
# (Steps 2 and 6), just for a multivariate distance-based test instead of
# a per-gene linear model:
#
#  - Condition is a BETWEEN-donor factor: every sample from one donor
#    shares the same label. Testing at the sample level would let a
#    donor with 5 timepoints cast 5 "votes" for the same label -- exactly
#    the pseudoreplication problem the whole DEG pipeline exists to
#    avoid. Fixed here by collapsing each donor to ONE point (mean
#    expression across that donor's own timepoints) first, so every
#    donor contributes exactly one independent observation.
#
#  - Time is a WITHIN-donor factor: each donor genuinely contributes
#    several different Time values, so collapsing to one point per donor
#    would destroy the thing being tested. Time is tested at the full
#    sample level instead, with strata = Donor_id restricting
#    permutations to within each donor -- the nested/blocked design
#    vegan's own adonis2 documentation describes for exactly this
#    situation (a factor that varies within a blocking factor).

run_permanova_condition <- function(mat, meta, panel_label) {
  mat <- as.matrix(mat)
  keep <- apply(mat, 1, function(v) sd(v, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]

  # Keep only the metadata rows for samples that actually exist as columns
  # in mat -- the metadata file has a few extra samples that never made it
  # into the filtered expression matrix, so without this line the code
  # below tries to average in samples that don't exist and errors out.
  meta <- meta[meta[[SAMPLE]] %in% colnames(mat), , drop = FALSE]

  donor_ids <- unique(as.character(meta[[SUBJECT]]))

  # One averaged expression vector per donor (genes x donors), then
  # transposed to the donors x genes shape dist() expects (distances are
  # computed between ROWS).
  donor_level_matrix <- sapply(donor_ids, function(donor_id) {
    samples_for_this_donor <- meta[[SAMPLE]][as.character(meta[[SUBJECT]]) == donor_id]
    rowMeans(mat[, samples_for_this_donor, drop = FALSE])
  })
  donor_level_matrix <- t(donor_level_matrix)

  donor_condition <- meta[[CONDITION]][match(donor_ids, as.character(meta[[SUBJECT]]))]
  donor_metadata  <- setNames(data.frame(donor_ids, donor_condition, stringsAsFactors = FALSE),
                               c(SUBJECT, CONDITION))

  distance_matrix   <- dist(donor_level_matrix, method = "euclidean")
  formula_condition <- stats::as.formula(paste("distance_matrix ~", CONDITION))
  permanova_result  <- vegan::adonis2(formula_condition, data = donor_metadata, permutations = 999)

  data.frame(
    Panel       = panel_label,
    Stage       = stage_for_panel(panel_label),
    Test        = "Condition (donor-level, 1 point/donor)",
    N           = length(donor_ids),
    N_genes     = nrow(mat),
    R2          = permanova_result$R2[1],
    F_statistic = permanova_result$F[1],
    P_value     = permanova_result$`Pr(>F)`[1]
  )
}

run_permanova_time <- function(mat, meta, panel_label) {
  mat <- as.matrix(mat)
  keep <- apply(mat, 1, function(v) sd(v, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]

  meta_ordered <- meta[match(colnames(mat), meta[[SAMPLE]]), , drop = FALSE]

  distance_matrix <- dist(t(mat), method = "euclidean")
  formula_time    <- stats::as.formula(paste("distance_matrix ~", TIME))
  permanova_result <- vegan::adonis2(
    formula_time, data = meta_ordered,
    strata = as.character(meta_ordered[[SUBJECT]]),
    permutations = 999
  )

  data.frame(
    Panel       = panel_label,
    Stage       = stage_for_panel(panel_label),
    Test        = "Time (donor-blocked, sample-level)",
    N           = ncol(mat),
    N_genes     = nrow(mat),
    R2          = permanova_result$R2[1],
    F_statistic = permanova_result$F[1],
    P_value     = permanova_result$`Pr(>F)`[1]
  )
}


# ==========================================================================
# THE WHOLE FEATURE-VALIDATION PIPELINE, AS ONE FUNCTION
# ==========================================================================
# Same reasoning as 02_DGE_Analysis's run_deg_pipeline(): sourcing this file
# only DEFINES run_feature_validation() and loads packages/paths -- it
# doesn't run anything until run_feature_validation() is actually called
# (see the bottom of this file).

run_feature_validation <- function() {

  # =====================================================================
  # STEP 1: LOAD DATA
  # =====================================================================
  message("[1/4] Loading corrected expression + metadata (same inputs 02_DGE_Analysis used)...")
  expr <- load_expression(EXPR_FILE)
  meta <- load_metadata(META_FILE)
  meta[[CONDITION]] <- factor(meta[[CONDITION]], levels = c("Healthy", "SLE"))
  meta[[TIME]]      <- factor(meta[[TIME]], levels = TIME_POINT_LEVELS)
  meta[[SUBJECT]]   <- as.character(meta[[SUBJECT]])
  message("       ", ncol(expr), " samples, ", nrow(expr), " genes, ",
          length(unique(meta[[SUBJECT]])), " donors.")

  # =====================================================================
  # STEP 2: VOLCANO + MA -- Method_All_Union_Metrics (spans all 3 methods,
  # not one method's own table), restricted to DEG_Union and Final_Panel.
  # =====================================================================
  message("[2/4] Volcano + MA plots (Condition: SLE vs Healthy)...")
  all_methods_metrics <- data.table::fread(ALL_METHODS_METRICS_FILE)
  all_methods_metrics <- all_methods_metrics |> dplyr::rename(logFC = pooled_logFC)

  # AveExpr isn't in this file -- it only ever comes from a single
  # method's own topTable output, which is exactly what we're avoiding
  # here. Computed directly from the raw expression matrix instead, so
  # it's tied to the data, not to any one method.
  average_expression_per_gene <- rowMeans(as.matrix(expr))
  all_methods_metrics$AveExpr <- average_expression_per_gene[all_methods_metrics$Gene]

  for (panel_name in c("DEG_Union", "Final_Panel")) {
    panel_genes <- get_deg_genes_from_list(GENE_PANELS[[panel_name]], expr, panel_name)
    panel_table <- all_methods_metrics[all_methods_metrics$Gene %in% panel_genes, ]

    plot_volcano(panel_table, file.path(V_DIR, paste0("Volcano_", panel_name, ".png")),
                 paste0("Volcano -- Condition (SLE vs Healthy): ", gsub("_", " ", panel_name)))
    plot_MA(panel_table, file.path(MA_DIR, paste0("MA_", panel_name, ".png")),
            paste0("MA -- Condition (SLE vs Healthy): ", gsub("_", " ", panel_name)))
    message("       ", panel_name, ": ", nrow(panel_table), " genes plotted.")
  }

  # =====================================================================
  # STEP 3: PCA -- BEFORE (all genes) vs. AFTER (DEG_Union, Final_Panel)
  # =====================================================================
  message("[3/4] PCA plots...")
  panel_expression_by_name <- list()
  for (panel_name in names(GENE_PANELS)) {
    panel_path <- GENE_PANELS[[panel_name]]
    panel_expression_by_name[[panel_name]] <- if (is.na(panel_path)) {
      expr
    } else {
      expr[get_deg_genes_from_list(panel_path, expr, panel_name), , drop = FALSE]
    }
    run_pca(
      panel_expression_by_name[[panel_name]], meta,
      title   = paste0("PCA -- ", stage_for_panel(panel_name), ": ", gsub("_", " ", panel_name),
                        " (", nrow(panel_expression_by_name[[panel_name]]), " genes)"),
      outfile = file.path(PCA_DIR, paste0("PCA_", panel_name, ".png"))
    )
    message("       ", panel_name, " [", stage_for_panel(panel_name), "]: ",
            nrow(panel_expression_by_name[[panel_name]]), " genes plotted.")
  }

  # =====================================================================
  # STEP 4: PERMANOVA -- same before/after panels as Step 3, but as a
  # formal R2/p-value instead of a PCA plot to eyeball. Compare each
  # panel's Condition R2 against All_Genes' to see the size of the jump.
  # =====================================================================
  message("[4/4] PERMANOVA (Condition, donor-level; Time, donor-blocked)...")
  permanova_rows <- list()
  for (panel_name in names(panel_expression_by_name)) {
    panel_expr <- panel_expression_by_name[[panel_name]]
    permanova_rows[[paste0(panel_name, "_Condition")]] <- run_permanova_condition(panel_expr, meta, panel_name)
    permanova_rows[[paste0(panel_name, "_Time")]]      <- run_permanova_time(panel_expr, meta, panel_name)
  }
  permanova_summary <- do.call(rbind, permanova_rows)
  data.table::fwrite(permanova_summary, file.path(PERMANOVA_DIR, "PERMANOVA_Summary.tsv"), sep = "\t")
  message("       PERMANOVA summary saved: ", file.path(PERMANOVA_DIR, "PERMANOVA_Summary.tsv"))
  print(permanova_summary, row.names = FALSE)

  message("Done. Results saved under: ", OUTPUT_DIRECTORY)
}


# ==========================================================================
# ACTUALLY RUN IT
# ==========================================================================
run_feature_validation()
