# ==========================================
# LIMMA TIMESERIES DEG ANALYSIS
# ==========================================

# ---- INSTALLATIONS ------------------
# =====================================
# turn into comment, onse you have installed then once

# install.packages(c("tidyverse","data.table"), dependencies = TRUE)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install(c("limma","ComplexHeatmap","circlize","UpSetR","matrixStats",
#                      "variancePartition","BiocParallel"),
#                      update = TRUE, ask = FALSE)
# install.packages("metafor", dependencies = TRUE)

suppressPackageStartupMessages({
  library(limma)             # linear modeling and eBayes
  library(ComplexHeatmap)    # advanced heatmaps
  library(circlize)          # color mapping for heatmaps
  library(UpSetR)            # intersection visualization
  library(tidyverse)         # general data wrangling
  library(data.table)        # fast data import/export
  library(matrixStats)       # row/column variance calculations
  library(metafor)           # random-effects meta-analysis across timepoints
  library(variancePartition) # dream(), gene-specific mixed model
  library(BiocParallel)      # parallel backend dream() uses internally
})

# dream() fits one mixed model per gene. Two options for HOW it does this:
#
# SAFE DEFAULT (below, active): SerialParam -- one gene at a time, one CPU
# core, but with a visible progress bar so you can actually watch it move
# and know it's alive, instead of staring at a silent console for hours
# with no way to tell if it's working or stuck.
BiocParallel::register(BiocParallel::SerialParam(progressbar = TRUE))
#
# FASTER BUT RISKIER (commented out): SnowParam splits the work across
# several CPU cores at once -- meaningfully faster when it works, but on
# some Windows setups (especially inside RStudio rather than a plain
# terminal) it can hang indefinitely with zero warning and zero visible
# progress, which is what very likely happened on the run that sat frozen
# for 4 hours. Only try this AFTER a SerialParam run has completed
# successfully at least once, so you know the rest of the script works and
# any future problem is isolated to parallelization specifically:
# BiocParallel::register(BiocParallel::SnowParam(max(1, parallel::detectCores() - 1), progressbar = TRUE))


# ---- FILE PATHS --------------------------
# ==========================================

setwd("C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/02_DGE_Analysis")

EXPR_FILE <- "C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/02_DGE_Analysis/00_final_filtered_expression_non_NP_donors.tsv"
META_FILE <- "C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/02_DGE_Analysis/00_non_NP_metadata.csv" 

OUTPUT_DIRECTORY <- "02_Results/"

# Setting Program B: the pooled Condition-only model (Time marginalized out,still blocked by Donor)
RUN_B_POOLED_MODEL <- TRUE

# Program C answers the SAME pooled question as Program B (Condition only,time removed), 
# but via dream() instead of duplicateCorrelation, doing the fitting per gene not per whole genome (within each donor)
RUN_C_DREAM_MODEL <- TRUE


# ---- VARIABLE DEFINITIONS -----------------
# ===========================================
SAMPLE    <- "Sample"
CONDITION <- "Condition"   # Healthy / SLE
TIME      <- "time_point"  # PP, <16 weeks, etc
SUBJECT   <- "Donor_id"

TIME_POINT_LEVELS <- c("<16 weeks", "16-23 weeks", "24-31 weeks", "32-40 weeks", "PP")

FDR_thresholds <- c(0.05, 0.01)
LOG_FOLD_CHANGE_UP_THRESHOLD   <-  1
LOG_FOLD_CHANGE_DOWN_THRESHOLD <- -1

dir.create(OUTPUT_DIRECTORY, showWarnings = FALSE, recursive = TRUE)


# ==========================================================================
# THE WHOLE PIPELINE, AS ONE FUNCTION
# ==========================================================================
# Everything from Step 1 to Step 13 is now wrapped inside this one function,
# run_deg_pipeline(). Nothing below this point has changed -- same steps,
# same comments, same order, same variable names, same indentation. R
# doesn't require code inside a function to be indented differently from
# code outside one (unlike Python) -- wrapping in { } is the only thing
# that matters for it to work, so nothing else needed touching. Source this
# whole file once (which installs/loads packages, sets your file paths, and
# DEFINES this function, but doesn't run the analysis yet), then call
# run_deg_pipeline() to actually run all 13 steps in one go and wait for it
# to finish (or error).
#
# One real difference worth knowing: everything Steps 1-13 create (expression,
# metadata, fitted_model_A, meta_analysis_results_A, etc.) now lives INSIDE
# this function while it runs, and disappears once it finishes -- you won't
# see those objects sitting in your R environment afterward the way you
# would running the steps directly. The printed messages, saved output
# files, and plots are unaffected -- this only changes whether you can
# inspect intermediate R objects in the console after a run completes.

run_deg_pipeline <- function() {

# =====================================================================
# STEP 1: LOADING INPUTS
# =====================================================================

# Every step here exists to prevent a SILENT sample mismatch, so that the script doesn't compute on misaligned data. 

expression <- data.table::fread(EXPR_FILE) |> as.data.frame() # reading expression data as a dataframe
rownames(expression) <- expression[[1]]
expression[[1]] <- NULL

metadata <- data.table::fread(META_FILE) |> as.data.frame() # reading metadata as a dataframe
metadata[[SAMPLE]] <- trimws(metadata[[SAMPLE]])
colnames(expression) <- trimws(colnames(expression))

common <- intersect(colnames(expression), metadata[[SAMPLE]]) # checking if they both have the same samples, and that nothing is missing
expression <- expression[, common, drop = FALSE]
metadata <- metadata[metadata[[SAMPLE]] %in% common, , drop = FALSE]
metadata <- metadata[match(colnames(expression), metadata[[SAMPLE]]), , drop = FALSE]

metadata[[CONDITION]] <- factor(metadata[[CONDITION]], levels = c("Healthy","SLE")) # verifying donor info, condition and timepoint levels
metadata[[TIME]] <- factor(
  metadata[[TIME]],
  levels = TIME_POINT_LEVELS
)
metadata[[SUBJECT]] <- factor(metadata[[SUBJECT]])

message("[1/13] Inputs loaded: ", ncol(expression), " samples, ", nrow(expression), " genes.")

# =====================================================================
# STEP 2: METHOD A - BUILDING DESIGN MATRIX AND FIT BASE MODEL 
# =====================================================================

# Here a SAMPLE-means design (~0 + Condition:Time) is set: each of the 10 "Condition (N=2) x Timepoint (N=5)" combinations 
# gets its own column, with no shared "baseline" the way an intercept model would have. This makes each column's coefficient 
# directly interpretable as that sample's mean expression.
#
# duplicateCorrelation() + block=Donor_id exists because samples from the same donor across different timepoints are NOT statistically independent, 
# as they share that donor's baseline biology, genetics, and technical batch quirks. Treating them as independent would understate the true
# uncertainty and inflate false positives (pseudoreplication). This estimates ONE consensus correlation value across the whole genome and
# feeds it into lmFit() as a generalized-least-squares correction (Smyth, Michaud & Scott, 2005)

design_matrix_A <- model.matrix(~ 0 + metadata[[CONDITION]]:metadata[[TIME]])
colnames(design_matrix_A) <- make.names(colnames(design_matrix_A))

donor_correlation_estimate_A <- duplicateCorrelation(expression, design_matrix_A, block = metadata[[SUBJECT]])
base_model_fit_A   <- lmFit(expression, design_matrix_A, block = metadata[[SUBJECT]], correlation = donor_correlation_estimate_A$consensus.correlation)
message("[2/13] Within-subject correlation (base model): ", round(donor_correlation_estimate_A$consensus.correlation, 3))

# Building a design_matrix, a plain named list, indexed directly wherever needed below. It is needed because
# model.matrix()'s auto-generated column names aren't readable strings
condition_timepoint_to_column_map_A <- list()
for (column_index in seq_len(ncol(design_matrix_A))) {
  idx <- which(design_matrix_A[, column_index] != 0)
  cond <- unique(as.character(metadata[[CONDITION]][idx]))
  tp   <- unique(as.character(metadata[[TIME]][idx]))
  if (length(cond) != 1 || length(tp) != 1) {
    stop("Design column ", colnames(design_matrix_A)[column_index],
         " does not map uniquely to a single (Condition, Timepoint).")
  }
  key <- paste(cond, tp, sep = "||")
  condition_timepoint_to_column_map_A[[key]] <- colnames(design_matrix_A)[column_index]
}


# =====================================================================
# STEP 3: METHOD A - Comparing SLE vs Healthy AT EACH TIMEPOINT
# =====================================================================

# METHOD A tests the SLE-vs-Healthy Condition effect separately within each of the 5 gestational timepoints, using limma's
# duplicateCorrelation to correct for repeated donor samples (one genome-wide correlation value applied to every gene). This follows the
# "conditional DEG analysis" axis used across the project's literature review (TiSA/Lefol et al. 2023; Chowdhury et al. 2023). 
# The practical payoff: running this 5 times, to directly see whether the SLE-vs-Healthy signal is stable across gestation 
# or concentrated at one stage, rather than logically assuming stability.

contrasts_A <- list()
for (tp in TIME_POINT_LEVELS) {
  sle_column <- condition_timepoint_to_column_map_A[[paste("SLE", tp, sep = "||")]]
  healthy_column  <- condition_timepoint_to_column_map_A[[paste("Healthy", tp, sep = "||")]]
  contrasts_A[[paste0("SLE_vs_Healthy_at_", tp)]] <- paste0(sle_column, " - ", healthy_column)
}

contrast_matrix_A   <- makeContrasts(contrasts = unlist(contrasts_A), levels = design_matrix_A)
fitted_model_A <- contrasts.fit(base_model_fit_A, contrast_matrix_A) |> eBayes()

message("[3/13] Method A fit: ", ncol(fitted_model_A), " timepoint contrasts.")


# ================================================================================
# STEP 4: METHOD A - Sving the DE TABLES, and collecting UP/DOWN SETS AT FDR 0.05
# ================================================================================

# The FULL table (every gene, unfiltered, ranked by p-value) is saved alongside the filtered UP/DOWN tables because a hard significance
# cutoff hides genes that are close but don't quite clear it, and the full ranked table is what lets a borderline case be inspected 
# directly rather than assumed irrelevant.

output_directory_A <- file.path(OUTPUT_DIRECTORY, "A_SLE_vs_Healthy")
dir.create(output_directory_A, showWarnings = FALSE, recursive = TRUE)

significant_upregulated_genes_by_timepoint_A   <- list()
significant_downregulated_genes_by_timepoint_A <- list()

for (coefficient_index in seq_len(ncol(fitted_model_A))) {
  contrast_title <- colnames(fitted_model_A)[coefficient_index]
  top_table <- topTable(fitted_model_A, coef = coefficient_index, number = Inf, sort.by = "P")
  top_table$Gene <- rownames(top_table)

  # full, unfiltered table
  fwrite(top_table, file.path(output_directory_A, paste0("A_FULL_", make.names(contrast_title), ".tsv")), sep = "\t")

  # UP/DOWN at FDR 0.05 -- the 0.01 tables were dropped: nothing downstream
  # in this script ever reads them (the final gene lists all come from the
  # meta-analysis's own 0.05 threshold applied to pooled statistics, not
  # from these raw per-contrast tables), so they were pure unused reference
  # material. FDR_thresholds is kept as a named constant below in case you
  # want a stricter table again later -- just loop over it here as before.
  upregulated_subset   <- subset(top_table, adj.P.Val < 0.05 & logFC >  LOG_FOLD_CHANGE_UP_THRESHOLD)
  downregulated_subset <- subset(top_table, adj.P.Val < 0.05 & logFC <  LOG_FOLD_CHANGE_DOWN_THRESHOLD)
  fwrite(upregulated_subset,   file.path(output_directory_A, paste0("A_UP_FDR0.05_logFC", LOG_FOLD_CHANGE_UP_THRESHOLD,  "_", make.names(contrast_title), ".tsv")), sep = "\t")
  fwrite(downregulated_subset, file.path(output_directory_A, paste0("A_DN_FDR0.05_logFC", LOG_FOLD_CHANGE_DOWN_THRESHOLD, "_", make.names(contrast_title), ".tsv")), sep = "\t")

  # collect sets at FDR 0.05 specifically, for the heatmap/UpSet steps below
  significant_upregulated_genes_by_timepoint_A[[contrast_title]]   <- subset(top_table, adj.P.Val < 0.05 & logFC >  LOG_FOLD_CHANGE_UP_THRESHOLD)$Gene
  significant_downregulated_genes_by_timepoint_A[[contrast_title]] <- subset(top_table, adj.P.Val < 0.05 & logFC <  LOG_FOLD_CHANGE_DOWN_THRESHOLD)$Gene
}

message("[4/13] Program A tables saved to: ", output_directory_A)
message("       A-UP sizes: ",   paste(sprintf("%s=%d", names(significant_upregulated_genes_by_timepoint_A),   lengths(significant_upregulated_genes_by_timepoint_A)),   collapse = " | "))
message("       A-DOWN sizes: ", paste(sprintf("%s=%d", names(significant_downregulated_genes_by_timepoint_A), lengths(significant_downregulated_genes_by_timepoint_A)), collapse = " | "))


# ================================================================================
# STEP 5: METHOD B - POOLED SLE vs Healthy Comparison, using duplicateCorrelation
# ================================================================================
# METHOD B, using limma's duplicateCorrelation. This is NOT a simple average of expression
# values: every sample becomes one data point in a generalized-least-
# squares fit (not an arithmetic mean), with donor-blocking correcting the
# STANDARD ERRORS for repeated measures.
#
# duplicateCorrelation estimates ONE consensus correlation value, across the WHOLE GENOME and applies it to every gene equally, 
# regardless of whether that specific gene's real within-donor correlation is higher or lower than average. It also corrects the 
# VARIANCE estimate for repeated donor samples, but does not rebalance each donor's INFLUENCE on the group mean. 
# So, a donor with 5 timepoints still pulls 5x harder on their Condition's average than a donor with 1 sample would. 

significant_upregulated_genes_B   <- list()
significant_downregulated_genes_B <- list()
fitted_model_B <- NULL
full_results_B <- NULL

if (isTRUE(RUN_B_POOLED_MODEL)) {

  design_matrix_B <- model.matrix(~ 0 + metadata[[CONDITION]])
  colnames(design_matrix_B) <- levels(metadata[[CONDITION]])

  donor_correlation_estimate_B <- duplicateCorrelation(expression, design_matrix_B, block = metadata[[SUBJECT]])
  message("[5/13] Pooled model (Method B) within-subject correlation: ",
          round(donor_correlation_estimate_B$consensus.correlation, 3))

  base_model_fit_B <- lmFit(expression, design_matrix_B, block = metadata[[SUBJECT]],
                       correlation = donor_correlation_estimate_B$consensus.correlation)

  contrast_matrix_B  <- makeContrasts(SLE_vs_Healthy_Pooled = SLE - Healthy, levels = design_matrix_B)
  fitted_model_B <- contrasts.fit(base_model_fit_B, contrast_matrix_B) |> eBayes()

  output_directory_B <- file.path(OUTPUT_DIRECTORY, "B_Pooled_SLE_vs_Healthy")
  dir.create(output_directory_B, showWarnings = FALSE, recursive = TRUE)

  for (coefficient_index in seq_len(ncol(fitted_model_B))) {
    contrast_title <- colnames(fitted_model_B)[coefficient_index]
    top_table <- topTable(fitted_model_B, coef = coefficient_index, number = Inf, sort.by = "P", confint = TRUE)
    top_table$Gene <- rownames(top_table)
    full_results_B <- top_table   # only one coefficient in B, so this is the whole model

    fwrite(top_table, file.path(output_directory_B, paste0("B_FULL_", make.names(contrast_title), ".tsv")), sep = "\t")

    upregulated_subset   <- subset(top_table, adj.P.Val < 0.05 & logFC >  LOG_FOLD_CHANGE_UP_THRESHOLD)
    downregulated_subset <- subset(top_table, adj.P.Val < 0.05 & logFC <  LOG_FOLD_CHANGE_DOWN_THRESHOLD)
    fwrite(upregulated_subset,   file.path(output_directory_B, paste0("B_UP_FDR0.05_logFC", LOG_FOLD_CHANGE_UP_THRESHOLD,  "_", make.names(contrast_title), ".tsv")), sep = "\t")
    fwrite(downregulated_subset, file.path(output_directory_B, paste0("B_DN_FDR0.05_logFC", LOG_FOLD_CHANGE_DOWN_THRESHOLD, "_", make.names(contrast_title), ".tsv")), sep = "\t")

    significant_upregulated_genes_B[[contrast_title]]   <- subset(top_table, adj.P.Val < 0.05 & logFC >  LOG_FOLD_CHANGE_UP_THRESHOLD)$Gene
    significant_downregulated_genes_B[[contrast_title]] <- subset(top_table, adj.P.Val < 0.05 & logFC <  LOG_FOLD_CHANGE_DOWN_THRESHOLD)$Gene
  }

  message("       Program B tables saved to: ", output_directory_B)
  message("       B-UP sizes: ",   paste(sprintf("%s=%d", names(significant_upregulated_genes_B),   lengths(significant_upregulated_genes_B)),   collapse = " | "))
  message("       B-DOWN sizes: ", paste(sprintf("%s=%d", names(significant_downregulated_genes_B), lengths(significant_downregulated_genes_B)), collapse = " | "))
}


# =====================================================================
# STEP 6: METHOD C - POOLED SLE vs Healthy, using dream
# =====================================================================
# METHOD C, answers the EXACT SAME pooled question as Method B, but it estimates the within-donor correlation SEPARATELY FOR EACH GENE via 
# a true linear mixed model (REML), instead of collapsing every gene onto one shared genome-wide number the way duplicateCorrelation does. 

# dream() is explicitly built on top of limma, but with variancePartition (Hoffman & Roussos 2021) instead.
# Because a mixed model estimates an extra variance parameter per gene, it also adjusts the degrees of freedom used for significance testing

significant_upregulated_genes_C   <- list()
significant_downregulated_genes_C <- list()
fitted_model_C <- NULL
full_results_C <- NULL

if (isTRUE(RUN_C_DREAM_MODEL)) {

  # No-intercept, cell-means formula -- same style as Method A and B's
  # design matrices, so this is a genuine two-sided contrast (SLE minus
  # Healthy) instead of naming only one side.
  model_formula_C <- ~ 0 + Condition + (1 | Donor_id)
  metadata_dataframe_C <- data.frame(
    Condition = metadata[[CONDITION]],   # factor, Healthy = reference level
    Donor_id  = metadata[[SUBJECT]]
  )
  # Explicit row names (matching expression's sample names) so dream can
  # actually verify sample alignment instead of assuming it from position.
  rownames(metadata_dataframe_C) <- metadata[[SAMPLE]]

  contrast_matrix_C <- variancePartition::makeContrastsDream(
    model_formula_C, metadata_dataframe_C,
    contrasts = c(SLE_vs_Healthy_Dream = "ConditionSLE - ConditionHealthy")
  )

  message("[6/13] Fitting Method C (dream) -- one mixed model per gene, this takes longer than Methods A/B ...")
  # BPPARAM here is written out explicitly (rather than only relying on the
  # earlier register() call up top) so it's obvious right at this line what
  # controls speed/parallelism, without having to remember a setting from
  # 300 lines earlier. bpparam() just returns whatever was most recently
  # registered -- SerialParam with a progress bar, by default in this file.
  mixed_model_fit_C <- variancePartition::dream(expression, model_formula_C, metadata_dataframe_C, L = contrast_matrix_C,
                                                  BPPARAM = BiocParallel::bpparam())

  # dream() can genuinely fail to converge for individual genes (singular
  # fits surface this instead of letting affected genes silently vanish from downstream results.
  dream_gene_errors <- attr(mixed_model_fit_C, "errors")
  if (!is.null(dream_gene_errors) && length(dream_gene_errors) > 0) {
    message("       WARNING: dream() failed to converge for ", length(dream_gene_errors), " gene(s): ",
            paste(head(names(dream_gene_errors), 10), collapse = ", "),
            if (length(dream_gene_errors) > 10) ", ..." else "")
    data.table::fwrite(
      data.frame(Gene = names(dream_gene_errors), Error = unlist(dream_gene_errors)),
      file.path(OUTPUT_DIRECTORY, "MethodC_dream_convergence_failures.tsv"), sep = "\t"
    )
  } else {
    message("       dream() converged for all genes.")
  }

  fitted_model_C    <- eBayes(mixed_model_fit_C)

  output_directory_C <- file.path(OUTPUT_DIRECTORY, "C_Dream_SLE_vs_Healthy")
  dir.create(output_directory_C, showWarnings = FALSE, recursive = TRUE)

  # BUG FIX: dream()'s L= argument keeps the raw cell-means coefficients
  # (ConditionHealthy, ConditionSLE) alongside the requested contrast
  # (SLE_vs_Healthy_Dream) rather than replacing them -- fitted_model_C
  # ends up with 3 columns, not 1. Looping over every column like Methods
  # A and B do is WRONG here: "ConditionHealthy"/"ConditionSLE" each test
  # whether that group's raw average log-expression differs from zero,
  # which is true for almost every gene on a log scale and is not a
  # meaningful comparison. The earlier version of this loop unioned
  # significant genes from all 3 columns, which silently inflated Method
  # C's "significant" gene list to nearly the entire dataset (confirmed:
  # ~9,055 of 9,055 genes). Only the real contrast column is used now.
  dream_contrast_columns <- which(colnames(fitted_model_C) == "SLE_vs_Healthy_Dream")

  for (coefficient_index in dream_contrast_columns) {
    contrast_title <- colnames(fitted_model_C)[coefficient_index]
    top_table <- topTable(fitted_model_C, coef = coefficient_index, number = Inf, sort.by = "P", confint = TRUE)
    top_table$Gene <- rownames(top_table)
    full_results_C <- top_table   # only one coefficient in C, so this is the whole model

    fwrite(top_table, file.path(output_directory_C, paste0("C_FULL_", make.names(contrast_title), ".tsv")), sep = "\t")

    upregulated_subset   <- subset(top_table, adj.P.Val < 0.05 & logFC >  LOG_FOLD_CHANGE_UP_THRESHOLD)
    downregulated_subset <- subset(top_table, adj.P.Val < 0.05 & logFC <  LOG_FOLD_CHANGE_DOWN_THRESHOLD)
    fwrite(upregulated_subset,   file.path(output_directory_C, paste0("C_UP_FDR0.05_logFC", LOG_FOLD_CHANGE_UP_THRESHOLD,  "_", make.names(contrast_title), ".tsv")), sep = "\t")
    fwrite(downregulated_subset, file.path(output_directory_C, paste0("C_DN_FDR0.05_logFC", LOG_FOLD_CHANGE_DOWN_THRESHOLD, "_", make.names(contrast_title), ".tsv")), sep = "\t")

    significant_upregulated_genes_C[[contrast_title]]   <- subset(top_table, adj.P.Val < 0.05 & logFC >  LOG_FOLD_CHANGE_UP_THRESHOLD)$Gene
    significant_downregulated_genes_C[[contrast_title]] <- subset(top_table, adj.P.Val < 0.05 & logFC <  LOG_FOLD_CHANGE_DOWN_THRESHOLD)$Gene
  }

  message("       Program C tables saved to: ", output_directory_C)
  message("       C-UP sizes: ",   paste(sprintf("%s=%d", names(significant_upregulated_genes_C),   lengths(significant_upregulated_genes_C)),   collapse = " | "))
  message("       C-DOWN sizes: ", paste(sprintf("%s=%d", names(significant_downregulated_genes_C), lengths(significant_downregulated_genes_C)), collapse = " | "))
}


# =====================================================================
# STEP 7: HEATMAP (Method A - DEG union)
# =====================================================================

# Union (not intersection) of UP+DOWN across all 5 Program A contrasts:
# this heatmap is a descriptive overview of everything Program A flagged
# anywhere, not a filtered "final" gene set.

heatmap_output_file <- file.path(OUTPUT_DIRECTORY, "Heatmap_A_DEG_Union_FDR0.05.png")
heatmap_gene_list <- unique(c(unlist(significant_upregulated_genes_by_timepoint_A), unlist(significant_downregulated_genes_by_timepoint_A)))
heatmap_gene_list <- heatmap_gene_list[heatmap_gene_list %in% rownames(expression)]

if (length(heatmap_gene_list) < 2) {
  message("[7/13] No genes for heatmap -- skipping: ", heatmap_output_file)
} else {

  # Z-scoring per gene (not raw expression) is what makes genes with very different absolute expression levels 
  # visually comparable on the same color scale
  raw_expression_matrix  <- as.matrix(expression[heatmap_gene_list, , drop = FALSE])
  z_scored_expression_matrix <- t(scale(t(raw_expression_matrix)))
  z_scored_expression_matrix[is.na(z_scored_expression_matrix)] <- 0

  # Column order is set explicitly (SLE block, then Healthy block; 5 timepoints ascending within each) 
  # rather than left to automatic clustering, so the temporal structure stays visually readable
  metadata_for_heatmap <- metadata
  metadata_for_heatmap[[TIME]]      <- factor(metadata_for_heatmap[[TIME]], levels = TIME_POINT_LEVELS)
  metadata_for_heatmap[[CONDITION]] <- factor(metadata_for_heatmap[[CONDITION]], levels = c("SLE","Healthy"))

  samples  <- colnames(z_scored_expression_matrix)
  metadata_ordered_for_heatmap <- metadata_for_heatmap[match(samples, metadata_for_heatmap[[SAMPLE]]), , drop = FALSE]
  column_order_index      <- order(metadata_ordered_for_heatmap[[CONDITION]], metadata_ordered_for_heatmap[[TIME]])
  z_scored_expression_matrix     <- z_scored_expression_matrix[, column_order_index, drop = FALSE]
  metadata_ordered_for_heatmap <- metadata_ordered_for_heatmap[column_order_index, , drop = FALSE]

  column_split_levels <- c(paste("SLE", TIME_POINT_LEVELS), paste("Healthy", TIME_POINT_LEVELS))
  column_split_groups <- factor(paste(metadata_ordered_for_heatmap[[CONDITION]], metadata_ordered_for_heatmap[[TIME]]), levels = column_split_levels)

  # Row clusters: k-means with a fixed seed, then label each cluster by comparing its mean expression in SLE columns specifically vs Healthy columns 
  # specifically. This is a deliberate fix for a real bug found earlier in this project: labeling by a cluster's raw row mean doesn't work, 
  # because z-scoring forces every individual gene's OWN row mean to ~0 by mathematical construction
  set.seed(7)
  kmeans_result <- stats::kmeans(z_scored_expression_matrix, centers = 3, nstart = 25)
  cluster_assignments <- kmeans_result$cluster

  sle_column_indices <- which(metadata_ordered_for_heatmap[[CONDITION]] == "SLE")
  healthy_column_indices  <- which(metadata_ordered_for_heatmap[[CONDITION]] == "Healthy")

  cluster_sle_minus_healthy_difference <- vapply(sort(unique(cluster_assignments)), function(cluster_id) {
    rows_in_this_cluster <- which(cluster_assignments == cluster_id)
    mean(z_scored_expression_matrix[rows_in_this_cluster, sle_column_indices, drop = FALSE]) - mean(z_scored_expression_matrix[rows_in_this_cluster, healthy_column_indices, drop = FALSE])
  }, numeric(1))
  names(cluster_sle_minus_healthy_difference) <- sort(unique(cluster_assignments))

  cluster_labels <- character(length(cluster_sle_minus_healthy_difference))
  names(cluster_labels) <- names(cluster_sle_minus_healthy_difference)
  for (cluster_id in names(cluster_sle_minus_healthy_difference)) {
    difference_value <- cluster_sle_minus_healthy_difference[[cluster_id]]
    cluster_labels[cluster_id] <- if (difference_value > 0.15) "Up-regulated in SLE"
                    else if (difference_value < -0.15) "Down-regulated in SLE"
                    else "Mixed/time-dependent genes"
  }
  row_split_groups <- factor(cluster_labels[as.character(cluster_assignments)],
                      levels = c("Up-regulated in SLE",
                                 "Mixed/time-dependent genes",
                                 "Down-regulated in SLE"))

  # Write one gene-list TSV per cluster (raw output for manual review / individual-gene annotation later).
  cluster_output_directory <- file.path(OUTPUT_DIRECTORY, "Gene_Clusters")
  dir.create(cluster_output_directory, showWarnings = FALSE, recursive = TRUE)
  for (cluster_label in levels(row_split_groups)) {
    filesystem_safe_cluster_label <- gsub("[^A-Za-z0-9]+", "_", cluster_label)
    genes_in_cluster <- rownames(z_scored_expression_matrix)[row_split_groups == cluster_label]
    if (length(genes_in_cluster) == 0) next
    write.table(data.frame(Gene = genes_in_cluster),
                file.path(cluster_output_directory, paste0("Heatmap_Cluster_", filesystem_safe_cluster_label, "_genes.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  condition_colors <- c(SLE = "#d73027", Healthy = "#1a9850")
  timepoint_colors   <- setNames(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"), TIME_POINT_LEVELS)

  annotation_dataframe <- data.frame(
    Condition  = droplevels(metadata_ordered_for_heatmap[[CONDITION]]),
    time_point = droplevels(metadata_ordered_for_heatmap[[TIME]])
  )
  rownames(annotation_dataframe) <- metadata_ordered_for_heatmap[[SAMPLE]]

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df  = annotation_dataframe,
    col = list(Condition = condition_colors, time_point = timepoint_colors),
    annotation_name_side = "left",
    simple_anno_size = grid::unit(3, "mm")
  )

  color_mapping_function <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  # ward.D2 for row clustering: a variance-minimizing hierarchical method that tends to produce compact, 
  # visually coherent groupings for expression heatmaps
  png(heatmap_output_file, width = 1800, height = 1200, res = 170)
  heatmap_object <- ComplexHeatmap::Heatmap(
    z_scored_expression_matrix,
    name = "Expression (z)",
    col  = color_mapping_function,
    top_annotation = top_annotation,
    column_split = column_split_groups,
    gap = grid::unit(1.5, "mm"),
    row_split = row_split_groups,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    clustering_method_rows = "ward.D2",
    column_title = "Columns = SAMPLES (SLE \u2192 Healthy; timepoints left\u2192right)",
    row_title    = "Rows = GENES (z-scored per gene)",
    heatmap_legend_param = list(
      at = c(-2, 0, 2),
      labels = c("Down-regulated", "Mean", "Up-regulated")
    )
  )
  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
  dev.off()

  message("[7/13] Heatmap saved: ", heatmap_output_file,
          "  (genes: ", nrow(z_scored_expression_matrix), ", samples: ", ncol(z_scored_expression_matrix), ")")
}


# =====================================================================
# STEP 8: UPSET PLOTS for Method A Timepoint contrasts
# =====================================================================

# -- UP genes --
nonempty_upregulated_sets <- significant_upregulated_genes_by_timepoint_A[lengths(significant_upregulated_genes_by_timepoint_A) > 0]
if (length(nonempty_upregulated_sets) >= 2) {
  png(file.path(OUTPUT_DIRECTORY, "UpSet_A_UP_FDR0.05.png"), width = 1400, height = 900, res = 150)
  upset_plot_object <- UpSetR::upset(
    UpSetR::fromList(nonempty_upregulated_sets),
    order.by = "freq",
    mainbar.y.label = "Intersections (A-UP)",
    sets.x.label    = "Set size",
    nsets           = min(length(nonempty_upregulated_sets), 15),
    nintersects     = 30
  )
  print(upset_plot_object)
  dev.off()
} else {
  message("[8/13] Skip UpSet A-UP (need >=2 non-empty sets).")
}

# -- DOWN genes --
nonempty_downregulated_sets <- significant_downregulated_genes_by_timepoint_A[lengths(significant_downregulated_genes_by_timepoint_A) > 0]
if (length(nonempty_downregulated_sets) >= 2) {
  png(file.path(OUTPUT_DIRECTORY, "UpSet_A_DN_FDR0.05.png"), width = 1400, height = 900, res = 150)
  upset_plot_object <- UpSetR::upset(
    UpSetR::fromList(nonempty_downregulated_sets),
    order.by = "freq",
    mainbar.y.label = "Intersections (A-DOWN)",
    sets.x.label    = "Set size",
    nsets           = min(length(nonempty_downregulated_sets), 15),
    nintersects     = 30
  )
  print(upset_plot_object)
  dev.off()
} else {
  message("[8/13] Skip UpSet A-DOWN (need >=2 non-empty sets).")
}

message("[8/13] UpSet plots done.")


# =====================================================================
# STEP 9: META-ANALYSIS - pool Method A timepoints per gene
# =====================================================================
# This replaces an earlier, simpler approach (taking the union of genes significant at >=1 of the 5 timepoints) that turned out to have a real
# statistical weakness: reducing each timepoint to a binary "significant or not" call throws away the actual effect size and precision of every
# estimate. This kind of binary reduction is a documented weak method in the meta-analysis literature ("vote counting") criticized specifically 
# for ignoring sample size and effect magnitude and having low statistical power relative to methods that use the real numbers (Light &
# Smith 1971; Hedges & Olkin 1980).
#
# Instead: each of Program A's 5 timepoint contrasts is treated as a "study" estimating the same underlying SLE-vs-Healthy effect, and their
# logFC + SE per gene are combined using a random-effects model (REML), inverse-variance weighted. Meaning, a timepoint with more samples (a more
# precise estimate) contributes more to the pooled effect, exactly as standard in gene-expression meta-analysis (Choi et al. 2003; Ramasamy et al. 2008; Viechtbauer 2010).
#
# RANDOM-effects specifically (not fixed-effects): a fixed-effects model assumes all 5 timepoints are noisy measurements of one single true
# number. Random-effects instead allows the real biological effect to genuinely differ a bit across gestational stages
# 
# I2/QEp quantify exactly how much heterogeneity exists between timepoints for a given gene: a strong pooled effect AND low heterogeneity (QEp >
# 0.05, i.e. no significant inconsistency detected by Cochran's Q test) is a good time-invariant candidate for a diagnostic marker meant to work
# regardless of when in pregnancy blood is drawn; a gene that's only "significant" because of noisy averaging across genuinely inconsistent
# timepoints will show up with high QEp/I2 instead.

log_fold_change_matrix <- fitted_model_A$coefficients
t_statistic_matrix     <- fitted_model_A$t
standard_error_matrix    <- log_fold_change_matrix / t_statistic_matrix
gene_names_A   <- rownames(fitted_model_A)

meta_analysis_results_A <- data.frame(
  Gene = gene_names_A, pooled_logFC = NA_real_, pooled_pval = NA_real_,
  I2 = NA_real_, tau2 = NA_real_, QEp = NA_real_
)

for (gene_index in seq_along(gene_names_A)) {
  yi  <- log_fold_change_matrix[gene_index, ]
  sei <- standard_error_matrix[gene_index, ]
  valid_data_mask  <- is.finite(yi) & is.finite(sei) & sei > 0
  if (sum(valid_data_mask) < 2) next
  meta_analysis_model_for_gene <- tryCatch(metafor::rma(yi = yi[valid_data_mask], sei = sei[valid_data_mask], method = "REML"),
                error = function(e) NULL)
  if (is.null(meta_analysis_model_for_gene)) next
  meta_analysis_results_A$pooled_logFC[gene_index] <- as.numeric(meta_analysis_model_for_gene$b[1])
  meta_analysis_results_A$pooled_pval[gene_index]  <- meta_analysis_model_for_gene$pval
  meta_analysis_results_A$I2[gene_index]           <- meta_analysis_model_for_gene$I2
  meta_analysis_results_A$tau2[gene_index]         <- meta_analysis_model_for_gene$tau2
  meta_analysis_results_A$QEp[gene_index]          <- meta_analysis_model_for_gene$QEp
}

# BH correction again here (same reasoning as the Variable Definitions
# section above): ~9,000 genes tested, so the pooled p-values need their
# own multiple-testing correction just like the per-timepoint ones did.
meta_analysis_results_A$adj.P.Val <- p.adjust(meta_analysis_results_A$pooled_pval, method = "BH")
meta_analysis_results_A$Low_heterogeneity <- !is.na(meta_analysis_results_A$QEp) & meta_analysis_results_A$QEp > 0.05

message("[9/13] Meta-analysis complete: ", sum(!is.na(meta_analysis_results_A$pooled_pval)), " genes pooled across 5 timepoints.")


# =====================================================================
# STEP 10: COMPARE FINAL A/B/C, EXPORT
# =====================================================================
# Comparing three genuinely different statistical approaches to the same underlying question is deliberate.

output_directory_differential_expression_genes <- file.path(OUTPUT_DIRECTORY, "DE_Genes")
dir.create(output_directory_differential_expression_genes, showWarnings = FALSE, recursive = TRUE)

final_significant_genes_A <- meta_analysis_results_A$Gene[!is.na(meta_analysis_results_A$adj.P.Val) &
                        meta_analysis_results_A$adj.P.Val < 0.05 &
                        abs(meta_analysis_results_A$pooled_logFC) > LOG_FOLD_CHANGE_UP_THRESHOLD]
final_significant_genes_B <- unique(c(unlist(significant_upregulated_genes_B), unlist(significant_downregulated_genes_B)))
final_significant_genes_C <- unique(c(unlist(significant_upregulated_genes_C), unlist(significant_downregulated_genes_C)))

genes_significant_in_all_three_methods   <- Reduce(intersect, list(final_significant_genes_A, final_significant_genes_B, final_significant_genes_C))
genes_significant_in_both_B_and_C     <- intersect(final_significant_genes_B, final_significant_genes_C)
union_of_all_significant_genes <- Reduce(union, list(final_significant_genes_A, final_significant_genes_B, final_significant_genes_C))

message("[10/13] Final A (meta-analysis): ", length(final_significant_genes_A),
        " | Final B (dupCor pooled): ", length(final_significant_genes_B),
        " | Final C (dream pooled): ", length(final_significant_genes_C),
        " | All three agree: ", length(genes_significant_in_all_three_methods),
        " | B and C agree: ", length(genes_significant_in_both_B_and_C))

gene_method_membership_table <- data.frame(
  Gene = union_of_all_significant_genes,
  In_A = union_of_all_significant_genes %in% final_significant_genes_A,
  In_B = union_of_all_significant_genes %in% final_significant_genes_B,
  In_C = union_of_all_significant_genes %in% final_significant_genes_C
)
gene_method_membership_table$N_methods <- rowSums(gene_method_membership_table[, c("In_A","In_B","In_C")])
gene_method_membership_table <- merge(gene_method_membership_table, meta_analysis_results_A[, c("Gene","pooled_logFC","adj.P.Val","I2","QEp","Low_heterogeneity")],
                     by = "Gene", all.x = TRUE)
# Naming convention throughout this section: Method_<A/B/C/All/B_C>_<Union/Intersection/Final>_<optional: Metrics>
# All_Genes_Metrics deliberately does NOT follow that pattern -- it isn't a
# combination of anything, it's the complete unfiltered reference (every
# gene tested, ~9,055 rows), so forcing it into the same template would
# misleadingly imply it's another filtered result.
data.table::fwrite(gene_method_membership_table, file.path(output_directory_differential_expression_genes, "Method_All_Union_Metrics.tsv"), sep = "\t")
data.table::fwrite(meta_analysis_results_A, file.path(output_directory_differential_expression_genes, "All_Genes_Metrics.tsv"), sep = "\t")

# Export gene lists + ML feature matrices. Deliberately NO variance-based filtering here: total sample variance doesn't indicate Condition-
# relevance. A gene can be highly variable due to timepoint, donor baseline differences, or plain technical noise, none of which reflect
# whether it actually distinguishes SLE from Healthy. Layering a raw-variance filter on top of genes that already passed a proper
# Condition-specific statistical test doesn't add information; it risks discarding a real marker that happens to have a small, tight, low-noise
# effect (precisely the profile a good diagnostic marker should have) while keeping a noisy gene whose high variance has nothing to do with disease
# status. The statistical test itself is the appropriate filter.

expression_matrix_for_machine_learning <- as.matrix(expression)

candidate_gene_lists <- list(
  # Method_B_Final, Method_C_Final, and Method_B_C_Intersection are all
  # currently identical (44 genes each) -- B and C converge completely as
  # of this run. Kept as 3 separate files anyway, by deliberate choice, so
  # each one is traceable by name to its own step in the pipeline diagram
  # rather than needing to remember "B and C happen to match right now."
  # Same story for Method_A_Final and Method_All_Intersection (both 42).
  Method_A_Union          = unique(c(unlist(significant_upregulated_genes_by_timepoint_A),
                                       unlist(significant_downregulated_genes_by_timepoint_A))),
  Method_A_Final          = final_significant_genes_A,
  Method_B_Final          = final_significant_genes_B,
  Method_C_Final          = final_significant_genes_C,
  Method_B_C_Intersection = genes_significant_in_both_B_and_C,
  Method_All_Intersection = genes_significant_in_all_three_methods,
  # Same 44 genes as Method_All_Union_Metrics.tsv, as a plain gene list
  # (no membership flags or metrics columns) -- the union-with-no-filtering
  # counterpart to Method_All_Intersection, same pairing pattern as
  # Method_A_Union/Method_A_Union_Metrics above.
  Method_All_Union        = union_of_all_significant_genes
)

for (candidate_list_name in names(candidate_gene_lists)) {
  genes <- candidate_gene_lists[[candidate_list_name]]
  genes <- genes[genes %in% rownames(expression)]
  if (length(genes) == 0) { message("       Skipping empty list: ", candidate_list_name); next }

  write.table(data.frame(Gene = genes), file.path(output_directory_differential_expression_genes, paste0(candidate_list_name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  machine_learning_feature_matrix_log2 <- t(expression_matrix_for_machine_learning[genes, , drop = FALSE])
  write.table(machine_learning_feature_matrix_log2, file.path(output_directory_differential_expression_genes, paste0("ML_FeatureMatrix_", candidate_list_name, "_log2.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
}

message("       Gene lists and ML feature matrices written to: ", output_directory_differential_expression_genes)

# Method_A_Union with metrics attached -- same 65 genes as Method_A_Union.tsv,
# but with Method A's meta-analysis statistics merged in (pooled_logFC,
# raw and adjusted p-value, I2, tau2, QEp, Low_heterogeneity) for direct
# per-gene interpretation without needing to cross-reference two files.
method_a_union_metrics <- meta_analysis_results_A[meta_analysis_results_A$Gene %in% candidate_gene_lists[["Method_A_Union"]], ]
data.table::fwrite(method_a_union_metrics, file.path(output_directory_differential_expression_genes, "Method_A_Union_Metrics.tsv"), sep = "\t")

# ---- Python ML pipeline compatibility export ----------------------------
# sle_pipeline.py's load_dataset() expects ONE file containing both the
# gene expression columns AND 6 specific metadata columns together (it
# checks for their exact presence and errors if any are missing). The plain
# ML_FeatureMatrix_*_log2.tsv files above are expression-only, with no
# metadata -- they were never compatible with that loader as-is. Built here,
# not reconstructed externally, specifically because this script already
# has verified, correct Sample/Donor_id/Condition/time_point data in memory
# -- guessing at that mapping from outside this script risks getting donor
# assignments wrong. Uses Method_All_Intersection (the strictest, most-validated
# panel) by default; change ML_PANEL_FOR_PYTHON below if you want a
# different candidate list fed to the Python pipeline instead.
ML_PANEL_FOR_PYTHON <- "Method_All_Intersection"

if (ML_PANEL_FOR_PYTHON %in% names(candidate_gene_lists)) {
  genes_for_python <- candidate_gene_lists[[ML_PANEL_FOR_PYTHON]]
  genes_for_python <- genes_for_python[genes_for_python %in% rownames(expression)]

  python_ready_matrix <- data.frame(
    SampleID         = metadata[[SAMPLE]],
    DonorID          = as.character(metadata[[SUBJECT]]),
    Condition        = ifelse(metadata[[CONDITION]] == "SLE", 1L, 0L),   # sle_pipeline.py requires exactly 0/1 integers
    Condition_label  = as.character(metadata[[CONDITION]]),
    Timepoint        = as.integer(factor(metadata[[TIME]], levels = TIME_POINT_LEVELS)),  # ordinal position, 1-5
    Time_label       = as.character(metadata[[TIME]]),
    t(expression_matrix_for_machine_learning[genes_for_python, , drop = FALSE])
  )

  write.table(python_ready_matrix,
              file.path(output_directory_differential_expression_genes,
                        paste0("ML_FeatureMatrix_", ML_PANEL_FOR_PYTHON, "_with_Metadata.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("       Python-pipeline-ready file (with metadata) written for panel: ", ML_PANEL_FOR_PYTHON)
}

# The "final" gene set used by Step 11's summary heatmap below. Genuinely a judgment call, not a statistically-derived answer

final_panel_genes <- genes_significant_in_all_three_methods
if (length(final_panel_genes) < 2) {
  message("       Method_All_Intersection has <2 genes -- falling back to Method_B_C_Intersection for the Step 11 heatmap.")
  final_panel_genes <- genes_significant_in_both_B_and_C
}


# =====================================================================
# STEP 11: FINAL-PANEL HEATMAP (Condition only, no timepoint split)
# =====================================================================
# Deliberately different from Step 7's heatmap in two ways: 
# (1) uses final_panel_genes
# (2) columns are split by Condition ONLY, with no timepoint sub-splitting

final_panel_heatmap_output_file <- file.path(OUTPUT_DIRECTORY, "Heatmap_FinalPanel_ConditionOnly.png")

if (length(final_panel_genes) < 2) {
  message("[11/13] Not enough genes in final_panel_genes for a heatmap -- skipping.")
} else {

  raw_expression_matrix_final_panel  <- as.matrix(expression[final_panel_genes, , drop = FALSE])
  z_scored_expression_matrix_final_panel <- t(scale(t(raw_expression_matrix_final_panel)))
  z_scored_expression_matrix_final_panel[is.na(z_scored_expression_matrix_final_panel)] <- 0

  metadata_ordered_final_panel <- metadata[match(colnames(z_scored_expression_matrix_final_panel), metadata[[SAMPLE]]), , drop = FALSE]
  metadata_ordered_final_panel[[CONDITION]] <- factor(metadata_ordered_final_panel[[CONDITION]], levels = c("SLE","Healthy"))
  column_order_index_final_panel <- order(metadata_ordered_final_panel[[CONDITION]])
  z_scored_expression_matrix_final_panel <- z_scored_expression_matrix_final_panel[, column_order_index_final_panel, drop = FALSE]
  metadata_ordered_final_panel <- metadata_ordered_final_panel[column_order_index_final_panel, , drop = FALSE]

  # Condition-only split -- no time_point anywhere in this annotation or split.
  column_split_groups_final_panel <- droplevels(metadata_ordered_final_panel[[CONDITION]])

  annotation_dataframe_final_panel <- data.frame(Condition = droplevels(metadata_ordered_final_panel[[CONDITION]]))
  rownames(annotation_dataframe_final_panel) <- metadata_ordered_final_panel[[SAMPLE]]

  top_annotation_final_panel <- ComplexHeatmap::HeatmapAnnotation(
    df  = annotation_dataframe_final_panel,
    col = list(Condition = c(SLE = "#d73027", Healthy = "#1a9850")),
    annotation_name_side = "left",
    simple_anno_size = grid::unit(3, "mm")
  )

  color_mapping_function_final_panel <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  png(final_panel_heatmap_output_file, width = 1400, height = 1000, res = 170)
  heatmap_object_final_panel <- ComplexHeatmap::Heatmap(
    z_scored_expression_matrix_final_panel,
    name = "Expression (z)",
    col  = color_mapping_function_final_panel,
    top_annotation = top_annotation_final_panel,
    column_split = column_split_groups_final_panel,
    gap = grid::unit(1.5, "mm"),
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = FALSE,
    cluster_columns = FALSE,
    clustering_method_rows = "ward.D2",
    column_title = "Columns = SAMPLES (SLE vs Healthy only -- no timepoint split)",
    row_title    = "Rows = FINAL PANEL GENES (z-scored per gene)",
    heatmap_legend_param = list(
      at = c(-2, 0, 2),
      labels = c("Down-regulated", "Mean", "Up-regulated")
    )
  )
  ComplexHeatmap::draw(heatmap_object_final_panel, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE)
  dev.off()

  message("[11/13] Final-panel heatmap saved: ", final_panel_heatmap_output_file,
          "  (genes: ", nrow(z_scored_expression_matrix_final_panel), ", samples: ", ncol(z_scored_expression_matrix_final_panel), ")")
}


# =====================================================================
# STEP 12: GLOBAL TIME x CONDITION F-TEST (TiSA-style)
# =====================================================================
# The purpose is to formally test whether the SLE-vs-Healthy effect ITSELF changes shape across gestation (a true statistical interaction), which is a
# different and more rigorous question than anything Steps 3-9 answer Method A shows the effect at each timepoint separately without directly
# testing whether those 5 estimates differ from each other in a statistically defensible way, and the meta-analysis's heterogeneity
# stats (I2/QEp) are per-gene, not one combined test across all interaction terms simultaneously. This step follows the "temporal DEG analysis"
# axis from TiSA (Lefol et al. 2023), which explicitly evaluates within-group temporal change alongside the conditional comparison.

condition_time_dataframe <- data.frame(
  Condition = droplevels(metadata[[CONDITION]]),
  Time      = droplevels(metadata[[TIME]])
)

message("[12/13] Sample counts per Condition x Time:")
print(addmargins(table(condition_time_dataframe$Condition, condition_time_dataframe$Time)))

design_matrix_interaction_test <- model.matrix(~ 0 + Condition + Time + Condition:Time, data = condition_time_dataframe)
colnames(design_matrix_interaction_test) <- make.names(colnames(design_matrix_interaction_test))

message("        Design columns (first 10): ", paste(head(colnames(design_matrix_interaction_test), 10), collapse = ", "))
term_assignment_vector  <- attr(design_matrix_interaction_test, "assign")
term_labels <- attr(attr(design_matrix_interaction_test, "terms"), "term.labels")
message("        Term labels: ", paste(term_labels, collapse = " | "))
message("        Assign counts by term index: ", paste(table(term_assignment_vector), collapse = " | "))

# Reuses the correlation estimate from Step 2's base model rather than
# recomputing duplicateCorrelation a third time -- the repeated-measures
# structure (which donors contributed which samples) hasn't changed, only
# the fixed-effect formula has, so the same correlation estimate applies.
base_model_fit_interaction_test <- lmFit(
  expression, design_matrix_interaction_test,
  block       = metadata[[SUBJECT]],
  correlation = donor_correlation_estimate_A$consensus.correlation
) |> eBayes()

interaction_term_index <- which(term_labels == "Condition:Time")
interaction_column_indices  <- which(term_assignment_vector == interaction_term_index)
interaction_column_names     <- colnames(design_matrix_interaction_test)[interaction_column_indices]

if (length(interaction_column_names) == 0) {
  interaction_column_names <- grep("^Condition.*\\.Time", colnames(design_matrix_interaction_test), value = TRUE)
}
message("        Detected interaction columns: ",
        ifelse(length(interaction_column_names) == 0, "<none>", paste(interaction_column_names, collapse = ", ")))

if (length(interaction_column_names) == 0) {
  stop("No interaction columns detected. This means either (a) the model formula wasn't used, ",
       "(b) a factor has a single level, or (c) some Condition x Time cells are empty. ",
       "Check the balance table above.")
}

# Joint F-test across ALL interaction coefficients simultaneously (not one
# t-test per coefficient) -- this is what actually tests "does the
# Condition effect vary across time" as a single combined question, rather
# than four separate, individually-uncorrected sub-questions.
interaction_contrast_matrix <- matrix(0, nrow = ncol(design_matrix_interaction_test), ncol = length(interaction_column_names),
            dimnames = list(colnames(design_matrix_interaction_test), interaction_column_names))
for (interaction_column_index in seq_along(interaction_column_names)) {
  interaction_contrast_matrix[interaction_column_names[interaction_column_index], interaction_column_index] <- 1
}

fitted_interaction_model <- contrasts.fit(base_model_fit_interaction_test, interaction_contrast_matrix)
fitted_interaction_model <- eBayes(fitted_interaction_model)

interaction_top_table <- topTable(fitted_interaction_model, number = Inf, sort.by = "F")
# BUG FIX: every other topTable() output in this script explicitly adds
# gene names as a column before saving (rownames aren't written to a TSV
# by fwrite() unless you do this) -- this one was missing that step, so
# the saved file had no way to identify which row belonged to which gene.
interaction_top_table$Gene <- rownames(interaction_top_table)
data.table::fwrite(interaction_top_table,
                   file.path(OUTPUT_DIRECTORY, "Global_Time_by_Condition_Ftest.tsv"),
                   sep = "\t")

message("        Global Time x Condition F-test saved.")


# =====================================================================
# STEP 13: FOREST PLOT - METHOD B vs METHOD C, TOP GENES
# =====================================================================
# Directly visualizes the Method B vs Method C comparison, by plotting their pooled logFC with 95% confidence intervals side by side for the top genes -- agreement (bars
# overlapping, similar point estimates. Disagreement is exactly where gene-specific correlation estimationwould show up visually.

if (!is.null(full_results_B) && !is.null(full_results_C)) {

  top_genes_for_forest_plot <- union(
    head(full_results_B$Gene[order(full_results_B$adj.P.Val)], 15),
    head(full_results_C$Gene[order(full_results_C$adj.P.Val)], 15)
  )

  forest_plot_data_B <- full_results_B[full_results_B$Gene %in% top_genes_for_forest_plot,
                                        c("Gene","logFC","CI.L","CI.R")]
  forest_plot_data_B$Method <- "B (duplicateCorrelation)"

  forest_plot_data_C <- full_results_C[full_results_C$Gene %in% top_genes_for_forest_plot,
                                        c("Gene","logFC","CI.L","CI.R")]
  forest_plot_data_C$Method <- "C (dream)"

  forest_plot_data <- rbind(forest_plot_data_B, forest_plot_data_C)

  # Order genes by Method B's effect size, purely so the plot reads
  # top-to-bottom in a sensible order rather than alphabetically.
  gene_order_for_plot <- forest_plot_data_B$Gene[order(forest_plot_data_B$logFC)]
  forest_plot_data$Gene <- factor(forest_plot_data$Gene, levels = gene_order_for_plot)

  forest_plot_object <- ggplot(forest_plot_data, aes(x = logFC, y = Gene, color = Method)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = CI.L, xmax = CI.R), height = 0.3,
                    position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    scale_color_manual(values = c("B (duplicateCorrelation)" = "#1b9e77", "C (dream)" = "#d95f02")) +
    labs(
      title = "Method B vs Method C -- pooled SLE-vs-Healthy effect, top genes",
      x = "logFC (SLE - Healthy) with 95% CI",
      y = NULL
    ) +
    theme_minimal(base_size = 11)

  forest_plot_output_file <- file.path(OUTPUT_DIRECTORY, "ForestPlot_MethodB_vs_MethodC.png")
  ggsave(forest_plot_output_file, forest_plot_object, width = 8, height = max(4, length(top_genes_for_forest_plot) * 0.3), dpi = 170)

  message("[13/13] Forest plot saved: ", forest_plot_output_file,
          "  (genes shown: ", length(top_genes_for_forest_plot), ")")

} else {
  message("[13/13] Skipping forest plot -- Method B and/or Method C did not run (check RUN_B_POOLED_MODEL / RUN_C_DREAM_MODEL).")
}

message("Done.")

}  # end of run_deg_pipeline()


# ==========================================================================
# ACTUALLY RUN IT
# ==========================================================================
# Sourcing this file up to this point only DEFINED run_deg_pipeline() -- it
# hasn't run yet. This line is what actually starts it. Comment this one
# line out (put a # in front of it) if you ever want to source the file
# just to have the function available without immediately running all 13
# steps.
run_deg_pipeline()
