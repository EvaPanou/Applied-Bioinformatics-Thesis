# ============================================================================
# 02_ALASCA_Analysis.R
#
# Repeated-measures (donor x timepoint x Condition) analysis of the cleaned,
# filtered expression matrix produced by 01_RawData_Analysis.ipynb, using the
# RM-ASCA+ method via the ALASCA package.
#
# WHAT CHANGED vs. the previous version of this script:
# The previous run used separate_effects = TRUE with formula
# value ~ time_point * Condition + (1|Donor_id), which -- per ALASCA's own
# documentation (Jarmund, Madssen & Giskeodegard 2022, Frontiers Mol.
# Biosci.; https://andjar.github.io/ALASCA/reference/ALASCA.html) -- is
# DOCUMENTED, INTENDED behavior for a two-way interaction model: it always
# splits into exactly TWO matrices, the first main effect alone, and
# "second main effect + interaction" bundled together. That's not a mistake
# -- it's what separate_effects=TRUE does by design for time*Condition. But
# it means "Condition" and "time_point:Condition" were never separable in
# that run.
#
# This version instead uses ALASCA's `effects` argument (documented at the
# link above: "The effect matrices can be specified with effects, e.g.,
# c('time','time+group+time:group','group+time:group')") to explicitly
# request THREE separate matrices: time_point alone, Condition alone, and
# the interaction alone. Treating Condition as cleanly separable from the
# interaction is reasonable here specifically because 02_DGE_Analysis's own
# Global_Time_by_Condition_Ftest.tsv already found 0 of 9,055 genes with a
# significant Condition x Time interaction -- there's essentially no real
# interaction being papered over by asking for a "pure" Condition effect.
#
# CAVEAT, stated plainly: the exact syntax for the `effects` argument (term
# strings, separator characters) is confirmed to exist from the package
# docs, but the precise behavior with participant_column/random effects
# combined with a fully-separated 3-term `effects` vector has NOT been
# tested against a live install in this session (no R available here).
# Run Step 4 below, then immediately run the `str()`/`print()` calls before
# trusting effect indices 1/2/3 downstream -- rename/reindex if what comes
# back doesn't match the ordering assumed here.

# ============================================================================
# Step 0: Which donor cohort this script is reading
# ============================================================================
# Set this by hand to match the run_label used when 01_RawData_Analysis.ipynb
# was last run and saved its outputs.

run_label <- "non_NP"   # set this to match whichever notebook run you want to read

# ============================================================================
# Step 1: Installing/updating packages (only needs to be done once, not every run)
# ============================================================================
# Run this block once, manually, the first time -- or any time
# packageVersion("ggplot2") comes back below 3.5.0, since create.dir= in
# the ggsave() calls below (Step 5/6) needs 3.5.0 or newer:
#
# install.packages("devtools")
# devtools::install_github("andjar/ALASCA")
# install.packages("ggplot2")

# import libraries as usual for every run
library(ALASCA)
library(dplyr)
library(ggplot2)


# ============================================================================
# Step 2: Loading the files the Python notebook already cleaned and filtered
# ============================================================================
# Loading the ouput files of the process that already happened in 
# 01_RawData_Analysis.ipynb. We only read the two files specific for ALASCA.
# adjust path if this script sits elsewhere

input_folder <- file.path("C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/01_RawData/02_Analysis_Results")  
output_folder <- file.path("C:/Users/evapa/Desktop/AppBio/Thesis/Git_Repo/Experiment/01_RawData/04_ALASCA_Output")
# Every ggsave() call below passes create.dir = TRUE, which -- since
# ggplot2 3.5.0 -- tells ggsave to just create a missing directory itself,
# with no interactive prompt at all under any circumstance. That's what
# ate the rest of the script as literal text last run.

expression_path <- file.path(input_folder, "15_final_filtered_expression_non_NP_donors.tsv")
metadata_path <- file.path(input_folder, "19_metadata_for_alasca_non_NP_donors.tsv")

final_filtered_expression_df <- read.delim(expression_path, sep = "\t", row.names = 1, check.names = FALSE)
dim(final_filtered_expression_df) # check shape
metadata_for_alasca <- read.delim(metadata_path, sep = "\t", row.names = 1, check.names = FALSE)
dim(metadata_for_alasca) # check shape

# quick sanity check: every sample in the expression matrix must have metadata,
# and vice versa, before we go any further
samples_in_expression <- colnames(final_filtered_expression_df)
samples_in_metadata <- rownames(metadata_for_alasca)
cat("Samples in expression matrix but missing from metadata:",
    setdiff(samples_in_expression, samples_in_metadata), "\n")
cat("Samples in metadata but missing from expression matrix:",
    setdiff(samples_in_metadata, samples_in_expression), "\n")

# ============================================================================
# Step 3: Reshaping into the "long format" ALASCA expects
# ============================================================================
# ALASCA needs one row per (sample, gene) combination, with columns identifying
# the donor, the timepoint, the group (Condition), and the measured value.

expression_long_format <- as.data.frame(t(final_filtered_expression_df))
expression_long_format$Sample <- rownames(expression_long_format)

expression_long_format <- merge(expression_long_format,
                                 data.frame(Sample = rownames(metadata_for_alasca), metadata_for_alasca,
                                            row.names = NULL, check.names = FALSE),
                                 by = "Sample")

gene_columns <- setdiff(colnames(expression_long_format), c("Sample", "Donor_id", "time_point", "Condition"))

expression_long_format_reshaped <- reshape(expression_long_format,
                                           varying = gene_columns,
                                           v.names = "value",
                                           timevar = "variable",
                                           times = gene_columns,
                                           idvar = "Sample",
                                           direction = "long")

cat("Long-format data has", nrow(expression_long_format_reshaped), "rows (samples x genes).\n")


# ============================================================================
# Step 4: Fit the RM-ASCA+ model with THREE explicitly-separated effects
# ============================================================================
# Same formula as before (time_point * Condition + (1|Donor_id)), but now
# passing `effects` explicitly instead of relying on separate_effects=TRUE's
# automatic 2-way guess. This asks ALASCA for three independent matrices:
#   1. time_point            -- the pure gestational-time effect, shared
#                                across both groups
#   2. Condition              -- the pure SLE-vs-Healthy main effect, with
#                                no time-interaction folded in
#   3. time_point:Condition   -- the interaction alone: does the Condition
#                                effect change shape across gestation? (this
#                                is the multivariate analogue of Step 12's
#                                Global_Time_by_Condition_Ftest -- expect it
#                                to look close to flat/near-zero separation,
#                                consistent with that test finding 0/9,055
#                                genes significant)
#
# NOTE: the exact expected column names and `effects` string syntax should
# still be double-checked against the current ALASCA package documentation
# / str(alasca_model) output before trusting this blindly -- see the caveat
# in the header comment above.

alasca_model <- ALASCA(df = expression_long_format_reshaped,
                        formula = value ~ time_point * Condition + (1 | Donor_id),
                        effects = c("time_point", "Condition", "time_point:Condition"),
                        participant_column = "Donor_id")
                      # validate = TRUE,          # note to self: uncomment if you feel validation is appropriate
                      # n_validation_runs = 20)   # bootstrap/permutation resampling, needed for Step 5 below

print(alasca_model)
str(alasca_model, max.level = 1) # CHECK HERE which effect number ended up as time_point / Condition / interaction --
                                  # the header caveat means this ordering isn't guaranteed to match the comments below.


# ============================================================================
# Step 5: Plot all three separated effects
# ============================================================================
# Loops over the 3 effects instead of hardcoding two separate blocks like the
# previous version did -- add a 4th entry here later if you ever re-run with
# a different formula (e.g. adding a covariate) that produces more effects.

effect_definitions <- list(
  list(index = 1, slug = "time_point",         label = "Time (pure gestational-time effect)"),
  list(index = 2, slug = "condition",           label = "Condition (pure SLE vs Healthy main effect)"),
  list(index = 3, slug = "time_by_condition",   label = "Time x Condition interaction")
)

for (effect_definition in effect_definitions) {
  effect_plot <- plot(alasca_model, effect = effect_definition$index, type = "effect")
  outfile <- file.path(output_folder, paste0("01_", effect_definition$slug, "_effect_plot.png"))
  ggsave(outfile, plot = effect_plot, width = 8, height = 5, units = "in", dpi = 300, create.dir = TRUE)
  print(effect_plot)
  message("Saved: ", outfile, "  (", effect_definition$label, ")")
}


# ============================================================================
# Step 6: PC1/PC2 group scores -- for the Condition effect AND the
# interaction effect
# ============================================================================
# ASSUMPTION (still unverified, same as before): get_scores() returns a
# data.frame/data.table with a "score" column plus time_point/Condition
# columns. Run str(scores_pc1)/head(scores_pc1) first and adjust the
# column names below if the real output differs.
#
# Two score plots now, not one:
#  - Condition effect (2): the cleaned-up equivalent of the old 6-dot plot,
#    minus the interaction folded in -- this SHOULD look similar to before
#    (Healthy vs SLE clearly apart), since the interaction being removed is
#    apparently negligible.
#  - Interaction effect (3): NEW. If the Global_Time_by_Condition_Ftest
#    result holds at the multivariate level too, the 5 SLE-timepoint dots
#    here should cluster together with little spread, near the Healthy
#    dot -- i.e. NOT separated by Condition -- since a real interaction is
#    what would make them differ from each other and from Healthy.

plot_group_scores <- function(alasca_model, effect_index, plot_title, output_file) {
  # Only used for the interaction effect below -- the only one of the three
  # whose score table actually has a time_point column to facet by.
  scores_pc1 <- get_scores(alasca_model, effect = effect_index, component = 1)
  scores_pc2 <- get_scores(alasca_model, effect = effect_index, component = 2)

  str(scores_pc1)
  head(scores_pc1)

  scores_pc1 <- scores_pc1 %>% rename(PC1 = score) %>% select(time_point, Condition, PC1)
  scores_pc2 <- scores_pc2 %>% rename(PC2 = score) %>% select(time_point, Condition, PC2)

  scores_combined <- merge(scores_pc1, scores_pc2, by = c("time_point", "Condition"))

  group_scores_plot <- ggplot(scores_combined, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Healthy" = "#0072B2", "SLE" = "#E69F00")) +
    labs(title = plot_title, x = "PC1", y = "PC2") +
    theme_minimal()

  ggsave(output_file, plot = group_scores_plot, width = 8, height = 7, units = "in", dpi = 300, create.dir = TRUE)
  print(group_scores_plot)
  message("Saved: ", output_file)

  invisible(scores_combined)
}

# ---- Condition-only effect (2 points: Healthy, SLE -- no time dimension,
# so no time_point column to select/merge on. Kept as its own simple block
# rather than folded into plot_group_scores(), since that function assumes
# a time_point column the Condition effect's score table doesn't have. ----
condition_scores_pc1 <- get_scores(alasca_model, effect = 2, component = 1) %>%
  rename(PC1 = score) %>% select(Condition, PC1)
condition_scores_pc2 <- get_scores(alasca_model, effect = 2, component = 2) %>%
  rename(PC2 = score) %>% select(Condition, PC2)
condition_scores <- merge(condition_scores_pc1, condition_scores_pc2, by = "Condition")

condition_plot <- ggplot(condition_scores, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Healthy" = "#0072B2", "SLE" = "#E69F00")) +
  labs(title = "Group scores -- Condition main effect (interaction removed)", x = "PC1", y = "PC2") +
  theme_minimal()

condition_outfile <- file.path(output_folder, "02_condition_only_pca_by_condition.png")
ggsave(condition_outfile, plot = condition_plot, width = 8, height = 7, units = "in", dpi = 300, create.dir = TRUE)
print(condition_plot)
message("Saved: ", condition_outfile)

interaction_scores <- plot_group_scores(
  alasca_model, effect_index = 3,
  plot_title  = "Group scores -- Time x Condition interaction",
  output_file = file.path(output_folder, "03_interaction_pca_by_condition.png")
)


# ============================================================================
# Step 7: Save the ALASCA outputs
# ============================================================================
# Saved under a distinct filename (not overwriting the earlier 2-effect run)
# so you can compare both versions if useful.

saveRDS(alasca_model, file = file.path(output_folder, "03_alasca_model_no_NP_3effects.rds"))

print(paste("ALASCA model object saved to:",
            file.path(output_folder, "03_alasca_model_no_NP_3effects.rds")))
