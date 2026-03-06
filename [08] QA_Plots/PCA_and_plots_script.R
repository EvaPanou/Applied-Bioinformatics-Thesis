# ===============================
# Visualization & PCA Utilities
# (All A/B FULL contrasts + 3 PCA views)
# ===============================

# -----------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ---------------------------
# 0) CONFIG (edit as needed)
# ---------------------------
EXPR_FILE <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/POST_DE/Input/08_df_no_low_genes.tsv"
META_FILE <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/POST_DE/Input/GSE108497_updated_metadata_II.csv"
OUTDIR    <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/POST_DE"

# DEG lists for PCA subsets
DEG_LIST_LIMMA_TSV <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/POST_DE/Input/AllDEGs.tsv"
DEG_LIST_OLD_TSV   <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/POST_DE/Input/13_Common_Top500_3_Methods.tsv"

# Create subdirectories for organized output
MA_DIR  <- file.path(OUTDIR, "QA_Plots", "MA_Plots")
V_DIR   <- file.path(OUTDIR, "QA_Plots", "Volcano_Plots")
PCA_DIR <- file.path(OUTDIR, "QA_Plots", "PCA_Plots")
dir.create(MA_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(V_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(PCA_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1) Load data helpers
# ---------------------------
load_expression <- function(path) {
  x <- data.table::fread(path) |> as.data.frame()
  rownames(x) <- x[[1]]
  x[[1]] <- NULL
  x
}

load_metadata <- function(path) {
  m <- data.table::fread(path) |> as.data.frame()
  m$Sample <- trimws(m$Sample)
  m
}

sanitize_label <- function(path) {
  base <- basename(path)
  base <- sub("\\.tsv$", "", base)
  gsub("[<>:\"/\\\\|\\?\\*]", "_", base)
}

# ---------------------------
# 2) Volcano & MA plots
# ---------------------------
plot_volcano <- function(tt, outfile,
                         fdr_cut = 0.05, lfc_cut = 1, top_n_labels = 15) {
  if (!all(c("logFC", "adj.P.Val") %in% names(tt))) {
    stop("Table must contain 'logFC' and 'adj.P.Val' columns for volcano plot.")
  }
  if (!"Gene" %in% names(tt)) tt$Gene <- rownames(tt)
  
  tt <- tt |>
    mutate(sig = adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
  tt$neglog10FDR <- -log10(pmax(tt$adj.P.Val, 1e-300))
  
  base <- ggplot(tt, aes(x = logFC, y = neglog10FDR)) +
    geom_point(aes(alpha = sig), size = 1.2) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    labs(x = "log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot") +
    theme_minimal(base_size = 12) +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3), guide = "none")
  
  lab <- tt |> arrange(adj.P.Val) |> head(top_n_labels)
  p <- base + ggrepel::geom_text_repel(data = lab, aes(label = Gene), size = 3, max.overlaps = 50)
  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
}

plot_MA <- function(tt, outfile, lfc_cut = 1) {
  if (!all(c("AveExpr", "logFC") %in% names(tt))) {
    stop("Table must contain 'AveExpr' and 'logFC' columns for MA plot.")
  }
  
  p <- ggplot(tt, aes(x = AveExpr, y = logFC)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    labs(x = "Average Expression (log2)", y = "log2 Fold Change", title = "MA Plot") +
    theme_minimal(base_size = 12)
  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
}

# ---------------------------
# 3) PCA (generic)
# ---------------------------
run_pca <- function(mat, meta = NULL, title = "PCA", outfile) {
  mat <- as.matrix(mat)
  keep <- apply(mat, 1, function(v) sd(v, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  
  pcs <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  df <- as.data.frame(pcs$x[, 1:2, drop = FALSE])
  df$Sample <- rownames(df)
  
  if (!is.null(meta)) {
    meta$Sample <- trimws(meta$Sample)
    df <- df |> dplyr::left_join(meta, by = "Sample")
  }
  
  p <- ggplot(df, aes(x = PC1, y = PC2, color = Condition, shape = time_point)) +
    geom_point(size = 2, alpha = 0.9) +
    labs(title = title) +
    theme_minimal(base_size = 12)
  
  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
  invisible(list(pca = pcs, scores = df))
}

# ---------------------------
# 4) DEG helper for PCA
# ---------------------------
get_deg_genes_from_list <- function(path, expr, subset_label = "DEG subset") {
  if (!file.exists(path)) stop(paste0("File not found: ", path))
  dg <- data.table::fread(path)
  col_gene <- intersect(names(dg), c("Gene", "gene", "SYMBOL", "symbol"))
  if (length(col_gene) == 0) stop(paste0("Cannot find gene column in ", subset_label))
  genes <- unique(dg[[col_gene[1]]])
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < 2) stop(paste0("Too few genes for PCA in ", subset_label))
  genes
}

# ---------------------------
# 5) MAIN — batch volcano/MA + PCA plots
# ---------------------------
main <- function() {
  expr <- load_expression(EXPR_FILE)
  meta <- load_metadata(META_FILE)
  
  # ---- 5A) Volcano + MA for all contrasts ----
  full_files <- list.files(
    OUTDIR,
    pattern = "^[AB]_FULL_.*\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(full_files) == 0) stop("No A_FULL_*/B_FULL_* tables found under OUTDIR.")
  
  message("Found ", length(full_files), " FULL tables. Generating Volcano + MA plots...")
  
  for (ff in full_files) {
    label <- sanitize_label(ff)
    message("  - Processing: ", label)
    tt <- data.table::fread(ff)
    if (!"Gene" %in% names(tt)) tt$Gene <- rownames(tt)
    
    volcano_file <- file.path(V_DIR, paste0("Volcano_", label, ".png"))
    ma_file      <- file.path(MA_DIR, paste0("MA_", label, ".png"))
    
    plot_volcano(tt, volcano_file)
    plot_MA(tt, ma_file)
  }
  
  # ---- 5B) PCA plots ----
  message("Running PCA: All Genes...")
  run_pca(expr, meta, title = "PCA — All Genes",
          outfile = file.path(PCA_DIR, "PCA_AllGenes.png"))
  
  message("Running PCA: LIMMA DEG subset...")
  deg_limma <- get_deg_genes_from_list(DEG_LIST_LIMMA_TSV, expr, "LIMMA DEG subset")
  run_pca(expr[deg_limma, , drop = FALSE], meta,
          title = "PCA — LIMMA DEG Subset",
          outfile = file.path(PCA_DIR, "PCA_LIMMA_DEG_Subset.png"))
  
  message("Running PCA: Old DEG subset...")
  deg_old <- get_deg_genes_from_list(DEG_LIST_OLD_TSV, expr, "Old DEG subset")
  run_pca(expr[deg_old, , drop = FALSE], meta,
          title = "PCA — Old DEG Subset",
          outfile = file.path(PCA_DIR, "PCA_Old_DEG_Subset.png"))
  
  message("Done! Results saved under:\n",
          "- Volcano plots: ", V_DIR, "\n",
          "- MA plots:      ", MA_DIR, "\n",
          "- PCA plots:     ", PCA_DIR)
}

if (identical(environment(), globalenv())) main()

