# ==========================================
# LIMMA TIMESERIES DEG ANALYSIS - FUNCTIONS
# ==========================================

# ---- INSTALLATIONS ------------------
# =====================================

install.packages(c("tidyverse","data.table"), dependencies = TRUE)

# ---- Bioconductor setup ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# ---- Core analysis packages ----
BiocManager::install(c("limma","ComplexHeatmap","circlize","UpSetR","matrixStats"), 
                     update = TRUE, ask = FALSE)

# ---- Enrichment analysis packages ----
BiocManager::install(c("clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot"), 
                     update = TRUE, ask = FALSE)

# ---- Infrastructure dependencies (fixes GenomeInfoDb / AnnotationDbi errors) ----
BiocManager::install(c("GenomeInfoDb","S4Vectors","IRanges","XVector","BiocGenerics",
                       "Biobase","AnnotationDbi","DBI","RSQLite"), 
                     update = TRUE, ask = FALSE)

# ---- Database annotation packages ----
BiocManager::install(c("GO.db","reactome.db","KEGGREST"), 
                     update = TRUE, ask = FALSE)


# ---- LOAD PACKAGES --------------------
# =======================================

suppressPackageStartupMessages({
  library(limma)             # linear modeling and eBayes
  library(ComplexHeatmap)    # advanced heatmaps
  library(circlize)          # color mapping for heatmaps
  library(UpSetR)            # intersection visualization
  library(tidyverse)         # general data wrangling
  library(data.table)        # fast data import/export
  library(matrixStats)       # row/column variance calculations
  library(clusterProfiler)   # enrichment (GO, KEGG)
  library(org.Hs.eg.db)      # human gene annotation database
  library(ReactomePA)        # Reactome pathway enrichment
  library(enrichplot)        # enrichment visualization
  library(GO.db)             # GO term database
  library(KEGGREST)          # KEGG REST API access
  library(GenomeInfoDb)      # genomic data structures
  library(S4Vectors)         # Bioconductor data containers
  library(IRanges)           # genomic ranges
  library(Biobase)           # expression data containers
  library(BiocGenerics)      # core Bioconductor functions
  library(AnnotationDbi)     # annotation database interface
  library(DBI)               # database connectivity
  library(RSQLite)           # SQLite backend
})


`%+%` <- function(a,b) paste0(a,b)

# ---- FILE PATHS --------------------------
# ==========================================

setwd("C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis")
EXPR_FILE <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/08_df_no_low_genes.tsv"
META_FILE <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/GSE108497_updated_metadata_II.csv"
OUTDIR    <- "C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/DEG analysis/Results"

# ---- VARIABLE DEFINITIONS -----------------
# ===========================================

SAMPLE   <- "Sample"
CONDITION<- "Condition"   # Healthy / SLE
TIME     <- "time_point"   # PP, <16 weeks, etc
SUBJECT  <- "Donor_id"

FDRS      <- c(0.05, 0.01)
LFC_UP    <-  1
LFC_DOWN  <- -1

# temporal mode: "adjacent" or "all" (TiSA-style do_all_combinations)
TEMPORAL_MODE <- "adjacent"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- HELPER FUNCTIONS -------------------
# =========================================

load_inputs <- function(expr_file, meta_file) {
  expr <- data.table::fread(expr_file) |> as.data.frame()
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  
  meta <- data.table::fread(meta_file) |> as.data.frame()
  # trim whitespace on both
  meta[[SAMPLE]] <- trimws(meta[[SAMPLE]])
  colnames(expr)  <- trimws(colnames(expr))
  
  # Keep only common samples and align order
  common <- intersect(colnames(expr), meta[[SAMPLE]])
  expr <- expr[, common, drop = FALSE]
  meta <- meta[meta[[SAMPLE]] %in% common, , drop = FALSE]
  meta <- meta[match(colnames(expr), meta[[SAMPLE]]), , drop = FALSE]
  
  # Convert to factors
  meta[[CONDITION]] <- factor(meta[[CONDITION]], levels = c("Healthy","SLE"))
  meta[[TIME]] <- factor(
    meta[[TIME]],
    levels = c("<16 weeks", "16-23 weeks", "24-31 weeks", "32-40 weeks", "PP")
  )
  meta[[SUBJECT]] <- factor(meta[[SUBJECT]])
  
  list(expr = expr, meta = meta)
}

# ---- Build Design matrix for DEG analysis -------

build_design <- function(expr, meta) {
  design <- model.matrix(~ 0 + meta[[CONDITION]]:meta[[TIME]])
  colnames(design) <- make.names(colnames(design))
  corfit <- duplicateCorrelation(expr, design, block = meta[[SUBJECT]])
  fit0   <- lmFit(expr, design, block = meta[[SUBJECT]], correlation = corfit$consensus.correlation)
  message("Within-subject correlation: ", round(corfit$consensus.correlation,3))
  list(design=design, fit0=fit0, corfit=corfit)
}

# Build a lookup : (Condition, Timepoint) -> Design column
make_column_lookup <- function(design, meta) {
  lookup <- list()
  for (j in seq_len(ncol(design))) {
    # Find which samples belong to this column (non-zero)
    idx <- which(design[, j] != 0)
    
    # Identify their unique Condition and Timepoint
    cond <- unique(as.character(meta[[CONDITION]][idx]))
    tp   <- unique(as.character(meta[[TIME]][idx]))
    
    # Sanity check: exactly one condition and one timepoint
    if (length(cond) != 1 || length(tp) != 1) {
      stop("Design column ", colnames(design)[j],
           " does not map uniquely to a single (Condition, Timepoint).")
    }
    
    # Create a key like "SLE||<16 weeks>"
    key <- paste(cond, tp, sep = "||")
    lookup[[key]] <- colnames(design)[j]
  }
  lookup
}

# ---- Define Lookup Function ------

# Lookup function: given Condition + Timepoint -> return column name
column_id_lookup <- function(lookup, condition_label, time_label) {
  key <- paste(as.character(condition_label), as.character(time_label), sep = "||")
  if (!key %in% names(lookup)) {
    stop("No design column for ", condition_label, " @ ", time_label,
         ". Check Condition/Timepoint levels.")
  }
  lookup[[key]]
}

log_set_sizes <- function(sets, label) {
  sizes <- lengths(sets)
  nonzero <- sizes[sizes > 0]
  if (length(nonzero) == 0) {
    message("[", label, "] All sets are empty at current thresholds.")
  } else {
    msg <- paste(sprintf("%s=%d", names(nonzero), nonzero), collapse = " | ")
    message("[", label, "] Non-empty sets: ", msg)
  }
}

# ---- Run Conditional/Pairwise DEG analysis ------

run_conditional_DE <- function(fit0, design, meta, lookup) {
  tlevels <- levels(meta[[TIME]])
  cons <- list()
  
  for (tp in tlevels) {
    sle_col <- column_id_lookup(lookup, "SLE", tp)
    hl_col  <- column_id_lookup(lookup, "Healthy", tp)
    cons[[paste0("SLE_vs_Healthy_at_", tp)]] <- paste0(sle_col, " - ", hl_col)
  }
  
  C <- makeContrasts(contrasts = unlist(cons), levels = design)
  fit <- contrasts.fit(fit0, C) |> eBayes()
  fit
}

# ---- Run Temporal/Longitudinal DEG analysis ------

run_temporal_DE <- function(fit0, design, meta, lookup, mode = c("adjacent","all")) {
  mode <- match.arg(mode)
  tlevels <- levels(meta[[TIME]])
  cons <- list()
  
  for (cond in levels(meta[[CONDITION]])) {
    if (mode == "adjacent") {
      for (i in seq_len(length(tlevels) - 1)) {
        tp1 <- tlevels[i]; tp2 <- tlevels[i + 1]
        c1 <- column_id_lookup(lookup, cond, tp1)  # earlier (control)
        c2 <- column_id_lookup(lookup, cond, tp2)  # later   (experiment)
        cons[[paste0(cond, "_", tp2, "_vs_", tp1)]] <- paste0(c2, " - ", c1)
      }
    } else {
      for (i in 1:(length(tlevels) - 1)) {
        for (j in (i + 1):length(tlevels)) {
          tp1 <- tlevels[i]; tp2 <- tlevels[j]
          c1 <- column_id_lookup(lookup, cond, tp1)
          c2 <- column_id_lookup(lookup, cond, tp2)
          cons[[paste0(cond, "_", tp2, "_vs_", tp1)]] <- paste0(c2, " - ", c1)
        }
      }
    }
  }
  
  C <- makeContrasts(contrasts = unlist(cons), levels = design)
  fit <- contrasts.fit(fit0, C) |> eBayes()
  fit
}

# ---- Save outputs in tsv format -----

save_DE_tables <- function(fit, stem, outdir, fdrs=FDRS, lfc_up=LFC_UP, lfc_down=LFC_DOWN) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (coef_i in seq_len(ncol(fit))) {
    ct <- colnames(fit)[coef_i]
    tt <- topTable(fit, coef=coef_i, number=Inf, sort.by="P")
    tt$Gene <- rownames(tt)
    fwrite(tt, file.path(outdir, paste0(stem, "_FULL_", make.names(ct), ".tsv")), sep="\t")
    for (f in fdrs) {
      up   <- subset(tt, adj.P.Val < f & logFC >  lfc_up)
      down <- subset(tt, adj.P.Val < f & logFC <  lfc_down)
      fwrite(up,   file.path(outdir, paste0(stem, "_UP_FDR", f, "_logFC", lfc_up,  "_", make.names(ct), ".tsv")), sep="\t")
      fwrite(down, file.path(outdir, paste0(stem, "_DN_FDR", f, "_logFC", lfc_down, "_", make.names(ct), ".tsv")), sep="\t")
    }
  }
}

# ---- Plot Gene expression in Heatmap -----

collect_sets <- function(fit, fdr=0.05, lfc_up=LFC_UP, lfc_down=LFC_DOWN) {
  ups <- list(); dns <- list()
  for (coef_i in seq_len(ncol(fit))) {
    ct <- colnames(fit)[coef_i]
    tt <- topTable(fit, coef=coef_i, number=Inf, sort.by="P")
    tt$Gene <- rownames(tt)
    ups[[ct]] <- subset(tt, adj.P.Val < fdr & logFC >  lfc_up)$Gene
    dns[[ct]] <- subset(tt, adj.P.Val < fdr & logFC <  lfc_down)$Gene
  }
  list(up=ups, down=dns)
}

plot_heatmap_union <- function(expr, meta, gene_sets, outfile, row_k = 3) {
  # ---- 1) Union of DEGs (up + down) ----
  genes <- unique(c(unlist(gene_sets$up), unlist(gene_sets$down)))
  genes <- genes[genes %in% rownames(expr)]
  if (length(genes) < 2) {
    message("No genes for heatmap: ", outfile)
    return(invisible(NULL))
  }
  
  # ---- 2) Z-score per gene ----
  mat  <- as.matrix(expr[genes, , drop = FALSE])
  matz <- t(scale(t(mat)))
  matz[is.na(matz)] <- 0
  
  # ---- 3) Column order: SLE (left) → Healthy (right); timepoints left→right ----
  time_levels <- c("<16 weeks","16-23 weeks","24-31 weeks","32-40 weeks","PP")
  meta[[TIME]]      <- factor(meta[[TIME]], levels = time_levels)
  meta[[CONDITION]] <- factor(meta[[CONDITION]], levels = c("SLE","Healthy"))
  
  samples  <- colnames(matz)
  meta_ord <- meta[match(samples, meta[[SAMPLE]]), , drop = FALSE]
  ord      <- order(meta_ord[[CONDITION]], meta_ord[[TIME]])
  matz     <- matz[, ord, drop = FALSE]
  meta_ord <- meta_ord[ord, , drop = FALSE]
  
  split_levels <- c(paste("SLE", time_levels), paste("Healthy", time_levels))
  col_split <- factor(paste(meta_ord[[CONDITION]], meta_ord[[TIME]]), levels = split_levels)
  
  # ---- 4) Row clusters + automatic labels ----
  set.seed(7)
  km <- stats::kmeans(matz, centers = row_k, nstart = 25)
  cl  <- km$cluster
  cl_mean <- tapply(rowMeans(matz), cl, mean)
  
  label_for <- function(m) {
    if (m >  0.15) "Up-regulated genes"
    else if (m < -0.15) "Down-regulated genes"
    else "Mixed/time-dependent genes"
  }
  cl_labels <- setNames(vapply(cl_mean, label_for, ""), names(cl_mean))
  row_split <- factor(cl_labels[as.character(cl)],
                      levels = c("Up-regulated genes",
                                 "Mixed/time-dependent genes",
                                 "Down-regulated genes"))
  
  # ---- 5) Annotations ----
  cond_cols <- c(SLE = "#d73027", Healthy = "#1a9850")
  tp_cols   <- setNames(c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"), time_levels)
  
  ann_df <- data.frame(
    Condition  = droplevels(meta_ord[[CONDITION]]),
    time_point = droplevels(meta_ord[[TIME]])
  )
  rownames(ann_df) <- meta_ord[[SAMPLE]]
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    df  = ann_df,
    col = list(Condition = cond_cols, time_point = tp_cols),
    annotation_name_side = "left",
    simple_anno_size = grid::unit(3, "mm")
  )
  
  # ---- 6) Draw heatmap ----
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  png(outfile, width = 1800, height = 1200, res = 170)
  ht <- ComplexHeatmap::Heatmap(
    matz,
    name = "Expression (z)",
    col  = col_fun,
    top_annotation = top_ann,
    column_split = col_split,
    gap = grid::unit(1.5, "mm"),
    row_split = row_split,                # labeled Up/Mixed/Down groups
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,              # preserve time order
    clustering_method_rows = "ward.D2",
    column_title = "Columns = SAMPLES (SLE → Healthy; timepoints left→right)",
    row_title    = "Rows = GENES (z-scored per gene)",
    heatmap_legend_param = list(
      at = c(-2, 0, 2),
      labels = c("Down-regulated", "Mean", "Up-regulated")
    )
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
  dev.off()
  
  message("✅ Heatmap saved: ", outfile,
          "  (genes: ", nrow(matz), ", samples: ", ncol(matz), ")")
}

# ---- Plot Gene expression in Upset Plots -----

plot_upset <- function(set_list, outfile, title) {
  non_empty <- set_list[lengths(set_list) > 0]
  if (length(non_empty) < 2) {
    message("Skip UpSet (need ≥2 non-empty sets): ", outfile)
    return(invisible(NULL))
  }
  png(outfile, width=1400, height=900, res=150)
  p <- UpSetR::upset(
    UpSetR::fromList(non_empty),
    order.by = "freq",
    mainbar.y.label = title,
    sets.x.label    = "Set size",
    nsets           = min(length(non_empty), 15),
    nintersects     = 30
  )
  print(p)  # ensure it renders to the device
  dev.off()
}

# ---- Export Differential Expressed Gene Lists for the ML next steps -----

export_DEG_features <- function(expr, sets_A, sets_B, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  deg_union <- unique(c(unlist(sets_A$up), unlist(sets_A$down), unlist(sets_B$up), unlist(sets_B$down)))
  deg_union <- deg_union[deg_union %in% rownames(expr)]
  write.table(data.frame(Gene=deg_union), file.path(outdir, "AllDEGs.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  vs_matrix <- as.matrix(expr) # already variance-stabilized (log2 microarray)
  rv <- rowVars(vs_matrix)
  cutoff <- quantile(rv, 0.95, na.rm=TRUE)
  hi_var <- names(rv)[rv >= cutoff]
  write.table(data.frame(Gene=hi_var), file.path(outdir, "HighVariance_genes_top5pct.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  deg_hv <- intersect(deg_union, hi_var)
  write.table(data.frame(Gene=deg_hv), file.path(outdir, "DEG_HighVariance.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  
  ml_all <- t(vs_matrix[deg_union, , drop=FALSE])
  ml_hv  <- t(vs_matrix[deg_hv,  , drop=FALSE])
  ml_hv_z<- scale(ml_hv)
  write.table(ml_all, file.path(outdir, "ML_FeatureMatrix_AllDEGs_log2.tsv"),
              sep="\t", quote=FALSE, col.names=NA)
  write.table(ml_hv,  file.path(outdir, "ML_FeatureMatrix_DEG_HighVariance_log2.tsv"),
              sep="\t", quote=FALSE, col.names=NA)
  write.table(ml_hv_z,file.path(outdir, "ML_FeatureMatrix_DEG_HighVariance_Zscored.tsv"),
              sep="\t", quote=FALSE, col.names=NA)
}

# ---- Find Biological Significance -----

map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  suppressMessages({
    m <- clusterProfiler::bitr(symbols,
                               fromType = "SYMBOL",
                               toType   = "ENTREZID",
                               OrgDb    = org.Hs.eg.db)
  })
  unique(m$ENTREZID)
}

run_clusterprofiler_on_row_clusters <- function(outdir, background_symbols = NULL, p_cut = 0.05, q_cut = 0.10) {
  # background universe from all genes in expr, if provided
  universe <- NULL
  if (!is.null(background_symbols)) {
    universe <- map_symbols_to_entrez(background_symbols)
  }
  
  files <- list.files(outdir, pattern = "^Heatmap_RowCluster_\\d+_genes\\.tsv$", full.names = TRUE)
  if (length(files) == 0) { message("No row-cluster TSVs found in: ", outdir); return(invisible(NULL)) }
  
  for (f in files) {
    cl_lab <- sub("\\.tsv$", "", basename(f))  # e.g., Heatmap_RowCluster_3_genes
    df <- data.table::fread(f)
    gene_symbols <- unique(df$Gene)
    genes_entrez <- map_symbols_to_entrez(gene_symbols)
    if (length(genes_entrez) < 5) { message("Skip ", cl_lab, " (too few mapped genes)"); next }
    
    # GO (BP/CC/MF)
    for (ont in c("BP","CC","MF")) {
      ego <- enrichGO(gene          = genes_entrez,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = ont,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = p_cut,
                      qvalueCutoff  = q_cut,
                      universe      = universe,
                      readable      = TRUE)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
        data.table::fwrite(as.data.frame(ego), file.path(outdir, paste0(cl_lab, "_GO_", ont, ".tsv")), sep="\t")
        png(file.path(outdir, paste0(cl_lab, "_GO_", ont, "_dotplot.png")), width=1400, height=1000, res=150)
        print(enrichplot::dotplot(ego, showCategory = 20))
        dev.off()
      }
    }
    
    # KEGG (needs human 'hsa')
    ek <- enrichKEGG(gene          = genes_entrez,
                     organism      = "hsa",
                     pvalueCutoff  = p_cut,
                     pAdjustMethod = "BH",
                     qvalueCutoff  = q_cut,
                     universe      = universe)
    if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
      data.table::fwrite(as.data.frame(ek), file.path(outdir, paste0(cl_lab, "_KEGG.tsv")), sep="\t")
      png(file.path(outdir, paste0(cl_lab, "_KEGG_dotplot.png")), width=1400, height=1000, res=150)
      print(enrichplot::dotplot(ek, showCategory = 20))
      dev.off()
    }
    
    # Reactome
    er <- ReactomePA::enrichPathway(gene         = genes_entrez,
                                    organism     = "human",
                                    pvalueCutoff = p_cut,
                                    pAdjustMethod= "BH",
                                    qvalueCutoff = q_cut,
                                    universe     = universe,
                                    readable     = TRUE)
    if (!is.null(er) && nrow(as.data.frame(er)) > 0) {
      data.table::fwrite(as.data.frame(er), file.path(outdir, paste0(cl_lab, "_Reactome.tsv")), sep="\t")
      png(file.path(outdir, paste0(cl_lab, "_Reactome_dotplot.png")), width=1400, height=1000, res=150)
      print(enrichplot::dotplot(er, showCategory = 20))
      dev.off()
    }
  }
  message("Enrichment completed for row clusters in: ", outdir)
}

# ---- SET UP EVERYTHIN IN A FUNCTION -------------
# ================================

deg_script <- function() {
  message("[1/10] Loading inputs ...")
  dat <- load_inputs(EXPR_FILE, META_FILE)
  expr <- dat$expr
  meta <- dat$meta
  
  message("[2/10] Building design matrix and fitting base model ...")
  mod <- build_design(expr, meta)
  design <- mod$design
  fit0 <- mod$fit0
  lookup <- make_column_lookup(design, meta)
  
  message("[3/10] Running conditional DEG analysis (SLE vs Healthy per timepoint) ...")
  fit_A <- run_conditional_DE(fit0, design, meta, lookup)
  
  message("[4/10] Saving DEG tables for conditional contrasts ...")
  save_DE_tables(fit_A, "A", file.path(OUTDIR, "A_SLE_vs_Healthy"), fdrs=FDRS)
  sets_A_005 <- collect_sets(fit_A, fdr=0.05)
  
  message("[5/10] Running temporal DEG analysis (within condition, mode = ", TEMPORAL_MODE, ") ...")
  fit_B <- run_temporal_DE(fit0, design, meta, lookup, mode=TEMPORAL_MODE)
  
  message("[6/10] Saving DEG tables for temporal contrasts ...")
  save_DE_tables(fit_B, "B", file.path(OUTDIR, "B_Temporal_" %+% TEMPORAL_MODE), fdrs=FDRS)
  sets_B_005 <- collect_sets(fit_B, fdr=0.05)
  
  # sanity log for set sizes before plotting
  log_set_sizes(sets_A_005$up,   "A-UP @ FDR 0.05")
  log_set_sizes(sets_A_005$down, "A-DOWN @ FDR 0.05")
  log_set_sizes(sets_B_005$up,   "B-UP @ FDR 0.05")
  log_set_sizes(sets_B_005$down, "B-DOWN @ FDR 0.05")
  
  message("[7/10] Creating ComplexHeatmap (SLE left, Healthy right) ...")
  plot_heatmap_union(
    expr, meta, sets_A_005,
    file.path(OUTDIR, "Heatmap_A_DEG_Union_FDR0.05.png")
  )
  
  # (Optional but recommended) Run enrichment on exported row clusters.
  # Use the full expression gene set as background universe for proper multiple testing.
  message("[7b] Running ClusterProfiler on heatmap row clusters ...")
  run_clusterprofiler_on_row_clusters(
    outdir = file.path(OUTDIR, "Enrichment"),
    background_symbols = rownames(expr),  # universe = all profiled genes
    p_cut = 0.05,
    q_cut = 0.10
  )
  
  message("[8/10] Generating UpSet plots for overlaps ...")
  plot_upset(sets_A_005$up,   file.path(OUTDIR, "UpSet_A_UP_FDR0.05.png"),   "Intersections (A-UP)")
  plot_upset(sets_A_005$down, file.path(OUTDIR, "UpSet_A_DN_FDR0.05.png"),   "Intersections (A-DOWN)")
  plot_upset(sets_B_005$up,   file.path(OUTDIR, "UpSet_B_UP_FDR0.05.png"),   "Intersections (B-UP)")
  plot_upset(sets_B_005$down, file.path(OUTDIR, "UpSet_B_DN_FDR0.05.png"),   "Intersections (B-DOWN)")
  
  message("[9/10] Exporting DE genes for ML feature matrices and lists ...")
  export_DEG_features(expr, sets_A_005, sets_B_005, file.path(OUTDIR, "DE_Genes"))
  
  message("[10/10] Performing global Time × Condition F-test (TiSA-style) ...")
  
  # ---- Build clean factors (variables), NOT meta[[...]] inside the formula
  df_TI <- data.frame(
    Condition = droplevels(meta[[CONDITION]]),
    Time      = droplevels(meta[[TIME]])
  )
  
  # Optional: flip reference if you prefer Healthy vs SLE naming
  # df_TI$Condition <- relevel(df_TI$Condition, ref = "Healthy")
  
  # Diagnostics: show balance table
  message("Sample counts per Condition × Time:")
  print(addmargins(table(df_TI$Condition, df_TI$Time)))
  
  # ---- Correct design with interaction
  design_TI <- model.matrix(~ 0 + Condition + Time + Condition:Time, data = df_TI)
  colnames(design_TI) <- make.names(colnames(design_TI))
  
  # Diagnostics: show what was created
  message("Design columns (first 10): ", paste(head(colnames(design_TI), 10), collapse = ", "))
  assign_vec  <- attr(design_TI, "assign")
  term_labels <- attr(attr(design_TI, "terms"), "term.labels")
  message("Term labels: ", paste(term_labels, collapse = " | "))
  message("Assign counts by term index: ", paste(table(assign_vec), collapse = " | "))
  
  # ---- Fit
  fit0_TI <- lmFit(
    expr, design_TI,
    block       = meta[[SUBJECT]],
    correlation = mod$corfit$consensus.correlation
  ) |> eBayes()
  
  # ---- Find interaction columns (robust primary method)
  int_term_idx <- which(term_labels == "Condition:Time")
  int_col_idx  <- which(assign_vec == int_term_idx)
  int_cols     <- colnames(design_TI)[int_col_idx]
  
  # ---- Fallback (name-based) in case attributes vary across setups
  if (length(int_cols) == 0) {
    int_cols <- grep("^Condition.*\\.Time", colnames(design_TI), value = TRUE)
  }
  message("Detected interaction columns: ", ifelse(length(int_cols)==0, "<none>", paste(int_cols, collapse=", ")))
  
  if (length(int_cols) == 0) {
    stop("No interaction columns detected. This means either (a) the model formula wasn’t used, ",
         "(b) a factor has a single level, or (c) some Condition×Time cells are empty. ",
         "Check the balance table above.")
  }
  
  # ---- Joint F-test across interaction coefficients
  # Build full-size contrast matrix (10 rows, 4 contrasts)
  L <- matrix(0, nrow = ncol(design_TI), ncol = length(int_cols),
              dimnames = list(colnames(design_TI), int_cols))
  for (i in seq_along(int_cols)) {
    L[int_cols[i], i] <- 1
  }
  
  # Apply contrasts and compute F-test
  # Apply contrasts and compute F-test
  fit_int <- contrasts.fit(fit0_TI, L)
  fit_int <- eBayes(fit_int)
  
  # Modern limma call (instead of topTableF)
  tt_inter <- topTable(fit_int, number = Inf, sort.by = "F")
  data.table::fwrite(tt_inter,
                     file.path(OUTDIR, "Global_Time_by_Condition_Ftest.tsv"),
                     sep = "\t")
  
}  

# ---- RUN COMMAND --------------------
deg_script()

