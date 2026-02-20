#!/usr/bin/env Rscript
# ==============================================================================
# MILO Neighbourhood DA Analysis
# Layer 1: HC vs IS (cohort)
# Layer 2: g3==1 vs g3==2 (IS only)
# ==============================================================================

.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Try to load myR milo pipeline
tryCatch({
  source("/home/user3/data_user3/git_repo/mylit/myR/R/analysis/milo_pipeline.R")
  cat("  Loaded: run_milo_pipeline() from myR\n")
  USE_MYR <- TRUE
}, error = function(e) {
  cat("  myR milo_pipeline.R not loadable, using inline implementation\n")
  USE_MYR <<- FALSE
})

cat("=== MILO Neighbourhood DA Analysis ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
OUT <- file.path(BASE, "milo")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: run Milo analysis ----
run_milo_analysis <- function(sobj, target_var, cluster_var = "anno1",
                               patient_var = "patient_name",
                               reduction = "integrated.scvi",
                               output_dir, label) {
  cat("\n=== Milo:", label, "===\n")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Convert to SCE
  cat("  Converting to SCE...\n")
  sce <- as.SingleCellExperiment(sobj, assay = "RNA")

  # Use scVI embedding for KNN graph
  red_mat <- Embeddings(sobj, reduction = reduction)
  reducedDim(sce, "SCVI") <- red_mat

  # Also store UMAP for plotting
  umap_mat <- Embeddings(sobj, reduction = "umap.scvi")
  reducedDim(sce, "UMAP_SCVI") <- umap_mat

  # Create Milo object
  cat("  Creating Milo object...\n")
  milo <- Milo(sce)

  # Build KNN graph on scVI latent space
  cat("  Building KNN graph (k=30)...\n")
  milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "SCVI")

  # Define neighbourhoods
  cat("  Defining neighbourhoods...\n")
  milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE,
                      reduced_dims = "SCVI")

  cat("  Nhood sizes: median =", median(colSums(nhoods(milo))),
      ", range =", paste(range(colSums(nhoods(milo))), collapse = "-"), "\n")

  # Count cells per neighbourhood per sample
  cat("  Counting cells...\n")
  milo <- countCells(milo, meta.data = data.frame(colData(milo)),
                     sample = patient_var)

  # Design matrix
  cat("  Building design matrix...\n")
  meta <- data.frame(colData(milo))

  # Get sample-level metadata
  sample_meta <- meta %>%
    group_by(across(all_of(patient_var))) %>%
    summarise(
      !!target_var := first(.data[[target_var]]),
      !!cluster_var := names(which.max(table(.data[[cluster_var]]))),
      .groups = "drop"
    ) %>%
    as.data.frame()
  rownames(sample_meta) <- sample_meta[[patient_var]]

  # Ensure target is factor
  sample_meta[[target_var]] <- factor(sample_meta[[target_var]])

  formula_str <- paste0("~ ", target_var)
  design <- model.matrix(as.formula(formula_str), data = sample_meta)

  # miloR::testNhoods internally expects "PCA" in reducedDimNames
  # Use scVI latent space as the PCA stand-in
  reducedDim(milo, "PCA") <- reducedDim(milo, "SCVI")

  # Test DA
  cat("  Testing DA (edgeR)...\n")
  milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "SCVI")
  da_results <- testNhoods(milo, design = design, design.df = sample_meta,
                            model.contrasts = colnames(design)[2])

  cat("  DA results:", nrow(da_results), "neighbourhoods\n")
  cat("  Significant (SpatialFDR < 0.1):", sum(da_results$SpatialFDR < 0.1, na.rm = TRUE), "\n")
  cat("  Significant (SpatialFDR < 0.05):", sum(da_results$SpatialFDR < 0.05, na.rm = TRUE), "\n")

  # Annotate with cell types
  da_results <- annotateNhoods(milo, da_results, coldata_col = cluster_var)

  # Save results
  write.csv(da_results, file.path(output_dir, "milo_da_results.csv"), row.names = FALSE)
  saveRDS(milo, file.path(output_dir, "milo_object.rds"))

  # ---- Plots ----
  generate_milo_plots(milo, da_results, cluster_var, output_dir, label)

  cat("  Done:", label, "\n")
  return(da_results)
}

# ---- Helper: generate plots from milo + da_results ----
generate_milo_plots <- function(milo, da_results, cluster_var, output_dir, label) {
  cat("  Generating plots...\n")

  # 1. Beeswarm plot
  p_beeswarm <- plotDAbeeswarm(da_results, group.by = cluster_var, alpha = 0.1) +
    ggtitle(paste0("MILO DA: ", label)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(output_dir, "01_beeswarm_da.pdf"), p_beeswarm,
         width = 10, height = 6)
  ggsave(file.path(output_dir, "01_beeswarm_da.png"), p_beeswarm,
         width = 10, height = 6, dpi = 200)

  # 2. UMAP with NH graph coloured by logFC
  # buildNhoodGraph required for plotNhoodGraphDA
  milo <- buildNhoodGraph(milo)
  p_nhood <- plotNhoodGraphDA(milo, da_results, layout = "UMAP_SCVI", alpha = 0.1) +
    ggtitle(paste0("Nhood DA on UMAP: ", label))

  ggsave(file.path(output_dir, "02_umap_nhood_da.png"), p_nhood,
         width = 10, height = 8, dpi = 150)

  # 3. Volcano-like: logFC vs -log10(SpatialFDR) per cell type
  da_results$neg_log10_fdr <- -log10(da_results$SpatialFDR + 1e-300)
  frac_col <- paste0(cluster_var, "_fraction")
  p_volcano <- ggplot(da_results, aes(x = logFC, y = neg_log10_fdr,
                                       color = .data[[frac_col]] > 0.7)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"),
                       labels = c("Mixed", "Homogeneous"), name = "Nhood purity") +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
    facet_wrap(as.formula(paste0("~ ", cluster_var)), scales = "free_y") +
    theme_bw() +
    ggtitle(paste0("Nhood DA per cell type: ", label))

  ggsave(file.path(output_dir, "03_volcano_per_celltype.png"), p_volcano,
         width = 16, height = 12, dpi = 150)

  # 4. Summary barplot: fraction of DA nhoods per cell type
  da_sig <- da_results %>%
    filter(.data[[frac_col]] > 0.7) %>%
    group_by(across(all_of(cluster_var))) %>%
    summarise(
      n_nhoods = n(),
      n_up = sum(logFC > 0 & SpatialFDR < 0.1),
      n_down = sum(logFC < 0 & SpatialFDR < 0.1),
      n_ns = sum(SpatialFDR >= 0.1),
      frac_up = n_up / n_nhoods,
      frac_down = n_down / n_nhoods,
      median_logFC = median(logFC),
      .groups = "drop"
    ) %>%
    arrange(median_logFC)

  write.csv(da_sig, file.path(output_dir, "04_da_summary_by_celltype.csv"),
            row.names = FALSE)

  da_long <- da_sig %>%
    tidyr::pivot_longer(cols = c(frac_up, frac_down),
                        names_to = "direction", values_to = "fraction") %>%
    mutate(
      fraction = ifelse(direction == "frac_down", -fraction, fraction),
      direction = ifelse(direction == "frac_up", "Enriched", "Depleted")
    )

  ct_col <- cluster_var
  p_bar <- ggplot(da_long, aes(x = reorder(.data[[ct_col]], median_logFC),
                                y = fraction, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Enriched" = "#E63946", "Depleted" = "#457B9D")) +
    labs(x = NULL, y = "Fraction of DA neighbourhoods",
         title = paste0("DA Summary: ", label),
         fill = "Direction") +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(file.path(output_dir, "04_da_summary_barplot.png"), p_bar,
         width = 8, height = 6, dpi = 200)
}

# ---- Layer 1: HC vs IS ----
# L1 DA was previously completed — reload milo object and regenerate plots
l1_dir <- file.path(OUT, "L1_cohort")
l1_milo_path <- file.path(l1_dir, "milo_object.rds")
l1_da_path <- file.path(l1_dir, "milo_da_results.csv")

if (file.exists(l1_milo_path) && file.exists(l1_da_path)) {
  cat("[1/2] L1 DA results exist — regenerating plots only\n")
  cat("  Loading milo object...\n")
  milo_l1 <- readRDS(l1_milo_path)
  da_l1 <- read.csv(l1_da_path)
  generate_milo_plots(milo_l1, da_l1, "anno1", l1_dir, "HC vs Stroke (Layer 1)")
  rm(milo_l1); gc(verbose = FALSE)
  cat("  Done: L1 plots regenerated\n")
} else {
  cat("[1/2] Loading Layer 1 data...\n")
  sobj_l1 <- qread(file.path(BASE, "5_1_hc_is.qs"))
  cat("  Cells:", ncol(sobj_l1), "\n")

  da_l1 <- run_milo_analysis(
    sobj_l1,
    target_var = "cohort",
    cluster_var = "anno1",
    output_dir = l1_dir,
    label = "HC vs Stroke (Layer 1)"
  )
  rm(sobj_l1); gc(verbose = FALSE)
}

# ---- Layer 2: g3 (IS only) ----
cat("\n[2/2] Loading Layer 2 data...\n")
sobj_l2 <- qread(file.path(BASE, "5_2_is_g3.qs"))
cat("  Cells:", ncol(sobj_l2), "\n")

da_l2 <- run_milo_analysis(
  sobj_l2,
  target_var = "g3",
  cluster_var = "anno1",
  output_dir = file.path(OUT, "L2_g3"),
  label = "Good vs Bad (Layer 2, g3)"
)
rm(sobj_l2); gc(verbose = FALSE)

cat("\n=== MILO COMPLETE ===\n")
cat("Finished:", format(Sys.time()), "\n")
cat("Output:", OUT, "\n")
