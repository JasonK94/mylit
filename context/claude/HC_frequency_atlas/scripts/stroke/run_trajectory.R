#!/usr/bin/env Rscript
# ==============================================================================
# Phase 4: Compartment-separated Trajectory Analysis
# Uses Slingshot + Monocle3 for lineage-specific pseudotime
# Followed by GAM-based gene dynamics analysis (g3 comparison)
# ==============================================================================
# Input:  5_strokev8_clean.qs (205,277 cells, all HC+IS+others)
# Output: trajectory/{mono,cd4,cd8}/ - results per compartment
#
# Design:
#   - Monocyte: CD14+ Mono → Inflammatory Mono → CD16+ Mono (± ISG+ Myeloid)
#   - CD4+ T:   Naive/Memory → Treg / S100A8_stress
#   - CD8+ T:   Cytotoxic → Trm (± MAIT)
#
# Each compartment gets:
#   1. Subset extraction
#   2. Slingshot trajectory (PCA on HVGs within subset)
#   3. Monocle3 trajectory
#   4. Gene dynamics analysis for top genes along pseudotime (g3 comparison)
# ==============================================================================

# Setup library paths FIRST (monocle3 etc. live in renv)
.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(SingleCellExperiment)
  library(slingshot)
  library(monocle3)
  library(S4Vectors)
  library(SummarizedExperiment)
  library(tibble)
  library(mgcv)
})

# Source pseudotime functions directly
PSEUDOTIME_SRC <- "/data/user3/git_repo/mylit/myR/R/analysis/pseudotime.R"
cat("Sourcing trajectory functions from:", PSEUDOTIME_SRC, "\n")
source(PSEUDOTIME_SRC)

cat("=== Phase 4: Trajectory Analysis ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ---- Configuration ----
BASE <- "/data/user3/sobj/stroke_hc_v8_2"
INPUT <- file.path(BASE, "5_strokev8_clean.qs")
TRAJ_DIR <- file.path(BASE, "trajectory")
dir.create(TRAJ_DIR, showWarnings = FALSE, recursive = TRUE)

# Top genes to analyze along pseudotime (from FGS candidates + literature)
TARGET_GENES <- c(
  # Myeloid activation / inflammation
  "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD14", "FCGR3A",
  "IL1B", "CXCL8", "CCL2", "CCL3", "TNF", "IL6",
  # ISG signature
  "ISG15", "IFI6", "IFIT1", "MX1", "OAS1",
  # T cell activation / exhaustion
  "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  # Stroke-relevant
  "DDIT4", "TXNIP", "HLA-B", "HLA-C", "CCL4", "HIF1A"
)

# ---- 1. Load data ----
cat("[1/7] Loading clean Seurat object...\n")
sobj <- qread(INPUT)
cat("  Cells:", ncol(sobj), "| Genes:", nrow(sobj), "\n")
cat("  Reductions:", paste(names(sobj@reductions), collapse = ", "), "\n")
cat("  anno1 levels:", paste(sort(unique(sobj$anno1)), collapse = ", "), "\n\n")

# ---- 2. Define compartment subsets ----
cat("[2/7] Defining compartment subsets...\n")

compartments <- list(
  mono = list(
    name = "Monocyte",
    anno1_include = c("CD14+ Monocyte", "Inflammatory Monocyte",
                      "CD16+ Monocyte", "ISG+ Myeloid"),
    start_hint = "CD14+ Monocyte",
    bio_order = c("CD14+ Monocyte", "Inflammatory Monocyte",
                  "CD16+ Monocyte", "ISG+ Myeloid"),
    description = "Classical → Inflammatory → Non-classical monocyte differentiation"
  ),
  cd4 = list(
    name = "CD4+ T cell",
    anno1_include = c("CD4+ T_Naive/Memory", "Treg",
                      "CD4_S100A8_CD14", "ISG+ T_cell"),
    start_hint = "CD4+ T_Naive/Memory",
    bio_order = c("CD4+ T_Naive/Memory", "Treg",
                  "CD4_S100A8_CD14", "ISG+ T_cell"),
    description = "Naive/Memory → Treg / stress-activated differentiation"
  ),
  cd8 = list(
    name = "CD8+ T cell",
    anno1_include = c("CD8+ T_Cytotoxic", "CD8+ Trm", "MAIT"),
    start_hint = "CD8+ T_Cytotoxic",
    bio_order = c("CD8+ T_Cytotoxic", "CD8+ Trm", "MAIT"),
    description = "Cytotoxic → Tissue-resident memory differentiation"
  )
)

# ---- 3. Function: Create subset and run preprocessing ----
create_trajectory_subset <- function(sobj, comp, out_dir) {
  cat("  Subsetting for:", comp$name, "\n")

  # Filter by anno1
  cells_keep <- colnames(sobj)[sobj$anno1 %in% comp$anno1_include]
  cat("    Cells selected:", length(cells_keep), "\n")
  cat("    anno1 composition:\n")
  tbl <- table(sobj$anno1[sobj$anno1 %in% comp$anno1_include])
  for (nm in names(tbl)) cat("      ", nm, ":", tbl[nm], "\n")

  if (length(cells_keep) < 500) {
    cat("    WARNING: Too few cells (<500). Skipping.\n")
    return(NULL)
  }

  sub <- subset(sobj, cells = cells_keep)

  # Re-normalize and find HVGs within subset
  cat("    Normalizing and finding HVGs...\n")
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 30, verbose = FALSE)

  # New UMAP based on subset PCA
  cat("    Computing subset UMAP...\n")
  sub <- RunUMAP(sub, dims = 1:20, reduction = "pca",
                 reduction.name = "umap_subset", verbose = FALSE)

  # Clustering within subset (for Slingshot)
  sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
  sub <- FindClusters(sub, resolution = 0.5, verbose = FALSE)

  # Save subset
  qsave(sub, file.path(out_dir, "subset.qs"))
  cat("    Saved subset to:", file.path(out_dir, "subset.qs"), "\n")

  return(sub)
}

# ---- 4. Function: Run Slingshot ----
run_slingshot_analysis <- function(sub, comp, out_dir) {
  cat("  Running Slingshot...\n")

  # Identify start cluster as the cluster with most cells of start_hint anno1
  start_anno <- comp$start_hint
  meta <- sub@meta.data
  meta_start <- meta[meta$anno1 == start_anno, ]
  start_cluster_candidates <- as.data.frame(
    table(cluster = meta_start$seurat_clusters),
    stringsAsFactors = FALSE
  )
  colnames(start_cluster_candidates) <- c("seurat_clusters", "n")
  start_cluster_candidates <- start_cluster_candidates[
    order(-start_cluster_candidates$n), ]

  if (nrow(start_cluster_candidates) == 0) {
    cat("    WARNING: No cells found for start_hint '", start_anno, "'. Using cluster 0.\n")
    start_cl <- "0"
  } else {
    start_cl <- as.character(start_cluster_candidates$seurat_clusters[1])
    cat("    Start cluster:", start_cl, "(most cells of", start_anno, ")\n")
  }

  # Run slingshot (use umap_subset 2D, not full PCA - much faster)
  sce <- tryCatch({
    run_slingshot_from_seurat(
      seurat_obj = sub,
      cluster_col = "seurat_clusters",
      reduced_dim_name = "umap_subset",
      start_cluster = start_cl,
      counts_assay_name = "RNA",
      approx_points = 150
    )
  }, error = function(e) {
    cat("    ERROR in Slingshot:", e$message, "\n")
    return(NULL)
  })

  if (is.null(sce)) return(NULL)

  # Extract pseudotime
  pt_cols <- grep("slingPseudotime", colnames(colData(sce)), value = TRUE)
  cat("    Slingshot lineages found:", length(pt_cols), "\n")

  # Save SCE
  saveRDS(sce, file.path(out_dir, "slingshot_sce.rds"))
  cat("    Saved Slingshot SCE to:", file.path(out_dir, "slingshot_sce.rds"), "\n")

  # Create pseudotime summary plot
  if (length(pt_cols) > 0) {
    pt_df <- as.data.frame(colData(sce))
    pt_df$UMAP1 <- sub@reductions$umap_subset@cell.embeddings[rownames(pt_df), 1]
    pt_df$UMAP2 <- sub@reductions$umap_subset@cell.embeddings[rownames(pt_df), 2]

    plots <- list()
    # Cluster annotation plot
    plots[["anno1"]] <- ggplot(pt_df, aes(x = UMAP1, y = UMAP2, color = anno1)) +
      geom_point(size = 0.3, alpha = 0.5) +
      theme_minimal() + ggtitle(paste(comp$name, "- anno1")) +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(override.aes = list(size = 2)))

    # Pseudotime plots per lineage
    for (pt_col in pt_cols) {
      p <- ggplot(pt_df, aes(x = UMAP1, y = UMAP2, color = .data[[pt_col]])) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_color_viridis_c(na.value = "grey85") +
        theme_minimal() + ggtitle(paste("Slingshot:", pt_col)) +
        theme(legend.position = "bottom")
      plots[[pt_col]] <- p
    }

    # Combine plots
    combined <- wrap_plots(plots, ncol = min(length(plots), 3))
    ggsave(file.path(out_dir, "slingshot_pseudotime_umap.pdf"),
           combined, width = min(18, 6 * length(plots)), height = 6, limitsize = FALSE)
    cat("    Saved Slingshot UMAP plot\n")
  }

  return(sce)
}

# ---- 5. Function: Run Monocle3 ----
run_monocle3_analysis <- function(sub, comp, out_dir) {
  cat("  Running Monocle3...\n")

  # Determine root cluster from Seurat anno1
  start_anno <- comp$start_hint
  start_cells <- colnames(sub)[sub$anno1 == start_anno]
  if (length(start_cells) == 0) {
    cat("    WARNING: No cells for root annotation. Using NULL.\n")
    start_cells <- NULL
  } else {
    # Sample down to avoid extremely large root set
    if (length(start_cells) > 500) {
      set.seed(42)
      start_cells <- sample(start_cells, 500)
    }
    cat("    Root cells:", length(start_cells), "from", start_anno, "\n")
  }

  cds <- tryCatch({
    run_monocle3_from_seurat(
      seurat_obj = sub,
      counts_assay_name = "RNA",
      reduction_name = "umap_subset",
      root_cells = start_cells,
      preprocess_method = "PCA",
      num_dimensions = 20,
      reduction_method = "UMAP",
      use_partition = FALSE,
      resolution = 1e-4,
      verbose = TRUE
    )
  }, error = function(e) {
    cat("    ERROR in Monocle3:", e$message, "\n")
    return(NULL)
  })

  if (is.null(cds)) return(NULL)

  # Save CDS
  saveRDS(cds, file.path(out_dir, "monocle3_cds.rds"))
  cat("    Saved Monocle3 CDS to:", file.path(out_dir, "monocle3_cds.rds"), "\n")

  # Extract pseudotime
  pt <- tryCatch(monocle3::pseudotime(cds), error = function(e) NULL)
  if (!is.null(pt)) {
    pt_finite <- pt[is.finite(pt)]
    cat("    Pseudotime range:", round(min(pt_finite), 2), "-",
        round(max(pt_finite), 2), "\n")
    cat("    Inf pseudotime cells:", sum(is.infinite(pt)), "/", length(pt), "\n")
  }

  # Plot
  tryCatch({
    p1 <- monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                                cell_size = 0.5, show_trajectory_graph = TRUE) +
      ggtitle(paste(comp$name, "- Monocle3 Pseudotime"))

    # Color by anno1 if available
    p2 <- monocle3::plot_cells(cds, color_cells_by = "anno1",
                                cell_size = 0.5, show_trajectory_graph = TRUE,
                                label_groups_by_cluster = FALSE) +
      ggtitle(paste(comp$name, "- anno1"))

    combined <- p1 + p2
    ggsave(file.path(out_dir, "monocle3_pseudotime_umap.pdf"),
           combined, width = 14, height = 6)
    cat("    Saved Monocle3 UMAP plot\n")
  }, error = function(e) {
    cat("    Plot error:", e$message, "\n")
  })

  return(cds)
}

# ---- 6. Function: Gene dynamics analysis along pseudotime (g3 comparison) ----
run_gene_dynamics <- function(cds, comp, out_dir, genes = TARGET_GENES) {
  cat("  Running gene dynamics analysis...\n")

  dyn_dir <- file.path(out_dir, "gene_dynamics")
  dir.create(dyn_dir, showWarnings = FALSE, recursive = TRUE)

  # Check which target genes exist in cds
  available_genes <- intersect(genes, rownames(cds))
  cat("    Target genes available:", length(available_genes), "/", length(genes), "\n")

  if (length(available_genes) == 0) {
    cat("    No target genes found. Skipping dynamics.\n")
    return(NULL)
  }

  # Check if g3 is available for condition comparison
  cell_meta <- as.data.frame(monocle3::pData(cds))
  has_g3 <- "g3" %in% colnames(cell_meta) && sum(!is.na(cell_meta$g3)) > 100

  if (has_g3) {
    condition_col <- "g3"
    cat("    Using g3 as condition (", sum(!is.na(cell_meta$g3)), "cells with valid g3)\n")
  } else if ("cohort" %in% colnames(cell_meta)) {
    condition_col <- "cohort"
    cat("    g3 not available. Using cohort as condition.\n")
  } else {
    cat("    No condition variable available. Skipping dynamics.\n")
    return(NULL)
  }

  # Run gene dynamics for each gene
  results <- list()
  for (gene in available_genes) {
    cat("    Analyzing:", gene, "...\n")
    res <- tryCatch({
      analyze_gene_dynamics(
        gene_id = gene,
        cds_obj = cds,
        condition_col_name = condition_col,
        output_dir = dyn_dir,
        k_val = 6,
        min_cells_for_fit = 30,
        plot_split = FALSE,
        scale_DR = TRUE
      )
    }, error = function(e) {
      cat("      Error:", e$message, "\n")
      list(gene = gene, status = "error", metrics = NULL)
    })
    results[[gene]] <- res
  }

  # Summarize
  summary_df <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r$metrics)) {
      data.frame(gene = r$gene, status = r$status,
                 TV = NA, DR = NA, interaction_pval = NA,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(
        gene = r$gene, status = r$status,
        TV = ifelse(!is.null(r$metrics$TV), r$metrics$TV, NA),
        DR = ifelse(!is.null(r$metrics$DR), r$metrics$DR, NA),
        interaction_pval = ifelse(!is.null(r$metrics$interaction_pval),
                                  r$metrics$interaction_pval, NA),
        stringsAsFactors = FALSE
      )
    }
  }))

  write.csv(summary_df, file.path(dyn_dir, "gene_dynamics_summary.csv"),
            row.names = FALSE)
  cat("    Saved gene dynamics summary:", nrow(summary_df), "genes\n")

  return(summary_df)
}

# ---- 7. Main: Run per compartment ----
cat("\n[3/7] Running trajectory analysis per compartment...\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

results_summary <- list()

for (comp_id in names(compartments)) {
  comp <- compartments[[comp_id]]
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat(" COMPARTMENT:", comp$name, "\n")
  cat(" ", comp$description, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  comp_dir <- file.path(TRAJ_DIR, comp_id)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  # Step A: Create subset
  cat("[A] Creating subset...\n")
  sub <- create_trajectory_subset(sobj, comp, comp_dir)
  if (is.null(sub)) {
    cat("  Skipping compartment due to insufficient cells.\n")
    next
  }

  # Step B: Slingshot
  cat("[B] Slingshot trajectory...\n")
  sce <- run_slingshot_analysis(sub, comp, comp_dir)

  # Step C: Monocle3
  cat("[C] Monocle3 trajectory...\n")
  cds <- run_monocle3_analysis(sub, comp, comp_dir)

  # Step D: Gene dynamics (using Monocle3 CDS with pseudotime)
  if (!is.null(cds)) {
    cat("[D] Gene dynamics...\n")
    dyn <- run_gene_dynamics(cds, comp, comp_dir)
    results_summary[[comp_id]] <- dyn
  } else {
    cat("[D] Skipping gene dynamics (no CDS).\n")
  }

  cat("\n  Compartment", comp$name, "complete.\n")
}

# ---- 8. Final summary ----
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat(" TRAJECTORY ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Finished:", format(Sys.time()), "\n")

# List all output files
cat("\nOutput directory:", TRAJ_DIR, "\n")
for (comp_id in names(compartments)) {
  comp_dir <- file.path(TRAJ_DIR, comp_id)
  if (dir.exists(comp_dir)) {
    files <- list.files(comp_dir, recursive = TRUE)
    cat("\n", comp_id, ":\n")
    for (f in files) cat("  ", f, "\n")
  }
}

cat("\nDone.\n")
