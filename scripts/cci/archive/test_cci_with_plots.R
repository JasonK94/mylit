#!/usr/bin/env Rscript
# Test CCI analysis with all plots including Circos plot
# Focus on Monocyte, NK cell, T cell interactions

cat("=== CCI Tool Test with All Plots ===\n\n")

# Load packages
cat("1. Loading packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
})

cat("  ✓ Packages loaded\n\n")

# Load CCI functions
cat("2. Loading CCI functions...\n")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
cat("  ✓ CCI functions loaded\n\n")

# Load run_nichenet_analysis
cat("3. Loading run_nichenet_analysis...\n")
source("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")
cat("  ✓ run_nichenet_analysis loaded\n\n")

# Load test data
cat("4. Loading test data...\n")
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
cat("  ✓ Data loaded: ", ncol(sobj_test), " cells, ", nrow(sobj_test), " genes\n\n")

# Check available clusters
cat("5. Checking available clusters...\n")
clusters <- unique(sobj_test@meta.data[["anno3.scvi"]])
clusters <- clusters[!is.na(clusters)]
clusters <- sort(clusters)

cat("  Total clusters: ", length(clusters), "\n")
cat("  Cluster IDs: ", paste(head(clusters, 10), collapse = ", "), "...\n\n")

# Find Monocyte, NK cell, T cell clusters
monocyte_clusters <- grep("mono|macrophage|mac", clusters, ignore.case = TRUE, value = TRUE)
nk_clusters <- grep("nk", clusters, ignore.case = TRUE, value = TRUE)
tcell_clusters <- grep("t-cell|t cell|tcell|cd4|cd8", clusters, ignore.case = TRUE, value = TRUE)

cat("  Monocyte-related: ", paste(monocyte_clusters, collapse = ", "), "\n")
cat("  NK cell-related: ", paste(nk_clusters, collapse = ", "), "\n")
cat("  T cell-related: ", paste(tcell_clusters, collapse = ", "), "\n\n")

# Select receiver (prefer T cell or NK cell)
if (length(tcell_clusters) > 0) {
  receiver_cluster_test <- tcell_clusters[1]
  cat("  Selected receiver: ", receiver_cluster_test, " (T cell)\n")
} else if (length(nk_clusters) > 0) {
  receiver_cluster_test <- nk_clusters[1]
  cat("  Selected receiver: ", receiver_cluster_test, " (NK cell)\n")
} else {
  receiver_cluster_test <- clusters[1]
  cat("  Selected receiver: ", receiver_cluster_test, " (first cluster)\n")
}

# Select senders (Monocyte, NK cell, T cell)
sender_candidates <- c(monocyte_clusters, nk_clusters, tcell_clusters)
sender_candidates <- sender_candidates[sender_candidates != receiver_cluster_test]
sender_candidates <- unique(sender_candidates)

if (length(sender_candidates) > 0) {
  sender_clusters_test <- head(sender_candidates, min(5, length(sender_candidates)))
  cat("  Selected senders: ", paste(sender_clusters_test, collapse = ", "), "\n\n")
} else {
  sender_clusters_test <- NULL
  cat("  Will auto-identify senders\n\n")
}

# Check conditions
cat("6. Checking conditions...\n")
g3_vals <- unique(sobj_test@meta.data$g3)
g3_vals <- g3_vals[!is.na(g3_vals)]
cat("  - g3 values: ", paste(g3_vals, collapse = ", "), "\n")

if (length(g3_vals) >= 2) {
  condition_oi <- as.character(g3_vals[1])
  condition_ref <- as.character(g3_vals[2])
  cat("  - Condition OI: ", condition_oi, "\n")
  cat("  - Condition Ref: ", condition_ref, "\n")
  cat("  ✓ Conditions suitable for analysis\n\n")
} else {
  stop("Not enough conditions for analysis")
}

# Prepare DEG list (use actual expressed genes)
cat("7. Preparing DEG list...\n")
receiver_cells <- sobj_test@meta.data[["anno3.scvi"]] == receiver_cluster_test
receiver_expr <- Seurat::GetAssayData(sobj_test, assay = "RNA", layer = "counts")[, receiver_cells]
expressed_genes <- rownames(receiver_expr)[Matrix::rowSums(receiver_expr > 0) > 0]
actual_genes <- head(expressed_genes, 30)

deg_df_example <- data.frame(
  gene = actual_genes,
  cluster = rep(receiver_cluster_test, length(actual_genes)),
  avg_log2FC = c(rep(1.5, 15), rep(-1.2, length(actual_genes) - 15)),
  p_val_adj = c(0.001, 0.0001, 0.01, 0.05, 0.1, rep(0.01, length(actual_genes) - 5)),
  stringsAsFactors = FALSE
)
cat("  ✓ DEG list prepared with ", nrow(deg_df_example), " genes\n\n")

# Set up output directory
output_dir <- "/data/user3/sobj/cci_plots_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
cat("8. Output directory: ", output_dir, "\n\n")

# Prepare Seurat object
cat("9. Preparing Seurat object...\n")
if ("SCT" %in% names(sobj_test@assays)) {
  tryCatch({
    sobj_test <- Seurat::PrepSCTFindMarkers(sobj_test)
    cat("  - PrepSCTFindMarkers completed\n")
  }, error = function(e) {
    cat("  - PrepSCTFindMarkers warning: ", e$message, "\n")
  })
}

# Use existing NicheNet data directory
nichenet_data_dir <- "/data/user3/git_repo/human"
cat("  - Using NicheNet data from: ", nichenet_data_dir, "\n\n")

# Run CCI analysis with all plots
cat("10. Running full CCI analysis with all plots...\n")
cat("    This will generate:\n")
cat("    - Ligand-Target heatmap\n")
cat("    - Ligand-Receptor heatmap\n")
cat("    - Ligand activity histogram\n")
cat("    - Ligand AUPR heatmap\n")
cat("    - Circos plot (ligand-receptor interactions)\n\n")

tryCatch({
  results <- run_cci_analysis(
    sobj = sobj_test,
    cluster_col = "anno3.scvi",
    deg_df = deg_df_example,
    receiver_cluster = receiver_cluster_test,
    sender_clusters = sender_clusters_test,  # Specify senders
    condition_col = "g3",
    condition_oi = condition_oi,
    condition_ref = condition_ref,
    species = "human",
    assay_name = "RNA",
    verbose = TRUE,
    auto_save = TRUE,
    top_n_ligands = 20,  # More ligands for better visualization
    top_n_targets_per_ligand = 200,
    run_circos = TRUE,  # Enable circos plot
    output_dir = output_dir,  # Save plots
    nichenet_data_dir = nichenet_data_dir,
    p_val_adj_cutoff = 1.1,
    logfc_cutoff = 0.05
  )
  
  cat("\n  ✓ Full CCI analysis completed!\n")
  cat("  - Top ligands: ", length(results$nichenet_results$best_upstream_ligands), "\n")
  cat("  - Sender clusters: ", length(results$sender_clusters), "\n")
  cat("  - Top 10 ligands: ", paste(head(results$nichenet_results$best_upstream_ligands, 10), collapse = ", "), "\n")
  
  if (!is.null(results$saved_path)) {
    cat("  - Results saved to: ", results$saved_path, "\n")
  }
  
  # Check output directory for plots
  if (dir.exists(output_dir)) {
    plot_files <- list.files(output_dir, pattern = "\\.(png|pdf)$", full.names = TRUE)
    cat("  - Generated plot files: ", length(plot_files), "\n")
    for (f in plot_files) {
      cat("    * ", basename(f), "\n")
    }
  }
  
}, error = function(e) {
  cat("  ✗ Full analysis failed: ", e$message, "\n")
  traceback()
})

cat("\n=== Test Complete ===\n")

