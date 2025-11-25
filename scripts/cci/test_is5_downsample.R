#!/usr/bin/env Rscript
# Test CCI analysis with heavily downsampled IS5 data
# Downsample: max 10 cells per patient per cluster

# Set environment
if (file.exists("/home/user3/GJC_KDW_250721/start.R")) {
  source("/home/user3/GJC_KDW_250721/start.R")
}

library(Seurat)
library(dplyr)
library(qs)

# Source CCI functions
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")

cat("=== IS5 Heavy Downsampling Test ===\n\n")

# 1. Load IS5 data
cat("1. Loading IS5 data...\n")
is5_paths <- c(
  "/data/user3/sobj/IS_scvi_251107.qs",
  "/data/user3/sobj/IS5_g3NA_removal_251110.qs",
  "/data/user3/sobj/IS6_sex_added_251110.qs"
)

is5 <- NULL
for (path in is5_paths) {
  if (file.exists(path)) {
    cat("  Loading from:", path, "\n")
    is5 <- qs::qread(path)
    break
  }
}

if (is.null(is5)) {
  stop("IS5 data not found in any of the expected paths")
}

cat("  Original cells: ", ncol(is5), ", Genes: ", nrow(is5), "\n")

# 2. Check metadata columns
cat("\n2. Checking metadata...\n")
meta_cols <- colnames(is5@meta.data)
cat("  Available columns:", paste(head(meta_cols, 10), collapse = ", "), "...\n")

# Find patient and cluster columns
patient_col <- NULL
cluster_col <- NULL

if ("hos_no" %in% meta_cols) {
  patient_col <- "hos_no"
} else if ("patient" %in% meta_cols) {
  patient_col <- "patient"
} else if ("sample" %in% meta_cols) {
  patient_col <- "sample"
}

if ("anno3.scvi" %in% meta_cols) {
  cluster_col <- "anno3.scvi"
} else if ("seurat_clusters" %in% meta_cols) {
  cluster_col <- "seurat_clusters"
}

if (is.null(patient_col) || is.null(cluster_col)) {
  stop("Required metadata columns not found. Patient: ", patient_col, ", Cluster: ", cluster_col)
}

cat("  Using patient column:", patient_col, "\n")
cat("  Using cluster column:", cluster_col, "\n")

# 3. Heavy downsampling: max 10 cells per patient per cluster
cat("\n3. Heavy downsampling (max 10 cells per patient per cluster)...\n")
meta <- is5@meta.data
meta$cell_id <- colnames(is5)

# Group by patient and cluster, sample max 10 cells
set.seed(42)
downsampled_cells <- meta %>%
  dplyr::group_by(!!sym(patient_col), !!sym(cluster_col)) %>%
  dplyr::sample_n(min(10, n())) %>%
  dplyr::pull(cell_id)

is5_ds <- is5[, downsampled_cells]
cat("  Downsampled cells: ", ncol(is5_ds), "\n")
cat("  Reduction: ", round((1 - ncol(is5_ds)/ncol(is5)) * 100, 1), "%\n")

# 4. Save downsampled data
ds_output_path <- "/data/user3/sobj/IS5_ds10_per_patient_cluster.qs"
cat("\n4. Saving downsampled data to:", ds_output_path, "\n")
qs::qsave(is5_ds, ds_output_path)
cat("  ✓ Saved\n")

# 5. Test CCI analysis
cat("\n5. Testing CCI analysis...\n")
clusters <- unique(is5_ds@meta.data[[cluster_col]])
clusters <- clusters[!is.na(clusters)]
receiver_cluster <- clusters[1]

cat("  Receiver cluster:", receiver_cluster, "\n")
cat("  Total clusters:", length(clusters), "\n")

# Get sender clusters (all except receiver)
sender_clusters <- clusters[clusters != receiver_cluster]
cat("  Sender clusters:", length(sender_clusters), "\n")

# Compute DEGs using FindMarkers (method 1)
cat("\n  Computing DEGs using FindMarkers...\n")
if (!exists("marker_filter")) {
  # Try to source marker_filter
  if (file.exists("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_markers.R")) {
    source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_markers.R")
  } else if (file.exists("/data/user3/git_repo/mylit/myR/R/utils_markers.R")) {
    source("/data/user3/git_repo/mylit/myR/R/utils_markers.R")
  }
}

# Use FindMarkers with subset
deg_df <- tryCatch({
  Seurat::FindMarkers(
    subset(is5_ds, !!sym(cluster_col) == receiver_cluster),
    group.by = "g3",
    ident.1 = "2"
  ) %>%
    {if(exists("marker_filter")) marker_filter(.) else .} %>%
    dplyr::mutate(cluster = receiver_cluster)
}, error = function(e) {
  cat("  Error computing DEGs:", e$message, "\n")
  cat("  Trying alternative method...\n")
  NULL
})

if (is.null(deg_df) || nrow(deg_df) == 0) {
  stop("Failed to compute DEGs. Please check data and conditions.")
}

cat("  Found ", nrow(deg_df), " DEGs\n")
cat("  DEG columns:", paste(colnames(deg_df), collapse = ", "), "\n")
if ("avg_log2FC" %in% colnames(deg_df)) {
  cat("  avg_log2FC range:", range(deg_df$avg_log2FC, na.rm = TRUE), "\n")
  cat("  Upregulated (avg_log2FC > 0.1):", sum(deg_df$avg_log2FC > 0.1, na.rm = TRUE), "\n")
}
if ("p_val_adj" %in% colnames(deg_df)) {
  cat("  p_val_adj < 0.1:", sum(deg_df$p_val_adj < 0.1, na.rm = TRUE), "\n")
}

# Run CCI analysis
start_time <- Sys.time()
tryCatch({
  results <- run_cci_analysis(
    sobj = is5_ds,
    cluster_col = cluster_col,
    deg_df = deg_df,
    receiver_cluster = receiver_cluster,
    sender_clusters = sender_clusters,
    condition_col = if ("g3" %in% meta_cols) "g3" else NULL,
    condition_oi = if ("g3" %in% meta_cols) "2" else NULL,
    condition_ref = if ("g3" %in% meta_cols) "1" else NULL,
    species = "human",
    nichenet_data_dir = "/data/user3/git_repo/human",
    output_dir = "/data/user3/sobj/cci/test_ds10",
    verbose = TRUE,
    top_n_targets_per_ligand = 50,
    p_val_adj_cutoff = 1.0,  # Very lenient for downsampled data (effectively no filter)
    logfc_cutoff = 0.1  # More lenient for downsampled data
  )
  
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat("\n  ✓ CCI analysis completed successfully!\n")
  cat("  Time elapsed: ", round(elapsed, 1), " seconds\n")
  
  # Check if results are valid
  if (!is.null(results$nichenet_results)) {
    cat("  ✓ NicheNet results present\n")
    if (!is.null(results$nichenet_results$plot_circos)) {
      if (inherits(results$nichenet_results$plot_circos, "recordedplot")) {
        cat("  ✓ Circos plot recorded (length: ", length(results$nichenet_results$plot_circos), ")\n")
      } else {
        cat("  ⚠ Circos plot is not a recordedplot\n")
      }
    }
    if (!is.null(results$nichenet_results$plot_ligand_target_network)) {
      cat("  ✓ Ligand-target heatmap present\n")
    }
  }
  
}, error = function(e) {
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat("\n  ✗ CCI analysis failed after ", round(elapsed, 1), " seconds\n")
  cat("  Error: ", e$message, "\n")
  traceback()
  stop("Test failed")
})

cat("\n=== Test Complete ===\n")

