# User's test code with debugging and fixes
# Run this in R session after: cd /home/user3/GJC_KDW_250721 && R

cat("=== CCI Analysis Test (User Code) ===\n\n")

# 1. Load packages
cat("1. Loading packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
})
cat("  ✓ Packages loaded\n\n")

# 2. Load functions
cat("2. Loading CCI functions...\n")
devtools::load_all("/home/user3/data_user3/git_repo/_wt/cci/myR")

# Check if marker_filter is available
if (!exists("marker_filter")) {
  # Try to source it from utils_markers.R
  if (file.exists("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_markers.R")) {
    source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_markers.R")
    cat("  ✓ marker_filter loaded from utils_markers.R\n")
  } else {
    warning("marker_filter function not found. Proceeding without filtering.")
  }
}

# Source CCI.R (worktree version with optimizations)
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R")
cat("  ✓ CCI.R loaded\n")

# Source CCI module functions
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
cat("  ✓ All CCI functions loaded\n\n")

# 3. Check if is5 object exists
cat("3. Checking data object...\n")
if (!exists("is5")) {
  cat("  - is5 object not found. Loading from file...\n")
  if (file.exists("/data/user3/sobj/IS6_sex_added_251110.qs")) {
    is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")
    cat("  ✓ is5 loaded from file\n")
  } else {
    stop("is5 object not found and cannot load from /data/user3/sobj/IS6_sex_added_251110.qs")
  }
} else {
  cat("  ✓ is5 object found in environment\n")
}
cat("  - Cells: ", ncol(is5), ", Genes: ", nrow(is5), "\n\n")

# 4. Check metadata
cat("4. Checking metadata...\n")
if (!"anno3.scvi" %in% colnames(is5@meta.data)) {
  stop("anno3.scvi column not found in is5 metadata")
}
if (!"g3" %in% colnames(is5@meta.data)) {
  stop("g3 column not found in is5 metadata")
}
cat("  ✓ Required metadata columns found\n\n")

# 5. Prepare DEG list
cat("5. Preparing DEG list...\n")
cat("  - Running FindMarkers for Platelets / Megakaryocytes...\n")

# Check if PrepSCTFindMarkers is needed
if ("SCT" %in% names(is5@assays)) {
  tryCatch({
    Seurat::PrepSCTFindMarkers(is5)
    cat("  ✓ PrepSCTFindMarkers completed\n")
  }, error = function(e) {
    cat("  - PrepSCTFindMarkers warning: ", e$message, "\n")
  })
}

# Check if Platelets / Megakaryocytes cluster exists
platelet_cells <- sum(is5@meta.data$anno3.scvi == "Platelets / Megakaryocytes", na.rm = TRUE)
if (platelet_cells == 0) {
  stop("No cells found for 'Platelets / Megakaryocytes' cluster")
}
cat("  - Cells in cluster: ", platelet_cells, "\n")

# Check conditions
platelet_meta <- is5@meta.data[is5@meta.data$anno3.scvi == "Platelets / Megakaryocytes", ]
g3_counts <- table(platelet_meta$g3, useNA = "ifany")
cat("  - Cells per condition: ", paste(names(g3_counts), "=", g3_counts, collapse = ", "), "\n")

# Run FindMarkers
deg_df <- tryCatch({
  result <- Seurat::FindMarkers(
    subset(is5, anno3.scvi == "Platelets / Megakaryocytes"),
    group.by = "g3",
    ident.1 = "2",
    assay = "RNA",  # Use RNA assay to avoid SCT issues
    min.pct = 0.1,
    logfc.threshold = 0
  )
  
  # Add gene column if not present
  if (!"gene" %in% colnames(result)) {
    result$gene <- rownames(result)
  }
  
  # Apply marker_filter if available
  if (exists("marker_filter")) {
    result <- marker_filter(result)
  } else {
    cat("  - Warning: marker_filter not available, skipping gene filtering\n")
  }
  
  # Add cluster column
  result <- result %>%
    mutate(cluster = "Platelets / Megakaryocytes")
  
  result
}, error = function(e) {
  cat("  ✗ FindMarkers failed: ", e$message, "\n")
  stop("FindMarkers failed. Please check the error above.")
})

cat("  ✓ DEG list prepared: ", nrow(deg_df), " genes\n")
cat("  - Columns: ", paste(colnames(deg_df), collapse = ", "), "\n\n")

# 6. Set output directory
cat("6. Setting output directory...\n")
if (!exists("output_dir")) {
  output_dir <- "/data/user3/sobj/cci_plots_output"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  ✓ Created output directory: ", output_dir, "\n")
  } else {
    cat("  ✓ Using existing output directory: ", output_dir, "\n")
  }
} else {
  cat("  ✓ Using provided output_dir: ", output_dir, "\n")
}
cat("\n")

# 7. Run CCI analysis
cat("7. Running CCI analysis...\n")
cat("  This may take a while...\n\n")

CCI_results <- tryCatch({
  run_cci_analysis(
    sobj = is5,
    cluster_col = "anno3.scvi",
    deg_df = deg_df,
    receiver_cluster = "Platelets / Megakaryocytes",
    condition_col = "g3",
    condition_oi = "2",
    condition_ref = "1",
    run_circos = TRUE,
    circos_show_legend = TRUE,
    circos_legend_position = "topright",
    output_dir = output_dir,
    p_val_adj_cutoff = 1.1,
    nichenet_data_dir = "/data/user3/git_repo/human",  # Explicitly set NicheNet data path
    verbose = TRUE
  )
}, error = function(e) {
  cat("\n  ✗ CCI analysis failed: ", e$message, "\n")
  cat("  Error details:\n")
  print(e)
  stop("CCI analysis failed")
})

cat("\n  ✓ CCI analysis completed successfully!\n\n")

# 8. Check results
cat("8. Checking results...\n")
cat("  - Result components: ", paste(names(CCI_results), collapse = ", "), "\n")
if (!is.null(CCI_results$nichenet_results)) {
  cat("  - Top ligands: ", length(CCI_results$nichenet_results$best_upstream_ligands), "\n")
  if (length(CCI_results$nichenet_results$best_upstream_ligands) > 0) {
    cat("    - First 5: ", paste(head(CCI_results$nichenet_results$best_upstream_ligands, 5), collapse = ", "), "\n")
  }
}
cat("  - Output directory: ", CCI_results$output_path, "\n")
cat("\n")

cat("=== Test Complete ===\n")
cat("Results saved in: ", output_dir, "\n")

