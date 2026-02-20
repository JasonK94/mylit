#!/usr/bin/env Rscript
# Run CCI analysis for all receiver clusters in IS5 original data
# Save results to sobj/cci/whole/

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

cat("=== IS5 Whole Dataset CCI Analysis ===\n\n")

# 1. Load IS5 original data
cat("1. Loading IS5 original data...\n")
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

cat("  Cells: ", ncol(is5), ", Genes: ", nrow(is5), "\n")

# 2. Setup
cat("\n2. Setup...\n")
meta_cols <- colnames(is5@meta.data)
cluster_col <- if ("anno3.scvi" %in% meta_cols) "anno3.scvi" else "seurat_clusters"
condition_col <- if ("g3" %in% meta_cols) "g3" else NULL

cat("  Cluster column:", cluster_col, "\n")
cat("  Condition column:", condition_col, "\n")

# Get all clusters
clusters <- unique(is5@meta.data[[cluster_col]])
clusters <- clusters[!is.na(clusters)]
clusters <- sort(clusters)
cat("  Total clusters:", length(clusters), "\n")
cat("  Clusters:", paste(clusters, collapse = ", "), "\n")

# Setup output directory
output_base <- "/data/user3/sobj/cci/whole"
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
cat("  Output directory:", output_base, "\n")

# 3. Load marker_filter if available
if (!exists("marker_filter")) {
  if (file.exists("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_markers.R")) {
    source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/utils_markers.R")
  } else if (file.exists("/data/user3/git_repo/mylit/myR/R/utils_markers.R")) {
    source("/data/user3/git_repo/mylit/myR/R/utils_markers.R")
  }
}

# 4. Run analysis for each receiver cluster
cat("\n4. Running CCI analysis for each receiver cluster...\n")
cat("  Total receivers to process: ", length(clusters), "\n\n")

results_summary <- list()
failed_clusters <- character(0)

for (i in seq_along(clusters)) {
  receiver_cluster <- clusters[i]
  cat("\n", "=", rep("=", 60), "\n", sep = "")
  cat("Receiver ", i, "/", length(clusters), ": ", receiver_cluster, "\n", sep = "")
  cat("=", rep("=", 60), "\n\n", sep = "")
  
  # Get sender clusters (all except receiver)
  sender_clusters <- clusters[clusters != receiver_cluster]
  
  cat("  Receiver:", receiver_cluster, "\n")
  cat("  Senders:", length(sender_clusters), "clusters\n")
  
  # Create output directory for this receiver
  receiver_output_dir <- file.path(output_base, paste0("receiver_", receiver_cluster))
  dir.create(receiver_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Try to load pre-computed DEGs from consensus results
  deg_df <- NULL
  # Try multiple patterns for consensus file names
  cluster_name_clean <- tolower(gsub("[^A-Za-z0-9]", "_", receiver_cluster))
  # Remove common prefixes/suffixes and try variations
  cluster_name_short <- gsub("^platelets_|_megakaryocytes$|_cells$", "", cluster_name_clean)
  
  consensus_patterns <- c(
    paste0("deg_consensus_ds_*_final_result_", cluster_name_clean, ".qs"),
    paste0("deg_consensus_ds_*_final_result_", cluster_name_short, ".qs"),
    paste0("*_final_result_", cluster_name_clean, ".qs"),
    paste0("*_final_result_", cluster_name_short, ".qs")
  )
  
  consensus_files <- character(0)
  for (pattern in consensus_patterns) {
    consensus_path <- file.path("/data/user3/sobj/consensus/final", pattern)
    found_files <- Sys.glob(consensus_path)
    if (length(found_files) > 0) {
      consensus_files <- c(consensus_files, found_files)
      break
    }
  }
  
  if (length(consensus_files) > 0) {
    cat("  Loading pre-computed DEGs from:", consensus_files[1], "\n")
    tryCatch({
      deg_consensus <- qs::qread(consensus_files[1])
      if (!is.null(deg_consensus$consensus_scores)) {
        deg_df <- deg_consensus$consensus_scores %>%
          dplyr::mutate(
            p_val = meta_p,
            p_val_adj = meta_p_adj,
            avg_log2FC = mean_beta,
            cluster = receiver_cluster
          )
        cat("  Found ", nrow(deg_df), " pre-computed DEGs\n")
      }
    }, error = function(e) {
      cat("  Failed to load consensus DEGs:", e$message, "\n")
    })
  }
  
  # If no pre-computed DEGs, compute using FindMarkers
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    cat("  Computing DEGs using FindMarkers...\n")
    tryCatch({
      deg_df <- Seurat::FindMarkers(
        subset(is5, !!sym(cluster_col) == receiver_cluster),
        group.by = condition_col,
        ident.1 = "2"
      )
      if (exists("marker_filter")) {
        deg_df <- marker_filter(deg_df)
      }
      deg_df <- deg_df %>%
        dplyr::mutate(cluster = receiver_cluster)
      cat("  Found ", nrow(deg_df), " DEGs\n")
    }, error = function(e) {
      cat("  Error computing DEGs:", e$message, "\n")
      deg_df <- NULL
    })
  }
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    cat("  ✗ Cannot proceed without DEGs. Skipping this receiver.\n")
    failed_clusters <- c(failed_clusters, receiver_cluster)
    results_summary[[receiver_cluster]] <- list(
      status = "failed",
      error = "No DEGs available"
    )
    next
  }
  
  start_time <- Sys.time()
  
  tryCatch({
    results <- run_cci_analysis(
      sobj = is5,
      cluster_col = cluster_col,
      deg_df = deg_df,
      receiver_cluster = receiver_cluster,
      sender_clusters = sender_clusters,
      condition_col = condition_col,
      condition_oi = if (!is.null(condition_col)) "2" else NULL,
      condition_ref = if (!is.null(condition_col)) "1" else NULL,
      species = "human",
      nichenet_data_dir = "/data/user3/git_repo/human",
      output_dir = receiver_output_dir,
      verbose = TRUE,
      auto_save = TRUE,
      save_prepared_data = FALSE
    )
    
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat("\n  ✓ Analysis completed successfully!\n")
    cat("  Time: ", round(elapsed, 1), " minutes\n")
    cat("  Output: ", receiver_output_dir, "\n")
    
    # Check results
    if (!is.null(results$nichenet_results)) {
      cat("  ✓ NicheNet results present\n")
      if (inherits(results$nichenet_results$plot_circos, "recordedplot")) {
        cat("  ✓ Circos plot recorded\n")
      }
      if (!is.null(results$nichenet_results$plot_ligand_target_network)) {
        cat("  ✓ Ligand-target heatmap present\n")
      }
    }
    
    results_summary[[receiver_cluster]] <- list(
      status = "success",
      time_minutes = as.numeric(elapsed),
      output_dir = receiver_output_dir,
      n_ligands = length(results$nichenet_results$best_upstream_ligands),
      n_lr_pairs = nrow(results$nichenet_results$ligand_receptor_network_df)
    )
    
  }, error = function(e) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat("\n  ✗ Analysis failed after ", round(elapsed, 1), " minutes\n")
    cat("  Error: ", e$message, "\n")
    
    failed_clusters <- c(failed_clusters, receiver_cluster)
    results_summary[[receiver_cluster]] <- list(
      status = "failed",
      time_minutes = as.numeric(elapsed),
      error = e$message
    )
  })
}

# 5. Summary
cat("\n\n", "=", rep("=", 60), "\n", sep = "")
cat("Analysis Summary\n")
cat("=", rep("=", 60), "\n\n", sep = "")

successful <- sum(sapply(results_summary, function(x) x$status == "success"))
failed <- length(failed_clusters)

cat("Total receivers processed: ", length(clusters), "\n")
cat("Successful: ", successful, "\n")
cat("Failed: ", failed, "\n")

if (failed > 0) {
  cat("\nFailed clusters:\n")
  for (fc in failed_clusters) {
    cat("  - ", fc, "\n", sep = "")
  }
}

# Save summary
summary_path <- file.path(output_base, "analysis_summary.qs")
qs::qsave(results_summary, summary_path)
cat("\nSummary saved to: ", summary_path, "\n")

cat("\n=== Complete ===\n")

