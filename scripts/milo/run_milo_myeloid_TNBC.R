#!/usr/bin/env Rscript
# Script to run Milo analysis for TNBC myeloid data
# Comparing 'low' vs 'high' for intra_T and TIL_2

library(Seurat)
library(qs)
library(miloR)

# Source the milo pipeline function from main repo
milo_pipeline_r <- "/data/kjc2/git_repo/mylit/myR/R/milo_pipeline.R"
if (file.exists(milo_pipeline_r)) {
  source(milo_pipeline_r)
} else {
  stop("Could not find milo_pipeline.R at: ", milo_pipeline_r)
}

# File paths
seurat_qs_path <- "/data/kjc2/git_repo/_wt/h5ad2sobj/myeloid_v1_TNBC.qs"
output_base_dir <- "/data/kjc2/git_repo/_wt/h5ad2sobj/milo_results"

# Check if file exists
if (!file.exists(seurat_qs_path)) {
  stop("Seurat .qs file not found: ", seurat_qs_path)
}

# Create output directory
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
message("=", rep("=", 79))
message("Loading Seurat object...")
message("=", rep("=", 79))
seurat_obj <- qs::qread(seurat_qs_path)

message("Seurat object dimensions: ", nrow(seurat_obj), " genes x ", ncol(seurat_obj), " cells")
message("Available reductions: ", paste(Seurat::Reductions(seurat_obj), collapse = ", "))

# Check required metadata columns
required_cols <- c("batch", "Core", "intra_T", "TIL_2", "Annot3")
missing_cols <- required_cols[!required_cols %in% colnames(seurat_obj@meta.data)]
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Check available values
message("\n", "=", rep("=", 79))
message("Metadata Summary")
message("=", rep("=", 79))
message("batch values: ", paste(unique(seurat_obj@meta.data$batch), collapse = ", "))
message("Core (sample) count: ", length(unique(seurat_obj@meta.data$Core)))
message("intra_T values: ", paste(unique(seurat_obj@meta.data$intra_T), collapse = ", "))
message("TIL_2 values: ", paste(unique(seurat_obj@meta.data$TIL_2), collapse = ", "))
message("Annot3 (cluster) values: ", paste(unique(seurat_obj@meta.data$Annot3), collapse = ", "))

# Filter for low vs high comparisons
message("\n", "=", rep("=", 79))
message("Filtering data for low vs high comparisons")
message("=", rep("=", 79))

# Function to run Milo analysis for a specific target variable
run_milo_comparison <- function(target_var, comparison_name) {
  message("\n", paste(rep("=", 80), collapse = ""))
  message("Running Milo analysis for: ", comparison_name)
  message("Target variable: ", target_var)
  message(paste(rep("=", 80), collapse = ""))
  
  # Filter to low and high only (exclude mid if present)
  seurat_subset <- subset(seurat_obj, 
                         cells = which(seurat_obj@meta.data[[target_var]] %in% c("low", "high")))
  
  message("Filtered to ", ncol(seurat_subset), " cells (low + high only)")
  message("low: ", sum(seurat_subset@meta.data[[target_var]] == "low"), " cells")
  message("high: ", sum(seurat_subset@meta.data[[target_var]] == "high"), " cells")
  
  # Check if we have enough samples
  sample_counts <- table(seurat_subset@meta.data$Core, seurat_subset@meta.data[[target_var]])
  message("\nSample counts by Core and ", target_var, ":")
  print(sample_counts)
  
  # Check minimum samples per group
  min_samples_per_group <- min(apply(sample_counts, 2, function(x) sum(x > 0)))
  if (min_samples_per_group < 2) {
    warning("Warning: Some groups have less than 2 samples. Analysis may fail.")
  }
  
  # Determine which reduction to use for graph building
  # Prefer integrated.scvi if available, otherwise use pca
  available_reductions <- Seurat::Reductions(seurat_subset)
  if ("integrated.scvi" %in% available_reductions) {
    graph_reduction <- "integrated.scvi"
    layout_reduction <- "umap"  # or "umap.scvi" if available
    message("Using graph_reduction: ", graph_reduction)
    message("Using layout_reduction: ", layout_reduction)
  } else if ("pca" %in% available_reductions) {
    graph_reduction <- "pca"
    layout_reduction <- "umap"
    message("Using graph_reduction: ", graph_reduction)
    message("Using layout_reduction: ", layout_reduction)
  } else {
    stop("No suitable reduction found for graph building. Available: ", 
         paste(available_reductions, collapse = ", "))
  }
  
  # Run Milo pipeline
  output_dir <- file.path(output_base_dir, comparison_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\nStarting Milo pipeline...")
  milo_results <- run_milo_pipeline(
    seurat_obj = seurat_subset,
    patient_var = "Core",
    cluster_var = "Annot3",  # Use Annot3 as cluster variable
    target_var = target_var,
    batch_var = "batch",
    graph_reduction = graph_reduction,
    layout_reduction = layout_reduction,
    k = 30,
    d = 30,
    prop = 0.1,
    alpha = 0.1,
    output_dir = output_dir,
    prefix = paste0("milo_", comparison_name),
    verbose = TRUE,
    save = TRUE
  )
  
  message("\n", paste(rep("=", 80), collapse = ""))
  message("Milo analysis completed for ", comparison_name)
  message("Results saved to: ", output_dir)
  message(paste(rep("=", 80), collapse = ""))
  
  # Print summary of DA results
  if (!is.null(milo_results$da_results)) {
    da_results <- milo_results$da_results
    message("\nDA Results Summary:")
    message("  Total neighborhoods tested: ", nrow(da_results))
    message("  Significant (FDR < 0.1): ", sum(da_results$SpatialFDR < 0.1, na.rm = TRUE))
    message("  Significant (FDR < 0.05): ", sum(da_results$SpatialFDR < 0.05, na.rm = TRUE))
    message("  Mean logFC: ", round(mean(da_results$logFC, na.rm = TRUE), 4))
    message("  SD logFC: ", round(sd(da_results$logFC, na.rm = TRUE), 4))
  }
  
  return(milo_results)
}

# Run analysis for intra_T
message("\n", paste(rep("#", 80), collapse = ""))
message("ANALYSIS 1: intra_T low vs high")
message(paste(rep("#", 80), collapse = ""))
milo_intra_T <- run_milo_comparison("intra_T", "intra_T_low_vs_high")

# Run analysis for TIL_2
message("\n", paste(rep("#", 80), collapse = ""))
message("ANALYSIS 2: TIL_2 low vs high")
message(paste(rep("#", 80), collapse = ""))
milo_TIL_2 <- run_milo_comparison("TIL_2", "TIL_2_low_vs_high")

message("\n", paste(rep("=", 80), collapse = ""))
message("All Milo analyses completed!")
message("Results saved in: ", output_base_dir)
message(paste(rep("=", 80), collapse = ""))

