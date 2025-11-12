#' Test Script for Standardized Plot Functions
#'
#' Tests plot functions with TXNIP gene and 3-level grouping hierarchy
#' 
#' Hierarchy:
#' 1. Cell-level: cell x feature (group: anno3.scvi -> hos_no -> g3)
#' 2. Patient-level: patient x feature (split: g3)
#' 3. Group-level: group x feature (split: g3)

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Source plot functions
source("myR/R/utils_data.R")
source("myR/R/utils_aggregation.R")
source("myR/R/plots_scatter.R")
source("myR/R/plots_volcano.R")
source("myR/R/plots_heatmap.R")

# Test function
test_plots <- function(sobj, 
                       feature = "TXNIP",
                       output_dir = "test_output",
                       group.by = "anno3.scvi",
                       sample_col = "hos_no",
                       split.by = "g3") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Prepare features: gene + numeric metadata
  numeric_meta <- names(sobj@meta.data)[sapply(sobj@meta.data, is.numeric)]
  features <- c(feature, numeric_meta[1:min(3, length(numeric_meta))])  # Gene + up to 3 numeric metadata
  
  cat("Testing with features:", paste(features, collapse = ", "), "\n")
  
  # ============================================
  # 1. Cell-level: cell x feature
  # ============================================
  cat("\n1. Cell-level plots...\n")
  
  # Scatter plot (if we have a numeric covariate)
  if (length(numeric_meta) > 0) {
    tryCatch({
      p_scatter_cell <- plot_scatter(
        data = sobj,
        feature = feature,
        x_var = numeric_meta[1],
        group.by = group.by,
        aggregate = FALSE,  # Cell-level
        fitted_line = "linear"
      )
      ggsave(file.path(output_dir, "01_scatter_cell.png"), 
             p_scatter_cell, width = 10, height = 7, dpi = 300)
      cat("  ✓ Scatter plot (cell-level) saved\n")
    }, error = function(e) {
      cat("  ✗ Scatter plot (cell-level) failed:", e$message, "\n")
    })
  }
  
  # Heatmap (genes)
  tryCatch({
    p_heatmap_cell <- plot_heatmap_genes(
      data = sobj,
      features = c(feature),
      group.by = group.by,
      aggregate = FALSE,  # Cell-level
      normalize = TRUE,
      title = paste("Cell-level:", feature)
    )
    ggsave(file.path(output_dir, "02_heatmap_cell.png"), 
           p_heatmap_cell, width = 10, height = 7, dpi = 300)
    cat("  ✓ Heatmap (cell-level) saved\n")
  }, error = function(e) {
    cat("  ✗ Heatmap (cell-level) failed:", e$message, "\n")
  })
  
  # ============================================
  # 2. Patient-level: patient x feature (split: g3)
  # ============================================
  cat("\n2. Patient-level plots (aggregated by", sample_col, ")...\n")
  
  # Scatter plot
  if (length(numeric_meta) > 0) {
    tryCatch({
      p_scatter_patient <- plot_scatter(
        data = sobj,
        feature = feature,
        x_var = numeric_meta[1],
        group.by = sample_col,
        split.by = split.by,
        aggregate = TRUE,
        fitted_line = "linear"
      )
      ggsave(file.path(output_dir, "03_scatter_patient.png"), 
             p_scatter_patient, width = 12, height = 7, dpi = 300)
      cat("  ✓ Scatter plot (patient-level) saved\n")
    }, error = function(e) {
      cat("  ✗ Scatter plot (patient-level) failed:", e$message, "\n")
    })
  }
  
  # Heatmap
  tryCatch({
    p_heatmap_patient <- plot_heatmap_genes(
      data = sobj,
      features = c(feature),
      group.by = sample_col,
      split.by = split.by,
      aggregate = TRUE,
      aggregate_by = sample_col,
      normalize = TRUE,
      title = paste("Patient-level:", feature)
    )
    ggsave(file.path(output_dir, "04_heatmap_patient.png"), 
           p_heatmap_patient, width = 12, height = 7, dpi = 300)
    cat("  ✓ Heatmap (patient-level) saved\n")
  }, error = function(e) {
    cat("  ✗ Heatmap (patient-level) failed:", e$message, "\n")
  })
  
  # ============================================
  # 3. Group-level: group x feature (split: g3)
  # ============================================
  cat("\n3. Group-level plots (aggregated by", split.by, ")...\n")
  
  # Scatter plot
  if (length(numeric_meta) > 0) {
    tryCatch({
      p_scatter_group <- plot_scatter(
        data = sobj,
        feature = feature,
        x_var = numeric_meta[1],
        group.by = split.by,
        aggregate = TRUE,
        aggregate_by = split.by,
        fitted_line = "linear"
      )
      ggsave(file.path(output_dir, "05_scatter_group.png"), 
             p_scatter_group, width = 10, height = 7, dpi = 300)
      cat("  ✓ Scatter plot (group-level) saved\n")
    }, error = function(e) {
      cat("  ✗ Scatter plot (group-level) failed:", e$message, "\n")
    })
  }
  
  # Heatmap
  tryCatch({
    p_heatmap_group <- plot_heatmap_genes(
      data = sobj,
      features = c(feature),
      group.by = split.by,
      aggregate = TRUE,
      aggregate_by = split.by,
      normalize = TRUE,
      title = paste("Group-level:", feature)
    )
    ggsave(file.path(output_dir, "06_heatmap_group.png"), 
           p_heatmap_group, width = 10, height = 7, dpi = 300)
    cat("  ✓ Heatmap (group-level) saved\n")
  }, error = function(e) {
    cat("  ✗ Heatmap (group-level) failed:", e$message, "\n")
  })
  
  cat("\nAll plots saved to:", output_dir, "\n")
}

# Usage:
# test_plots(sobj, feature = "TXNIP", output_dir = "test_output")

