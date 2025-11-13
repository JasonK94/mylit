#' Test Script for Standardized Plot Functions
#'
#' Tests plot functions with TXNIP gene and 3-level grouping hierarchy
#' 
#' Hierarchy:
#' 1. Cell-level: cell x feature (group: anno3.scvi -> hos_no -> g3)
#' 2. Patient-level: patient x feature (split: g3)
#' 3. Group-level: group x feature (split: g3)

# Load required libraries with error checking
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "tidyr", "viridisLite", "RColorBrewer")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Source plot functions with error checking
source_files <- c(
  "/data/user3/git_repo/_wt/plots/myR/R/utils_data.R",
  "/data/user3/git_repo/_wt/plots/myR/R/utils_aggregation.R",
  "/data/user3/git_repo/_wt/plots/myR/R/plots_scatter.R",
  "/data/user3/git_repo/_wt/plots/myR/R/plots_volcano.R",
  "/data/user3/git_repo/_wt/plots/myR/R/plots_heatmap.R"
)

cat("Sourcing plot functions...\n")
for (f in source_files) {
  if (!file.exists(f)) {
    stop("File not found: ", f)
  }
  tryCatch({
    source(f)
    cat("  ✓ Loaded:", basename(f), "\n")
  }, error = function(e) {
    stop("Error loading ", f, ": ", e$message, "\nCall stack:\n", paste(capture.output(traceback()), collapse = "\n"))
  })
}

# Verify functions are loaded
required_functions <- c(".sobj_to_df", ".validate_df", ".prepare_plot_data", 
                       ".aggregate_cells", ".prepare_data_with_aggregation",
                       "plot_scatter", "plot_volcano", "plot_heatmap_genes")
missing_functions <- required_functions[!sapply(required_functions, exists)]

if (length(missing_functions) > 0) {
  stop("Missing required functions: ", paste(missing_functions, collapse = ", "))
}
cat("All functions loaded successfully!\n\n")

# Test function
test_plots <- function(sobj, 
                       feature = c("nih_change", "DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C", "TXNIP"),
                       output_dir = "test_output",
                       group.by = "anno3.scvi",
                       sample_col = "hos_no",
                       split.by = "g3") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Validate input
  if (!inherits(sobj, "Seurat")) {
    stop("sobj must be a Seurat object")
  }
  
  # Check required columns
  required_cols <- c(group.by, sample_col, split.by)
  missing_cols <- setdiff(required_cols, names(sobj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check features - separate genes and metadata
  assay <- Seurat::DefaultAssay(sobj)
  available_genes <- rownames(sobj[[assay]])
  available_meta <- names(sobj@meta.data)
  
  feature_genes <- feature[feature %in% available_genes]
  feature_meta <- feature[feature %in% available_meta]
  feature_missing <- setdiff(feature, c(available_genes, available_meta))
  
  if (length(feature_missing) > 0) {
    warning("Features not found: ", paste(feature_missing, collapse = ", "))
  }
  
  if (length(feature_genes) == 0 && length(feature_meta) == 0) {
    stop("No valid features found (genes or metadata)")
  }
  
  # Use first gene for scatter plots
  test_gene <- if (length(feature_genes) > 0) feature_genes[1] else feature_meta[1]
  test_features <- c(feature_genes, feature_meta)
  
  # Prepare numeric metadata for scatter plots (exclude test features)
  numeric_meta <- names(sobj@meta.data)[sapply(sobj@meta.data, is.numeric)]
  numeric_meta <- setdiff(numeric_meta, test_features)
  
  cat("Testing with features:\n")
  cat("  Genes:", length(feature_genes), paste(feature_genes, collapse = ", "), "\n")
  cat("  Metadata:", length(feature_meta), paste(feature_meta, collapse = ", "), "\n")
  cat("  Numeric metadata for x_var:", length(numeric_meta), "\n")
  if (length(numeric_meta) > 0) {
    cat("    Using:", numeric_meta[1], "\n")
  }
  
  # ============================================
  # 1. Cluster-level: cluster x feature (aggregated by cluster)
  # ============================================
  cat("\n1. Cluster-level plots (aggregated by", group.by, ")...\n")
  
  # Scatter plot (if we have a numeric covariate)
  if (length(numeric_meta) > 0) {
    tryCatch({
      p_scatter_cluster <- plot_scatter(
        data = sobj,
        feature = test_gene,
        x_var = numeric_meta[1],
        group.by = group.by,
        color_by = split.by,  # Color by g3 to see groups
        aggregate = TRUE,  # Aggregate by cluster
        aggregate_by = group.by,
        fitted_line = "linear",
        each_fit = TRUE,  # Fit separately for each g3 group
        show_stats = TRUE
      )
      ggsave(file.path(output_dir, "01_scatter_cluster.png"), 
             p_scatter_cluster, width = 10, height = 7, dpi = 300)
      cat("  ✓ Scatter plot (cluster-level) saved\n")
    }, error = function(e) {
      cat("  ✗ Scatter plot (cluster-level) failed:", e$message, "\n")
      if (getOption("show.error.messages", TRUE)) {
        cat("  Stack trace:\n")
        print(traceback())
      }
    })
  }
  
  # Heatmap (genes) - normalized at gene level
  tryCatch({
    p_heatmap_cluster <- plot_heatmap_genes(
      data = sobj,
      features = test_features[1:min(5, length(test_features))],
      group.by = group.by,
      aggregate = TRUE,
      aggregate_by = group.by,
      normalize = TRUE,
      normalize_by = "row",  # Normalize by gene (row)
      normalize_transpose = FALSE,  # FALSE = column normalization (sample/group level)
      show_group_separator = TRUE,  # Add vertical lines between groups
      title = paste("Cluster-level (column nmz):", paste(test_features[1:min(3, length(test_features))], collapse = ", "))
    )
    ggsave(file.path(output_dir, "02_heatmap_cluster.png"), 
           p_heatmap_cluster, width = 10, height = 7, dpi = 300)
    cat("  ✓ Heatmap (cluster-level) saved\n")
  }, error = function(e) {
    cat("  ✗ Heatmap (cluster-level) failed:", e$message, "\n")
    if (getOption("show.error.messages", TRUE)) {
      cat("  Stack trace:\n")
      print(traceback())
    }
  })
  
  # ============================================
  # 2. Patient-level: patient x feature (split by g3)
  # ============================================
  cat("\n2. Patient-level plots (aggregated by", sample_col, ", split by", split.by, ")...\n")
  
  # Scatter plot
  if (length(numeric_meta) > 0) {
    tryCatch({
      p_scatter_patient <- plot_scatter(
        data = sobj,
        feature = test_gene,
        x_var = numeric_meta[1],
        group.by = sample_col,
        color_by = split.by,  # Color by g3
        aggregate = TRUE,
        aggregate_by = sample_col,  # Aggregate by patient only
        fitted_line = "linear",
        each_fit = TRUE,  # Fit separately for each g3 group
        show_stats = TRUE
      )
      ggsave(file.path(output_dir, "03_scatter_patient.png"), 
             p_scatter_patient, width = 12, height = 7, dpi = 300)
      cat("  ✓ Scatter plot (patient-level) saved\n")
    }, error = function(e) {
      cat("  ✗ Scatter plot (patient-level) failed:", e$message, "\n")
      if (getOption("show.error.messages", TRUE)) {
        cat("  Stack trace:\n")
        print(traceback())
      }
    })
  }
  
  # Heatmap - split by g3 (g3==1 left, g3==2 right)
  tryCatch({
    p_heatmap_patient <- plot_heatmap_genes(
      data = sobj,
      features = test_features[1:min(5, length(test_features))],
      group.by = sample_col,
      split.by = split.by,  # Split by g3 for left/right panels
      aggregate = TRUE,
      aggregate_by = sample_col,  # Aggregate by patient
      normalize = TRUE,
      normalize_by = "row",  # Normalize by gene (row)
      normalize_transpose = FALSE,  # FALSE = column normalization (sample/group level)
      show_group_separator = TRUE,  # Add vertical lines between groups (including between g3=1 and g3=2)
      title = paste("Patient-level (split by", split.by, ", column nmz):", paste(test_features[1:min(3, length(test_features))], collapse = ", "))
    )
    ggsave(file.path(output_dir, "04_heatmap_patient.png"), 
           p_heatmap_patient, width = 14, height = 7, dpi = 300)
    cat("  ✓ Heatmap (patient-level, split by", split.by, ") saved\n")
  }, error = function(e) {
    cat("  ✗ Heatmap (patient-level) failed:", e$message, "\n")
    if (getOption("show.error.messages", TRUE)) {
      cat("  Stack trace:\n")
      print(traceback())
    }
  })
  
  # ============================================
  # 3. Group-level: group x feature (aggregated by g3)
  # ============================================
  cat("\n3. Group-level plots (aggregated by", split.by, ")...\n")
  
  # Scatter plot
  if (length(numeric_meta) > 0) {
    tryCatch({
      # Check if we have enough data points
      group_values <- unique(sobj@meta.data[[split.by]])
      group_values <- group_values[!is.na(group_values)]
      
      if (length(group_values) >= 2) {
        p_scatter_group <- plot_scatter(
          data = sobj,
          feature = test_gene,
          x_var = numeric_meta[1],
          group.by = split.by,
          aggregate = TRUE,
          aggregate_by = split.by,
          fitted_line = "linear",
          label = TRUE,  # Default label for group-level
          color_by = split.by,  # Color by group for legend
          each_fit = TRUE,  # Show each_fit option
          show_stats = TRUE,
          remove_na = FALSE  # Keep NA for categorical treatment
        )
        ggsave(file.path(output_dir, "05_scatter_group.png"), 
               p_scatter_group, width = 10, height = 7, dpi = 300)
        cat("  ✓ Scatter plot (group-level) saved\n")
      } else {
        cat("  ⊘ Scatter plot (group-level) skipped: only", length(group_values), "group(s) available\n")
      }
    }, error = function(e) {
      cat("  ✗ Scatter plot (group-level) failed:", e$message, "\n")
      if (getOption("show.error.messages", TRUE)) {
        cat("  Stack trace:\n")
        print(traceback())
      }
    })
  }
  
  # Heatmap
  tryCatch({
    p_heatmap_group <- plot_heatmap_genes(
      data = sobj,
      features = test_features[1:min(5, length(test_features))],
      group.by = split.by,
      aggregate = TRUE,
      aggregate_by = split.by,
      normalize = TRUE,
      normalize_by = "row",  # Normalize by gene (row) - default
      normalize_transpose = FALSE,  # Keep row normalization
      remove_na = TRUE,  # Remove NA values
      show_group_separator = TRUE,  # Add vertical lines between groups
      title = paste("Group-level:", paste(test_features[1:min(3, length(test_features))], collapse = ", "))
    )
    ggsave(file.path(output_dir, "06_heatmap_group.png"), 
           p_heatmap_group, width = 10, height = 7, dpi = 300)
    cat("  ✓ Heatmap (group-level) saved\n")
  }, error = function(e) {
    cat("  ✗ Heatmap (group-level) failed:", e$message, "\n")
    if (getOption("show.error.messages", TRUE)) {
      cat("  Stack trace:\n")
      print(traceback())
    }
  })
  
  cat("\nAll plots saved to:", output_dir, "\n")
}

# Usage:
# test_plots(sobj, feature = "TXNIP", output_dir = "test_output")

