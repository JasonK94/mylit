#' Debug Test Script for Standardized Plot Functions
#'
#' Enhanced test script with better error handling and debugging

# Load required libraries with error checking
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "rlang", "tidyr", "viridisLite", "RColorBrewer")
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
    stop("Error loading ", f, ": ", e$message)
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

# Test function with enhanced debugging
test_plots_debug <- function(sobj, 
                             feature = c("TXNIP", "DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C"),
                             output_dir = "test_output",
                             group.by = "anno3.scvi",
                             sample_col = "hos_no",
                             split.by = "g3") {
  
  cat("=== Starting Plot Tests ===\n")
  cat("Input validation...\n")
  
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
  cat("  ✓ All required columns present\n")
  
  # Check features
  cat("Checking features...\n")
  assay <- Seurat::DefaultAssay(sobj)
  available_genes <- rownames(sobj[[assay]])
  available_meta <- names(sobj@meta.data)
  
  # Separate genes and metadata
  feature_genes <- feature[feature %in% available_genes]
  feature_meta <- feature[feature %in% available_meta]
  feature_missing <- setdiff(feature, c(available_genes, available_meta))
  
  if (length(feature_missing) > 0) {
    warning("Features not found: ", paste(feature_missing, collapse = ", "))
  }
  
  if (length(feature_genes) == 0 && length(feature_meta) == 0) {
    stop("No valid features found (genes or metadata)")
  }
  
  cat("  ✓ Genes found:", length(feature_genes), paste(feature_genes, collapse = ", "), "\n")
  cat("  ✓ Metadata found:", length(feature_meta), paste(feature_meta, collapse = ", "), "\n")
  
  # Use first gene for scatter plots, all features for heatmaps
  test_gene <- if (length(feature_genes) > 0) feature_genes[1] else feature_meta[1]
  test_features <- c(feature_genes, feature_meta)
  
  # Prepare numeric metadata for scatter plots
  numeric_meta <- names(sobj@meta.data)[sapply(sobj@meta.data, is.numeric)]
  numeric_meta <- setdiff(numeric_meta, test_features)  # Exclude test features
  
  if (length(numeric_meta) == 0) {
    warning("No numeric metadata found for scatter plot x_var")
  }
  
  cat("  ✓ Numeric metadata for x_var:", length(numeric_meta), "\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  ✓ Created output directory:", output_dir, "\n")
  }
  
  cat("\n=== Test 1: Cell-level plots ===\n")
  
  # 1. Cell-level scatter plot
  if (length(numeric_meta) > 0 && length(feature_genes) > 0) {
    cat("  Testing scatter plot (cell-level)...\n")
    tryCatch({
      p_scatter_cell <- plot_scatter(
        data = sobj,
        feature = test_gene,
        x_var = numeric_meta[1],
        group.by = group.by,
        aggregate = FALSE,
        fitted_line = "linear"
      )
      ggsave(file.path(output_dir, "01_scatter_cell.png"), 
             p_scatter_cell, width = 10, height = 7, dpi = 300)
      cat("    ✓ Saved: 01_scatter_cell.png\n")
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
      cat("    Stack trace:\n")
      print(traceback())
    })
  } else {
    cat("    ⊘ Skipped (no numeric metadata or genes)\n")
  }
  
  # 2. Cell-level heatmap
  cat("  Testing heatmap (cell-level)...\n")
  tryCatch({
    p_heatmap_cell <- plot_heatmap_genes(
      data = sobj,
      features = test_features[1:min(5, length(test_features))],  # Limit to 5 features
      group.by = group.by,
      aggregate = FALSE,
      normalize = TRUE,
      title = paste("Cell-level:", paste(test_features[1:min(3, length(test_features))], collapse = ", "))
    )
    ggsave(file.path(output_dir, "02_heatmap_cell.png"), 
           p_heatmap_cell, width = 10, height = 7, dpi = 300)
    cat("    ✓ Saved: 02_heatmap_cell.png\n")
  }, error = function(e) {
    cat("    ✗ Error:", e$message, "\n")
    cat("    Stack trace:\n")
    print(traceback())
  })
  
  cat("\n=== Test 2: Patient-level plots ===\n")
  
  # 3. Patient-level scatter plot
  if (length(numeric_meta) > 0 && length(feature_genes) > 0) {
    cat("  Testing scatter plot (patient-level)...\n")
    tryCatch({
      p_scatter_patient <- plot_scatter(
        data = sobj,
        feature = test_gene,
        x_var = numeric_meta[1],
        group.by = sample_col,
        split.by = split.by,
        aggregate = TRUE,
        fitted_line = "linear"
      )
      ggsave(file.path(output_dir, "03_scatter_patient.png"), 
             p_scatter_patient, width = 12, height = 7, dpi = 300)
      cat("    ✓ Saved: 03_scatter_patient.png\n")
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
      cat("    Stack trace:\n")
      print(traceback())
    })
  } else {
    cat("    ⊘ Skipped (no numeric metadata or genes)\n")
  }
  
  # 4. Patient-level heatmap
  cat("  Testing heatmap (patient-level)...\n")
  tryCatch({
    p_heatmap_patient <- plot_heatmap_genes(
      data = sobj,
      features = test_features[1:min(5, length(test_features))],
      group.by = sample_col,
      split.by = split.by,
      aggregate = TRUE,
      aggregate_by = sample_col,
      normalize = TRUE,
      title = paste("Patient-level:", paste(test_features[1:min(3, length(test_features))], collapse = ", "))
    )
    ggsave(file.path(output_dir, "04_heatmap_patient.png"), 
           p_heatmap_patient, width = 12, height = 7, dpi = 300)
    cat("    ✓ Saved: 04_heatmap_patient.png\n")
  }, error = function(e) {
    cat("    ✗ Error:", e$message, "\n")
    cat("    Stack trace:\n")
    print(traceback())
  })
  
  cat("\n=== Test 3: Group-level plots ===\n")
  
  # 5. Group-level scatter plot
  if (length(numeric_meta) > 0 && length(feature_genes) > 0) {
    cat("  Testing scatter plot (group-level)...\n")
    tryCatch({
      p_scatter_group <- plot_scatter(
        data = sobj,
        feature = test_gene,
        x_var = numeric_meta[1],
        group.by = split.by,
        aggregate = TRUE,
        aggregate_by = split.by,
        fitted_line = "linear"
      )
      ggsave(file.path(output_dir, "05_scatter_group.png"), 
             p_scatter_group, width = 10, height = 7, dpi = 300)
      cat("    ✓ Saved: 05_scatter_group.png\n")
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
      cat("    Stack trace:\n")
      print(traceback())
    })
  } else {
    cat("    ⊘ Skipped (no numeric metadata or genes)\n")
  }
  
  # 6. Group-level heatmap
  cat("  Testing heatmap (group-level)...\n")
  tryCatch({
    p_heatmap_group <- plot_heatmap_genes(
      data = sobj,
      features = test_features[1:min(5, length(test_features))],
      group.by = split.by,
      aggregate = TRUE,
      aggregate_by = split.by,
      normalize = TRUE,
      title = paste("Group-level:", paste(test_features[1:min(3, length(test_features))], collapse = ", "))
    )
    ggsave(file.path(output_dir, "06_heatmap_group.png"), 
           p_heatmap_group, width = 10, height = 7, dpi = 300)
    cat("    ✓ Saved: 06_heatmap_group.png\n")
  }, error = function(e) {
    cat("    ✗ Error:", e$message, "\n")
    cat("    Stack trace:\n")
    print(traceback())
  })
  
  cat("\n=== Tests Complete ===\n")
  cat("Output directory:", output_dir, "\n")
  cat("Generated files:\n")
  files <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
  if (length(files) > 0) {
    for (f in files) {
      cat("  -", f, "\n")
    }
  } else {
    cat("  (no files generated)\n")
  }
}

# Quick test function
quick_test <- function(sobj) {
  cat("Quick function availability test...\n")
  
  # Test .sobj_to_df
  cat("Testing .sobj_to_df...\n")
  tryCatch({
    df <- .sobj_to_df(sobj, features = NULL, metadata_only = TRUE)
    cat("  ✓ .sobj_to_df works, returned", nrow(df), "rows\n")
  }, error = function(e) {
    cat("  ✗ .sobj_to_df failed:", e$message, "\n")
  })
  
  # Test .prepare_plot_data
  cat("Testing .prepare_plot_data...\n")
  tryCatch({
    plot_df <- .prepare_plot_data(sobj, features = NULL, group.by = "anno3.scvi")
    cat("  ✓ .prepare_plot_data works, returned", nrow(plot_df), "rows\n")
  }, error = function(e) {
    cat("  ✗ .prepare_plot_data failed:", e$message, "\n")
  })
  
  cat("\n")
}

