#!/usr/bin/env Rscript
#' Run Plot Function Tests
#'
#' This script loads test data and runs all plot function tests

# Set up error handling
options(error = function() {
  cat("\n=== ERROR OCCURRED ===\n")
  print(traceback())
  quit(status = 1)
})

# Load required libraries
cat("Loading required packages...\n")
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "tidyr", "viridisLite", "RColorBrewer")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Source plot functions
cat("\nSourcing plot functions...\n")
source_files <- c(
  "/data/user3/git_repo/_wt/plots/myR/R/utils_data.R",
  "/data/user3/git_repo/_wt/plots/myR/R/utils_aggregation.R",
  "/data/user3/git_repo/_wt/plots/myR/R/plots_scatter.R",
  "/data/user3/git_repo/_wt/plots/myR/R/plots_volcano.R",
  "/data/user3/git_repo/_wt/plots/myR/R/plots_heatmap.R"
)

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

# Source test function
source("/data/user3/git_repo/_wt/plots/test_plots.R")

# Try to find Seurat object
cat("\nLooking for Seurat object...\n")
sobj_paths <- c(
  "sobj.RData",
  "sobj.rds",
  "../sobj.RData",
  "../sobj.rds",
  "../../sobj.RData",
  "../../sobj.rds"
)

sobj <- NULL
for (path in sobj_paths) {
  if (file.exists(path)) {
    cat("Found Seurat object at:", path, "\n")
    if (grepl("\\.RData$", path)) {
      load(path)
      # Try to find Seurat object in environment
      env_objects <- ls()
      sobj_candidates <- env_objects[sapply(env_objects, function(x) inherits(get(x), "Seurat"))]
      if (length(sobj_candidates) > 0) {
        sobj <- get(sobj_candidates[1])
        cat("Loaded Seurat object:", sobj_candidates[1], "\n")
        break
      }
    } else if (grepl("\\.rds$", path)) {
      sobj <- readRDS(path)
      cat("Loaded Seurat object from RDS\n")
      break
    }
  }
}

if (is.null(sobj)) {
  cat("WARNING: No Seurat object found. Please provide one.\n")
  cat("Usage: Rscript run_tests.R [path_to_sobj.rds]\n")
  cat("Or set sobj variable in environment before running.\n\n")
  
  # Check if sobj exists in environment
  if (exists("sobj") && inherits(sobj, "Seurat")) {
    cat("Using sobj from environment...\n")
  } else {
    stop("No Seurat object available. Cannot run tests.")
  }
}

# Run tests
cat("\n=== Running Plot Tests ===\n")
cat("Seurat object info:\n")
cat("  Cells:", ncol(sobj), "\n")
cat("  Features:", nrow(sobj), "\n")
cat("  Assays:", paste(names(sobj@assays), collapse = ", "), "\n")
cat("  Default assay:", DefaultAssay(sobj), "\n")
cat("  Metadata columns:", length(names(sobj@meta.data)), "\n")
cat("  Sample columns:", paste(head(names(sobj@meta.data), 10), collapse = ", "), "\n\n")

# Check for required columns
group.by <- "anno3.scvi"
sample_col <- "hos_no"
split.by <- "g3"

missing_cols <- setdiff(c(group.by, sample_col, split.by), names(sobj@meta.data))
if (length(missing_cols) > 0) {
  cat("WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:", paste(head(names(sobj@meta.data), 20), collapse = ", "), "\n")
  cat("Please adjust group.by, sample_col, split.by in test_plots() call\n\n")
}

# Run test
test_plots(
  sobj = sobj,
  feature = c("TXNIP", "DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C"),
  output_dir = "test_output",
  group.by = group.by,
  sample_col = sample_col,
  split.by = split.by
)

cat("\n=== Tests Complete ===\n")

