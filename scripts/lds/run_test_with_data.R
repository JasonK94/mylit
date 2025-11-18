#!/usr/bin/env Rscript
# Run tests with provided Seurat objects
# Usage: Rscript run_test_with_data.R [downsampled|full]

cat("=== Plot Function Test with Provided Data ===\n\n")

# Check and load required packages first
cat("Checking required packages...\n")
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "tidyr", "viridisLite", "RColorBrewer", "qs")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), 
       "\nPlease install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
}

# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(qs)
cat("All packages loaded successfully!\n\n")

# Source test script
test_script <- "/data/user3/git_repo/_wt/plots/test_plots.R"
if (!file.exists(test_script)) {
  stop("Test script not found: ", test_script)
}

# Source without loading packages again (they're already loaded)
# Read the file and source only the test function
source(test_script, local = TRUE)

# Determine which file to use
args <- commandArgs(trailingOnly = TRUE)
use_full <- if (length(args) > 0 && args[1] == "full") TRUE else FALSE

if (use_full) {
  sobj_path <- "/data/user3/sobj/IS5_g3NA_removal_251110.qs"
  output_dir <- "test_output_full"
  cat("Using FULL dataset\n")
} else {
  sobj_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
  output_dir <- "test_output_ds2500"
  cat("Using DOWNSAMPLED dataset (2500 cells)\n")
}

cat("Loading Seurat object from:", sobj_path, "\n")

# Check if file exists
if (!file.exists(sobj_path)) {
  stop("File not found: ", sobj_path)
}

# Load qs file
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs package required. Install with: install.packages('qs')")
}

cat("Loading...\n")
sobj <- qs::qread(sobj_path)

# Validate sobj
cat("\nSeurat object validation:\n")
cat("  Cells:", ncol(sobj), "\n")
cat("  Features:", nrow(sobj), "\n")
cat("  Default assay:", Seurat::DefaultAssay(sobj), "\n")
cat("  Assays:", paste(names(sobj@assays), collapse = ", "), "\n")
cat("  Metadata columns:", length(names(sobj@meta.data)), "\n")
cat("  Sample metadata columns:\n")
meta_cols <- names(sobj@meta.data)
cat("    ", paste(head(meta_cols, 20), collapse = ", "), "\n")
if (length(meta_cols) > 20) {
  cat("    ... and", length(meta_cols) - 20, "more\n")
}

# Check required columns
group.by <- "anno3.scvi"
sample_col <- "hos_no"
split.by <- "g3"

required_cols <- c(group.by, sample_col, split.by)
available_cols <- names(sobj@meta.data)
missing_cols <- setdiff(required_cols, available_cols)

if (length(missing_cols) > 0) {
  cat("\nWARNING: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns:\n")
  cat(paste(head(available_cols, 30), collapse = ", "), "\n\n")
  cat("Trying to find alternatives...\n")
  
  # Try to find alternatives
  if (!group.by %in% available_cols) {
    cluster_cols <- grep("cluster|anno|cell_type|ident|scvi", available_cols, ignore.case = TRUE, value = TRUE)
    if (length(cluster_cols) > 0) {
      group.by <- cluster_cols[1]
      cat("  Using", group.by, "for group.by\n")
    }
  }
  
  if (!sample_col %in% available_cols) {
    sample_cols <- grep("sample|patient|hos|id", available_cols, ignore.case = TRUE, value = TRUE)
    if (length(sample_cols) > 0) {
      sample_col <- sample_cols[1]
      cat("  Using", sample_col, "for sample_col\n")
    }
  }
  
  if (!split.by %in% available_cols) {
    group_cols <- grep("group|condition|treatment|g[0-9]|g3", available_cols, ignore.case = TRUE, value = TRUE)
    if (length(group_cols) > 0) {
      split.by <- group_cols[1]
      cat("  Using", split.by, "for split.by\n")
    }
  }
  
  # Re-check
  required_cols <- c(group.by, sample_col, split.by)
  missing_cols <- setdiff(required_cols, available_cols)
  if (length(missing_cols) > 0) {
    stop("Cannot proceed without columns: ", paste(missing_cols, collapse = ", "))
  }
} else {
  cat("\n✓ All required columns found\n")
}

# Check features
assay <- Seurat::DefaultAssay(sobj)
available_genes <- rownames(sobj[[assay]])
test_features <- c("TXNIP", "DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C")
available_test_features <- intersect(test_features, available_genes)

cat("\nFeature check:\n")
cat("  Test features requested:", length(test_features), "\n")
cat("  Available in assay:", length(available_test_features), "\n")
if (length(available_test_features) > 0) {
  cat("  Found:", paste(available_test_features, collapse = ", "), "\n")
  test_features <- available_test_features
} else {
  cat("  WARNING: None of the test features found in assay.\n")
  cat("  Using first 5 genes as test features.\n")
  test_features <- head(available_genes, 5)
  cat("  Using:", paste(test_features, collapse = ", "), "\n")
}

# Check numeric metadata for scatter plots
numeric_meta <- names(sobj@meta.data)[sapply(sobj@meta.data, is.numeric)]
cat("  Numeric metadata columns:", length(numeric_meta), "\n")
if (length(numeric_meta) > 0) {
  cat("    Examples:", paste(head(numeric_meta, 5), collapse = ", "), "\n")
}

# Run tests
cat("\n=== Running Tests ===\n")
cat("Parameters:\n")
cat("  group.by:", group.by, "\n")
cat("  sample_col:", sample_col, "\n")
cat("  split.by:", split.by, "\n")
cat("  features:", paste(test_features, collapse = ", "), "\n")
cat("  output_dir:", output_dir, "\n\n")

start_time <- Sys.time()

tryCatch({
  test_plots(
    sobj = sobj,
    feature = test_features,
    output_dir = output_dir,
    group.by = group.by,
    sample_col = sample_col,
    split.by = split.by
  )
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat("\n✓ All tests completed successfully!\n")
  cat("  Time elapsed:", round(elapsed, 2), "seconds\n")
  cat("  Output directory:", output_dir, "\n")
  
  # List generated files
  if (dir.exists(output_dir)) {
    files <- list.files(output_dir, pattern = "\\.png$", full.names = FALSE)
    if (length(files) > 0) {
      cat("  Generated files:\n")
      for (f in files) {
        file_path <- file.path(output_dir, f)
        file_size <- file.info(file_path)$size
        cat("    -", f, paste0("(", round(file_size/1024, 1), " KB)\n"))
      }
    }
  }
  
}, error = function(e) {
  cat("\n✗ Test execution failed:\n")
  cat(e$message, "\n")
  cat("\nStack trace:\n")
  print(traceback())
  quit(status = 1)
})

cat("\n=== Test Complete ===\n")

