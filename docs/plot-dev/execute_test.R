#!/usr/bin/env Rscript
# Execute Plot Function Tests
# Usage: Rscript execute_test.R [path_to_sobj.rds] OR set sobj in environment

cat("=== Plot Function Test Execution ===\n\n")

# Source test script
test_script <- "/data/user3/git_repo/_wt/plots/test_plots.R"
if (!file.exists(test_script)) {
  stop("Test script not found: ", test_script)
}

source(test_script)

# Check for sobj in command line args or environment
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  # Load from file
  sobj_path <- args[1]
  cat("Loading Seurat object from:", sobj_path, "\n")
  if (!file.exists(sobj_path)) {
    stop("File not found: ", sobj_path)
  }
  
  if (grepl("\\.qs$", sobj_path, ignore.case = TRUE)) {
    # qs format
    if (!requireNamespace("qs", quietly = TRUE)) {
      stop("qs package required for .qs files. Install with: install.packages('qs')")
    }
    sobj <- qs::qread(sobj_path)
  } else if (grepl("\\.rds$", sobj_path, ignore.case = TRUE)) {
    sobj <- readRDS(sobj_path)
  } else if (grepl("\\.RData$", sobj_path, ignore.case = TRUE)) {
    load(sobj_path)
    # Find Seurat object in environment
    env_objects <- ls()
    sobj_candidates <- env_objects[sapply(env_objects, function(x) {
      tryCatch(inherits(get(x), "Seurat"), error = function(e) FALSE)
    })]
    if (length(sobj_candidates) > 0) {
      sobj <- get(sobj_candidates[1])
      cat("Loaded:", sobj_candidates[1], "\n")
    } else {
      stop("No Seurat object found in RData file")
    }
  } else {
    stop("Unsupported file format. Use .rds or .RData")
  }
} else {
  # Check environment
  if (exists("sobj") && inherits(sobj, "Seurat")) {
    cat("Using sobj from environment\n")
  } else {
    stop("No Seurat object provided.\n",
         "Usage: Rscript execute_test.R [path_to_sobj.rds]\n",
         "Or set sobj variable in R environment before running.")
  }
}

# Validate sobj
cat("\nSeurat object validation:\n")
cat("  Cells:", ncol(sobj), "\n")
cat("  Features:", nrow(sobj), "\n")
cat("  Default assay:", Seurat::DefaultAssay(sobj), "\n")
cat("  Metadata columns:", length(names(sobj@meta.data)), "\n\n")

# Check required columns
group.by <- "anno3.scvi"
sample_col <- "hos_no"
split.by <- "g3"

required_cols <- c(group.by, sample_col, split.by)
available_cols <- names(sobj@meta.data)
missing_cols <- setdiff(required_cols, available_cols)

if (length(missing_cols) > 0) {
  cat("WARNING: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns (first 30):\n")
  cat(paste(head(available_cols, 30), collapse = ", "), "\n\n")
  cat("Please adjust parameters or add missing columns.\n")
  cat("Trying with available columns...\n\n")
  
  # Try to find alternatives
  if (!group.by %in% available_cols) {
    # Look for cluster-like columns
    cluster_cols <- grep("cluster|anno|cell_type|ident", available_cols, ignore.case = TRUE, value = TRUE)
    if (length(cluster_cols) > 0) {
      group.by <- cluster_cols[1]
      cat("Using", group.by, "for group.by\n")
    }
  }
  
  if (!sample_col %in% available_cols) {
    # Look for sample-like columns
    sample_cols <- grep("sample|patient|hos|id", available_cols, ignore.case = TRUE, value = TRUE)
    if (length(sample_cols) > 0) {
      sample_col <- sample_cols[1]
      cat("Using", sample_col, "for sample_col\n")
    }
  }
  
  if (!split.by %in% available_cols) {
    # Look for group-like columns
    group_cols <- grep("group|condition|treatment|g[0-9]", available_cols, ignore.case = TRUE, value = TRUE)
    if (length(group_cols) > 0) {
      split.by <- group_cols[1]
      cat("Using", split.by, "for split.by\n")
    }
  }
  
  # Re-check
  required_cols <- c(group.by, sample_col, split.by)
  missing_cols <- setdiff(required_cols, available_cols)
  if (length(missing_cols) > 0) {
    stop("Cannot proceed without columns: ", paste(missing_cols, collapse = ", "))
  }
}

# Check features
assay <- Seurat::DefaultAssay(sobj)
available_genes <- rownames(sobj[[assay]])
test_features <- c("TXNIP", "DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C")
available_test_features <- intersect(test_features, available_genes)

if (length(available_test_features) == 0) {
  cat("WARNING: None of the test features found in assay.\n")
  cat("Using first 5 genes as test features.\n")
  test_features <- head(available_genes, 5)
} else {
  test_features <- available_test_features
  cat("Test features found:", length(test_features), "\n")
}

# Run tests
cat("\n=== Running Tests ===\n")
cat("Parameters:\n")
cat("  group.by:", group.by, "\n")
cat("  sample_col:", sample_col, "\n")
cat("  split.by:", split.by, "\n")
cat("  features:", paste(test_features, collapse = ", "), "\n\n")

tryCatch({
  test_plots(
    sobj = sobj,
    feature = test_features,
    output_dir = "test_output",
    group.by = group.by,
    sample_col = sample_col,
    split.by = split.by
  )
  cat("\n✓ All tests completed successfully!\n")
}, error = function(e) {
  cat("\n✗ Test execution failed:\n")
  cat(e$message, "\n")
  cat("\nStack trace:\n")
  print(traceback())
  quit(status = 1)
})

