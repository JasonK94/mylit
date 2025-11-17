#!/usr/bin/env Rscript

# ============================================================================
# Test script: run DEG consensus on FULL dataset inside start.R / renv environment
# ============================================================================

cat("=== Initialising start.R environment (FULL data) ===\n")

# 1. Move into the st project so that start.R detects the correct repo root
setwd("/home/user3/GJC_KDW_250721/st")

# 2. Source start.R (loads renv-like environment, core packages, and myR)
source("../start.R")

cat("=== Loading Seurat object (FULL) ===\n")

if (!requireNamespace("qs", quietly = TRUE)) {
  stop("Package 'qs' is required but not installed in this environment.")
}

is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

cat(sprintf("Loaded Seurat object (FULL): %d cells, %d genes\n", ncol(is5), nrow(is5)))

cat("=== Running DEG consensus pipeline (run_consensus_simple.R) on FULL data ===\n")

source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/run_consensus_simple.R")

cat("=== test_full_in_start_env.R completed ===\n")


