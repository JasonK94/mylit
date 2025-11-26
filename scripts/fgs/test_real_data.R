#!/usr/bin/env Rscript
# ==============================================================================
# Real Data Testing Script for FGS/TML Framework (Refactored)
# ==============================================================================
# Purpose: Validate recently patched features on actual single-cell data
#          Uses centralized config and smart downsampling.
# Author: Automated testing script
# Date: 2025-11-25
# ==============================================================================

cat("=================================================================\n")
cat("Real Data Testing Script - FGS/TML Framework\n")
cat("=================================================================\n")
cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Setup ---
cat(">>> Step 1: Environment Setup\n")
source("scripts/fgs/init_fgs_env.R")
source("/home/user3/GJC_KDW_250721/start.R")

library(Seurat)
library(qs)
library(dplyr)

# Load Config & Utils
# Note: init_fgs_env.R changes WD to GJC_KDW_250721 but saves 'original_wd'
script_dir <- file.path(original_wd, "scripts/fgs")
source(file.path(original_wd, "scripts/vars_config.R"))
source(file.path(script_dir, "test_utils.R"))

# [FIX] Explicitly set conflicts preferences
if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(dplyr::select)
    conflicted::conflicts_prefer(dplyr::summarise)
    conflicted::conflicts_prefer(dplyr::filter)
    conflicted::conflicts_prefer(dplyr::mutate)
    conflicted::conflicts_prefer(dplyr::arrange)
    conflicted::conflicts_prefer(base::intersect)
    conflicted::conflicts_prefer(base::setdiff)
    conflicted::conflicts_prefer(base::union)
    conflicted::conflicts_prefer(base::colMeans)
    conflicted::conflicts_prefer(base::colSums)
    conflicted::conflicts_prefer(base::rowMeans)
    conflicted::conflicts_prefer(base::rowSums)
}

# Load FGS/TML functions
cat("Loading FGS/TML functions...\n")
# Try to load from refactored location first
sig_path <- file.path(original_wd, "myR/R/signature.R")
if (file.exists(sig_path)) {
    source(sig_path)
} else {
    warning("signature.R not found at: ", sig_path)
}
source(file.path(original_wd, "myR/R/tml_utils.R")) # Ensure utilities are loaded
cat("✓ Functions loaded\n\n")

# --- Load Real Data ---
cat("\n>>> Step 2: Loading & Preparing Data\n")

# Get configuration for Stroke dataset
config <- get_dataset_config("stroke", level = "full")
cat("Dataset:", config$name, "\n")
cat("Source Path:", config$path, "\n")

# Load and smart subset
# max_cells_per_group=200 ensures we have enough cells for training
# but keeps runtime reasonable. Rare clusters are preserved.
sobj <- prepare_test_data(config, downsample = TRUE, max_cells_per_group = 200)

cat("Final dataset size:", ncol(sobj), "cells\n")

# Set variables from config
target_var <- config$target_var
group_var <- config$patient_var # hos_no

cat("Target variable:", target_var, "\n")
print(table(sobj@meta.data[[target_var]]))

cat("Group variable (Patient):", group_var, "\n")
if (group_var %in% colnames(sobj@meta.data)) {
    cat("Number of patients:", length(unique(sobj@meta.data[[group_var]])), "\n")
} else {
    warning(sprintf("Group variable '%s' not found in metadata!", group_var))
    # Fallback to orig.ident only if absolutely necessary (though user advised against it)
    if ("orig.ident" %in% colnames(sobj@meta.data)) {
        group_var <- "orig.ident"
        cat("Fallback to 'orig.ident' (Warning: Not recommended)\n")
    }
}

# ==============================================================================
# TEST 1: FGS with All Methods Including NMF
# ==============================================================================
cat("\n=================================================================\n")
cat("TEST 1: FGS Layer 1 - Multiple Methods Including NMF\n")
cat("=================================================================\n")

fgs_methods <- c(
    "random_forest_ranger", # Baseline
    "nmf_loadings", # Core method
    "lasso", # Regularization
    "ridge", # Regularization
    "pca_loadings" # Dim reduction
)

cat("Methods to test (Single Pass):", paste(fgs_methods, collapse = ", "), "\n\n")

fgs_results <- list()
timing_start <- Sys.time()

tryCatch(
    {
        # [OPTIMIZATION] Run all methods in ONE call to avoid redundant preprocessing
        fgs_result_all <- FGS(
            data = sobj,
            target_var = target_var,
            method = fgs_methods,
            n_features = 20,
            fgs_seed = 42
        )

        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))

        cat("✓ FGS Batch Completed in", round(elapsed, 2), "seconds\n")

        # Check results
        for (m in fgs_methods) {
            if (!is.null(fgs_result_all[[m]])) {
                cat("  -", m, ": Success (", length(fgs_result_all[[m]]$genes), "genes)\n")
                fgs_results[[m]] <- fgs_result_all[[m]]
            } else {
                cat("  -", m, ": Failed or Returned NULL\n")
            }
        }
    },
    error = function(e) {
        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))
        cat("✗ FGS Batch FAILED:", conditionMessage(e), "\n")
        cat("  Time elapsed:", round(elapsed, 2), "seconds\n")
    }
)

# Summary
cat("\n--- FGS Test Summary ---\n")
success_count <- length(fgs_results)
cat("Success rate:", success_count, "/", length(fgs_methods), "\n")

# Select signature for TML
# Use the first successful signature
if (length(fgs_results) > 0) {
    first_sig_name <- names(fgs_results)[1]
    fgs_signature <- fgs_results[[first_sig_name]]

    cat("Selected signature for TML:", first_sig_name, "\n")

    # [FIX] Wrap single signature in a named list so TML7 treats it correctly
    l1_sigs_for_tml <- list()
    l1_sigs_for_tml[[first_sig_name]] <- fgs_signature

    # ==============================================================================
    # TEST 2: TML Layer 2
    # ==============================================================================
    cat("\n=================================================================\n")
    cat("TEST 2: TML Layer 2 - All Methods\n")
    cat("=================================================================\n")

    l2_methods <- c("glm", "ranger", "xgbTree", "svmRadial", "mlp", "earth", "nnet", "glmnet")
    cat("L2 Methods to test:", paste(l2_methods, collapse = ", "), "\n\n")

    tryCatch(
        {
            tml_model <- TML7(
                l1_signatures = l1_sigs_for_tml,
                holdout_data = sobj,
                target_var = target_var,
                l2_methods = l2_methods,
                cv_folds = 5,
                cv_method = "cv", # Standard CV
                fgs_seed = 42
            )
            cat("✓ TML7 Completed Successfully\n")
        },
        error = function(e) {
            cat("✗ FAILED: TML7\n")
            cat("  Error:", conditionMessage(e), "\n")
        }
    )

    # ==============================================================================
    # TEST 3: LOGO CV
    # ==============================================================================
    cat("\n=================================================================\n")
    cat("TEST 3: LOGO CV\n")
    cat("=================================================================\n")

    if (!is.null(group_var)) {
        cat("Number of groups (patients):", length(unique(sobj@meta.data[[group_var]])), "\n")
        tryCatch(
            {
                tml_logo <- TML7(
                    l1_signatures = l1_sigs_for_tml,
                    holdout_data = sobj,
                    target_var = target_var,
                    l2_methods = c("glm", "ranger"), # Test with simple model first
                    cv_method = "LOGO",
                    cv_group_var = group_var,
                    fgs_seed = 42
                )
                cat("✓ LOGO CV Completed Successfully\n")

                # [NEW] Outlier Analysis
                if (exists("analyze_tml_outliers")) {
                    cat("\n--- Outlier Analysis (LOGO) ---\n")
                    try(
                        {
                            outliers <- analyze_tml_outliers(tml_logo)
                            if (!is.null(outliers)) print(head(outliers))
                        },
                        silent = TRUE
                    )
                }
            },
            error = function(e) {
                cat("✗ FAILED: LOGO CV\n")
                cat("  Error:", conditionMessage(e), "\n")
            }
        )
    } else {
        cat("Skipping LOGO CV (No group variable)\n")
    }
} else {
    cat("Skipping TML7 (No FGS results)\n")
}

cat("\n=================================================================\n")
cat("Testing Complete!\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
