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
    conflicted::conflicts_prefer(dplyr::select, quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::summarise, quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::filter, quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::mutate, quiet = TRUE)
    conflicted::conflicts_prefer(dplyr::arrange, quiet = TRUE)
    conflicted::conflicts_prefer(base::intersect, quiet = TRUE)
    conflicted::conflicts_prefer(base::setdiff, quiet = TRUE)
    conflicted::conflicts_prefer(base::union, quiet = TRUE)
    conflicted::conflicts_prefer(base::colMeans, quiet = TRUE)
    conflicted::conflicts_prefer(base::colSums, quiet = TRUE)
    conflicted::conflicts_prefer(base::rowMeans, quiet = TRUE)
    conflicted::conflicts_prefer(base::rowSums, quiet = TRUE)
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
if (!is.null(fgs_results$random_forest_ranger)) {
    fgsa <- fgs_results$random_forest_ranger
} else {
    fgsa <- fgs_results[[which(!sapply(fgs_results, is.null))[1]]]
}

# ==============================================================================
# TEST 2: TML7 with All L2 Methods
# ==============================================================================
cat("\n=================================================================\n")
cat("TEST 2: TML Layer 2 - All Methods\n")
cat("=================================================================\n")

l2_methods <- c(
    "glm", "ranger", "xgbTree", "svmRadial", "mlp", "earth"
)

cat("L2 Methods to test:", paste(l2_methods, collapse = ", "), "\n\n")

if (!is.null(fgsa)) {
    tryCatch(
        {
            tmla_full <- TML7(
                l1_signatures = fgsa,
                holdout_data = sobj,
                target_var = target_var,
                l2_methods = l2_methods,
                cv_folds = 5,
                cv_method = "cv",
                cv_group_var = group_var,
                fgs_seed = 42 # Correct argument name
            )
            cat("✓ SUCCESS: TML7 completed\n")

            # Metrics
            if (length(tmla_full$results) > 0) {
                cat("Performance (First Signature):\n")
                print(tmla_full$results[[1]]$models$glm$performance)
            }
        },
        error = function(e) {
            cat("✗ FAILED: TML7 -", conditionMessage(e), "\n")
            tmla_full <- NULL
        }
    )
} else {
    cat("Skipping TML7 (No FGS results)\n")
    tmla_full <- NULL
}

# ==============================================================================
# TEST 3: CV Method Comparison (LOGO)
# ==============================================================================
cat("\n=================================================================\n")
cat("TEST 3: LOGO CV\n")
cat("=================================================================\n")

if (!is.null(fgsa) && group_var %in% colnames(sobj@meta.data)) {
    # LOGO requires group variable
    n_groups <- length(unique(sobj@meta.data[[group_var]]))
    cat("Number of groups (patients):", n_groups, "\n")

    if (n_groups >= 3) { # Need at least a few groups for LOGO
        tryCatch(
            {
                tmla_logo <- TML7(
                    l1_signatures = fgsa,
                    holdout_data = sobj,
                    target_var = target_var,
                    l2_methods = c("glm", "ranger"),
                    cv_method = "LOGO",
                    cv_group_var = group_var,
                    fgs_seed = 42
                )
                cat("✓ SUCCESS: LOGO CV completed\n")

                if (!is.null(tmla_logo$results[[1]]$cv_folds)) {
                    cat("  Folds created:", length(tmla_logo$results[[1]]$cv_folds$index), "\n")
                }
            },
            error = function(e) {
                cat("✗ FAILED: LOGO CV -", conditionMessage(e), "\n")
            }
        )
    } else {
        cat("Skipping LOGO: Not enough groups (<3)\n")
    }
}

cat("\n=================================================================\n")
cat("Testing Complete!\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
