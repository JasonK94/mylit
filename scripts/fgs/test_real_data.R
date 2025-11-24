#!/usr/bin/env Rscript
# ==============================================================================
# Real Data Testing Script for FGS/TML Framework
# ==============================================================================
# Purpose: Validate recently patched features on actual single-cell data
# Author: Automated testing script
# Date: 2025-11-24
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

# Load FGS/TML functions
cat("Loading FGS/TML functions...\n")
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R")
source("/home/user3/data_user3/git_repo/_wt/fgs/myR/R/utils_tml.R")
cat("✓ Functions loaded\n\n")

# --- Load Real Data ---
cat("\n>>> Step 2: Loading Real Data\n")
data_path <- "/home/user3/data_user3/sobj/IS_scvi_251107_ds2500.qs"
cat("Data path:", data_path, "\n")

if (!file.exists(data_path)) {
    stop("Data file not found: ", data_path)
}

sobj_full <- qs::qread(data_path)
cat("Full dataset loaded:", ncol(sobj_full), "cells\n")

# Use all cells from downsampled data (already 2500 cells)
# No further subsetting needed for reasonable runtime
sobj <- sobj_full
cat("Using full downsampled dataset:", ncol(sobj), "cells\n")

# Check metadata
cat("\nAvailable metadata columns:\n")
print(colnames(sobj@meta.data))

# Determine target variable (g3 or alternative)
if ("g3" %in% colnames(sobj@meta.data)) {
    target_var <- "g3"
} else if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
    target_var <- "seurat_clusters"
    # Convert to binary if needed
    sobj@meta.data$target_binary <- ifelse(
        as.numeric(as.character(sobj@meta.data[[target_var]])) %% 2 == 0,
        "class_A", "class_B"
    )
    target_var <- "target_binary"
} else {
    stop("No suitable target variable found in metadata")
}

cat("Target variable:", target_var, "\n")
cat("Target distribution:\n")
print(table(sobj@meta.data[[target_var]]))

# Determine group variable (emrid or alternative)
if ("emrid" %in% colnames(sobj@meta.data)) {
    group_var <- "emrid"
} else if ("orig.ident" %in% colnames(sobj@meta.data)) {
    group_var <- "orig.ident"
} else {
    # Create synthetic groups
    sobj@meta.data$synthetic_group <- sample(
        paste0("Group_", 1:10),
        ncol(sobj),
        replace = TRUE
    )
    group_var <- "synthetic_group"
}

cat("Group variable:", group_var, "\n")
cat("Number of groups:", length(unique(sobj@meta.data[[group_var]])), "\n\n")

# ==============================================================================
# TEST 1: FGS with All Methods Including NMF
# ==============================================================================
cat("=================================================================\n")
cat("TEST 1: FGS Layer 1 - Multiple Methods Including NMF\n")
cat("=================================================================\n")

fgs_methods <- c(
    "random_forest_ranger", # Baseline (already validated)
    "nmf_loadings", # Core method we fixed
    "lasso", # Regularization method
    "ridge", # Regularization method
    "pca_loadings" # Dimensionality reduction
)

cat("Methods to test:", paste(fgs_methods, collapse = ", "), "\n\n")

fgs_results <- list()
fgs_timings <- list()

for (method in fgs_methods) {
    cat("\n>>> Testing method:", method, "\n")

    timing_start <- Sys.time()

    tryCatch(
        {
            fgs_result <- FGS(
                data = sobj,
                target_var = target_var,
                method = method,
                n_features = 20,
                fgs_seed = 42
            )

            timing_end <- Sys.time()
            elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))

            fgs_results[[method]] <- fgs_result
            fgs_timings[[method]] <- elapsed

            cat("✓ SUCCESS:", method, "completed in", round(elapsed, 2), "seconds\n")
            cat("  Signatures found:", length(fgs_result$signatures), "\n")

            if (length(fgs_result$signatures) > 0) {
                sig1 <- fgs_result$signatures[[1]]
                cat("  First signature genes:", length(sig1$genes), "\n")
                cat("  Sample genes:", paste(head(sig1$genes, 3), collapse = ", "), "\n")
            }
        },
        error = function(e) {
            timing_end <- Sys.time()
            elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))

            fgs_results[[method]] <- NULL
            fgs_timings[[method]] <- elapsed

            cat("✗ FAILED:", method, "\n")
            cat("  Error:", conditionMessage(e), "\n")
            cat("  Time elapsed:", round(elapsed, 2), "seconds\n")
        }
    )
}

# Summary of FGS test
cat("\n--- FGS Test Summary ---\n")
success_count <- sum(sapply(fgs_results, function(x) !is.null(x)))
cat("Success rate:", success_count, "/", length(fgs_methods), "\n")
cat("Average time per method:", round(mean(unlist(fgs_timings)), 2), "seconds\n\n")

# Use successful FGS result for TML tests (prefer ranger or nmf)
if (!is.null(fgs_results$random_forest_ranger)) {
    fgsa <- fgs_results$random_forest_ranger
    cat("Using 'random_forest_ranger' signatures for TML tests\n")
} else if (!is.null(fgs_results$nmf_loadings)) {
    fgsa <- fgs_results$nmf_loadings
    cat("Using 'nmf_loadings' signatures for TML tests\n")
} else {
    # Use first successful result
    fgsa <- fgs_results[[which(!sapply(fgs_results, is.null))[1]]]
    cat("Using first successful method for TML tests\n")
}

# ==============================================================================
# TEST 2: TML7 with All L2 Methods
# ==============================================================================
cat("\n=================================================================\n")
cat("TEST 2: TML Layer 2 - All Methods\n")
cat("=================================================================\n")

l2_methods <- c(
    "glm", # Baseline
    "ranger", # Already validated
    "xgbTree", # Already validated
    "svmRadial", # Recently fixed
    "mlp", # Recently fixed
    "earth" # Recently fixed
)

cat("L2 Methods to test:", paste(l2_methods, collapse = ", "), "\n\n")

cat(">>> Running TML7 with all L2 methods and standard CV...\n")
timing_start <- Sys.time()

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
            fgs_seed = 42
        )

        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "mins"))

        cat("✓ SUCCESS: TML7 completed in", round(elapsed, 2), "minutes\n")

        # Analyze results
        cat("\n--- TML7 Results Summary ---\n")
        cat("Number of signatures tested:", length(tmla_full$results), "\n")

        if (length(tmla_full$results) > 0) {
            first_sig <- tmla_full$results[[1]]
            cat("Models trained per signature:", length(first_sig$models), "\n")
            cat("CV folds recorded:", !is.null(first_sig$cv_folds), "\n")

            # Check performance metrics
            cat("\nPerformance Metrics (First Signature):\n")
            for (model_name in names(first_sig$models)) {
                perf <- first_sig$models[[model_name]]$performance
                cat(sprintf(
                    "  %s: ROC=%.3f, Sens=%.3f, Spec=%.3f\n",
                    model_name,
                    perf$ROC %||% NA,
                    perf$Sens %||% NA,
                    perf$Spec %||% NA
                ))
            }
        }

        # Visualize metrics
        cat("\n>>> Generating performance plots...\n")
        pdf("logs/fgs/test_real_data_metrics.pdf", width = 10, height = 6)
        plot_tml_metrics(tmla_full, metric = "ROC")
        plot_tml_metrics(tmla_full, metric = "Sens")
        plot_tml_metrics(tmla_full, metric = "Spec")
        dev.off()
        cat("Plots saved to: logs/fgs/test_real_data_metrics.pdf\n")

        # Analyze outliers
        cat("\n>>> Analyzing outliers (IQR method)...\n")
        outliers_roc <- analyze_tml_outliers(tmla_full, metric = "ROC", threshold_method = "iqr")
        outliers_sens <- analyze_tml_outliers(tmla_full, metric = "Sens", threshold_method = "iqr")

        cat("Outlier folds (ROC):\n")
        print(outliers_roc$summary)
        cat("\nOutlier folds (Sens):\n")
        print(outliers_sens$summary)
    },
    error = function(e) {
        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "mins"))

        cat("✗ FAILED: TML7 with all methods\n")
        cat("  Error:", conditionMessage(e), "\n")
        cat("  Time elapsed:", round(elapsed, 2), "minutes\n")

        tmla_full <- NULL
    }
)

# ==============================================================================
# TEST 3: CV Method Comparison
# ==============================================================================
cat("\n=================================================================\n")
cat("TEST 3: CV Method Integration Test\n")
cat("=================================================================\n")

# Use fast methods for CV comparison
fast_methods <- c("glm", "ranger")

# --- Test 3a: Standard K-Fold CV ---
cat("\n>>> Test 3a: Standard K-Fold CV\n")
timing_start <- Sys.time()

tryCatch(
    {
        tmla_kfold <- TML7(
            l1_signatures = fgsa,
            holdout_data = sobj,
            target_var = target_var,
            l2_methods = fast_methods,
            cv_folds = 5,
            cv_method = "cv",
            cv_group_var = group_var,
            fgs_seed = 42
        )

        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))

        cat("✓ SUCCESS: K-Fold CV completed in", round(elapsed, 2), "seconds\n")
        cat("  CV folds structure:\n")
        if (!is.null(tmla_kfold$results[[1]]$cv_folds)) {
            cat("    Number of folds:", length(tmla_kfold$results[[1]]$cv_folds$index), "\n")
        }
    },
    error = function(e) {
        cat("✗ FAILED: K-Fold CV -", conditionMessage(e), "\n")
        tmla_kfold <- NULL
    }
)

# --- Test 3b: LOGO CV ---
cat("\n>>> Test 3b: LOGO CV\n")

# Create subset with fewer groups for LOGO
n_groups <- length(unique(sobj@meta.data[[group_var]]))
cat("Total groups available:", n_groups, "\n")

if (n_groups < 20) {
    # Use full dataset
    sobj_logo <- sobj
    cat("Using full dataset for LOGO (", n_groups, "groups)\n")
} else {
    # Subset to ~15 groups
    set.seed(42)
    selected_groups <- sample(unique(sobj@meta.data[[group_var]]), 15)
    cells_logo <- rownames(sobj@meta.data)[sobj@meta.data[[group_var]] %in% selected_groups]
    sobj_logo <- subset(sobj, cells = cells_logo)
    cat("Subset created with", length(selected_groups), "groups\n")
}

timing_start <- Sys.time()

tryCatch(
    {
        tmla_logo <- TML7(
            l1_signatures = fgsa,
            holdout_data = sobj_logo,
            target_var = target_var,
            l2_methods = fast_methods,
            cv_method = "LOGO",
            cv_group_var = group_var,
            fgs_seed = 42
        )

        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))

        cat("✓ SUCCESS: LOGO CV completed in", round(elapsed, 2), "seconds\n")
        cat("  CV folds structure:\n")
        if (!is.null(tmla_logo$results[[1]]$cv_folds)) {
            cat("    Number of folds (groups):", length(tmla_logo$results[[1]]$cv_folds$index), "\n")
        }
    },
    error = function(e) {
        cat("✗ FAILED: LOGO CV -", conditionMessage(e), "\n")
        tmla_logo <- NULL
    }
)

# --- Test 3c: Repeated CV ---
cat("\n>>> Test 3c: Repeated CV\n")
timing_start <- Sys.time()

tryCatch(
    {
        tmla_repeated <- TML7(
            l1_signatures = fgsa,
            holdout_data = sobj,
            target_var = target_var,
            l2_methods = fast_methods,
            cv_folds = 5,
            cv_method = "repeatedcv",
            repeats = 3,
            cv_group_var = group_var,
            fgs_seed = 42
        )

        timing_end <- Sys.time()
        elapsed <- as.numeric(difftime(timing_end, timing_start, units = "secs"))

        cat("✓ SUCCESS: Repeated CV completed in", round(elapsed, 2), "seconds\n")
        cat("  CV folds structure:\n")
        if (!is.null(tmla_repeated$results[[1]]$cv_folds)) {
            cat(
                "    Number of folds (5 folds × 3 repeats = 15):",
                length(tmla_repeated$results[[1]]$cv_folds$index), "\n"
            )
        }

        # Analyze outliers for repeated CV
        cat("\n  Analyzing outliers in Repeated CV...\n")
        outliers_repeated <- analyze_tml_outliers(tmla_repeated, metric = "ROC", threshold_method = "iqr")
        cat("  Problematic folds found:", nrow(outliers_repeated$summary), "\n")
    },
    error = function(e) {
        cat("✗ FAILED: Repeated CV -", conditionMessage(e), "\n")
        tmla_repeated <- NULL
    }
)

# ==============================================================================
# Final Summary and Save Results
# ==============================================================================
cat("\n=================================================================\n")
cat("FINAL SUMMARY\n")
cat("=================================================================\n")

cat("\n1. FGS Layer 1 Test:\n")
cat("   Success rate:", success_count, "/", length(fgs_methods), "\n")
cat(
    "   Successful methods:",
    paste(names(fgs_results)[!sapply(fgs_results, is.null)], collapse = ", "), "\n"
)

cat("\n2. TML Layer 2 Test:\n")
if (!is.null(tmla_full)) {
    cat("   ✓ All L2 methods tested successfully\n")
    cat("   Methods:", paste(l2_methods, collapse = ", "), "\n")
} else {
    cat("   ✗ TML Layer 2 test failed\n")
}

cat("\n3. CV Method Integration:\n")
cat("   Standard K-Fold:", ifelse(!is.null(tmla_kfold), "✓", "✗"), "\n")
cat("   LOGO:", ifelse(!is.null(tmla_logo), "✓", "✗"), "\n")
cat("   Repeated CV:", ifelse(!is.null(tmla_repeated), "✓", "✗"), "\n")

# Save all results
cat("\n>>> Saving results...\n")
results_file <- "logs/fgs/test_real_data_results.rds"
saveRDS(
    list(
        fgs_results = fgs_results,
        fgs_timings = fgs_timings,
        tmla_full = tmla_full,
        tmla_kfold = tmla_kfold,
        tmla_logo = tmla_logo,
        tmla_repeated = tmla_repeated,
        test_info = list(
            data_path = data_path,
            n_cells = ncol(sobj),
            target_var = target_var,
            group_var = group_var,
            test_date = Sys.time()
        )
    ),
    results_file
)

cat("Results saved to:", results_file, "\n")

cat("\n=================================================================\n")
cat("Testing Complete!\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=================================================================\n")
