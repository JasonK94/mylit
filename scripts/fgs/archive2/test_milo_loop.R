#!/usr/bin/env Rscript
# Test script for run_milo_pipeline
# This script tests the pipeline with the provided Seurat object and cache files
# 
# Usage: Run from st/ directory after sourcing start.R
# Or: Rscript test_milo_loop.R (if renv is activated)

# Try to activate renv if available
# First try st/start.R which should activate renv
if (file.exists("../mylit/st/start.R")) {
    source("../mylit/st/start.R")
} else if (file.exists("../../mylit/st/start.R")) {
    source("../../mylit/st/start.R")
} else if (file.exists("/home/user3/data_user3/git_repo/mylit/st/start.R")) {
    source("/home/user3/data_user3/git_repo/mylit/st/start.R")
} else if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
} else if (file.exists("../renv/activate.R")) {
    source("../renv/activate.R")
}

# Check if myR is already loaded (start.R may have loaded it)
if (!"package:myR" %in% search()) {
    # start.R may have tried to load it but failed due to missing dependencies
    # Try to load myR with devtools, which may be more lenient
    myr_paths <- c(
        "/home/user3/data_user3/git_repo/mylit/myR",
        "../../mylit/myR",
        "../myR",
        "myR"
    )
    
    myr_found <- FALSE
    for (path in myr_paths) {
        if (file.exists(path) && file.exists(file.path(path, "R"))) {
            if (requireNamespace("devtools", quietly = TRUE)) {
                cat("Loading myR from:", path, "\n")
                # Use load_all with reset=FALSE to avoid re-checking imports
                tryCatch({
                    devtools::load_all(path, reset = FALSE, quiet = TRUE)
                    myr_found <- TRUE
                    break
                }, error = function(e) {
                    # If load_all fails, try to source the functions directly
                    cat("load_all failed, trying alternative loading method...\n")
                })
            }
        }
    }
    
    if (!myr_found) {
        # Last resort: try to source the milo_pipeline.R file directly
        milo_pipeline_path <- file.path("/home/user3/data_user3/git_repo/_wt/milo/myR/R/milo_pipeline.R")
        if (file.exists(milo_pipeline_path)) {
            cat("Sourcing milo_pipeline.R directly...\n")
            # First ensure required packages are loaded
            required_pkgs <- c("Seurat", "SingleCellExperiment", "miloR", "Matrix", "qs", "cli", "ggplot2", "patchwork", "scater")
            for (pkg in required_pkgs) {
                if (!requireNamespace(pkg, quietly = TRUE)) {
                    cat("Warning: Package", pkg, "not available. Some functions may fail.\n")
                } else {
                    # Load the namespace
                    tryCatch({
                        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
                    }, error = function(e) {
                        cat("Warning: Could not load package", pkg, "\n")
                    })
                }
            }
            source(milo_pipeline_path)
            myr_found <- TRUE
        } else {
            stop("myR package not found. Tried paths: ", paste(myr_paths, collapse = ", "))
        }
    }
} else {
    cat("myR package already loaded.\n")
}

# Set paths
seurat_qs_path <- "/data/user3/git_repo/_wt/milo/myeloid_v1_TNBC.qs"
output_dir <- "/data/user3/git_repo/_wt/milo"
cache_nhoods <- file.path(output_dir, "milo_01_nhoods_built.qs")
cache_distances <- file.path(output_dir, "milo_02_distances_calculated.qs")

cat("=== Milo Pipeline Test ===\n")
cat("Seurat object:", seurat_qs_path, "\n")
cat("Output directory:", output_dir, "\n")
cat("Cache files:\n")
cat("  - Nhoods:", cache_nhoods, if(file.exists(cache_nhoods)) "✓" else "✗", "\n")
cat("  - Distances:", cache_distances, if(file.exists(cache_distances)) "✓" else "✗", "\n")
cat("\n")

# Test parameters
test_params <- list(
    seurat_qs_path = seurat_qs_path,
    patient_var = "Core",
    cluster_var = "Annot3",
    target_var = "intra_T",
    batch_var = "batch",
    graph_reduction = "integrated.scvi",
    layout_reduction = "umap.scvi",
    k = 30,
    d = 10,
    output_dir = output_dir,
    prefix = "milo",
    save = TRUE,
    plotting = FALSE,  # Skip plotting for faster testing
    verbose = TRUE
)

# Run test
cat("Starting pipeline test...\n")
tryCatch({
    result <- do.call(run_milo_pipeline, test_params)
    
    cat("\n=== Test Results ===\n")
    cat("✓ Pipeline completed successfully!\n")
    cat("Milo object:", class(result$milo), "\n")
    cat("DA results:", nrow(result$da_results), "neighborhoods\n")
    cat("Plots:", if(is.null(result$plots)) "None" else length(result$plots), "\n")
    
    # Check cache files
    cat("\n=== Cache Files ===\n")
    cache_files <- list(
        nhoods = file.path(output_dir, "milo_01_nhoods_built.qs"),
        distances = file.path(output_dir, "milo_02_distances_calculated.qs"),
        da_milo = file.path(output_dir, "milo_03_tested.qs"),
        da_results = file.path(output_dir, "milo_03_da_results.qs")
    )
    
    for (name in names(cache_files)) {
        path <- cache_files[[name]]
        exists <- file.exists(path)
        size <- if(exists) file.info(path)$size else 0
        cat(sprintf("  %s: %s (%.2f MB)\n", name, if(exists) "✓" else "✗", size/1024/1024))
    }
    
}, error = function(e) {
    cat("\n=== Test Failed ===\n")
    cat("Error:", conditionMessage(e), "\n")
    traceback()
    quit(status = 1)
})

cat("\n=== Test Complete ===\n")

