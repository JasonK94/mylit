#!/usr/bin/env Rscript
# Test script for run_milo_pipeline
# Run this from st/ directory: cd st && Rscript ../_wt/milo/test_milo_from_st.R

# Source start.R to activate renv and load packages
if (file.exists("start.R")) {
    source("start.R")
} else {
    stop("Please run this script from the st/ directory")
}

# Load myR package (should be loaded by start.R, but ensure it's available)
if (!"package:myR" %in% search()) {
    if (file.exists("../myR")) {
        devtools::load_all("../myR")
    } else {
        stop("myR package not found")
    }
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
result <- NULL
error_occurred <- FALSE

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
    cat("\nStack trace:\n")
    traceback()
    error_occurred <<- TRUE
})

if (error_occurred) {
    quit(status = 1)
}

cat("\n=== Test Complete ===\n")

