#!/usr/bin/env Rscript
# Test script for run_milo_pipeline with is5 and is5s data
# Run from /home/user3/GJC_KDW_250721: Rscript /data/user3/git_repo/_wt/milo/scripts/fgs/test_milo_is5.R

# Source start.R to activate renv and load packages
start_r_paths <- c(
    "start.R",
    "../start.R",
    "/home/user3/GJC_KDW_250721/start.R"
)
start_r_found <- FALSE
for (path in start_r_paths) {
    if (file.exists(path)) {
        source(path)
        start_r_found <- TRUE
        break
    }
}

if (!start_r_found) {
    stop("start.R not found. Please run from /home/user3/GJC_KDW_250721 or ensure start.R is accessible.")
}

# Load myR package
if (!"package:myR" %in% search()) {
    myr_paths <- c(
        "../myR",
        "../../mylit/myR",
        "/data/user3/git_repo/mylit/myR",
        "/data/user3/git_repo/_wt/milo/myR"
    )
    myr_found <- FALSE
    for (path in myr_paths) {
        if (file.exists(path) && file.exists(file.path(path, "R"))) {
            devtools::load_all(path)
            myr_found <- TRUE
            break
        }
    }
    if (!myr_found) {
        stop("myR package not found. Tried paths: ", paste(myr_paths, collapse = ", "))
    }
}

# Test function
test_milo_on_dataset <- function(dataset_name, seurat_qs_path, output_dir) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("Testing Milo Pipeline on:", dataset_name, "\n")
    cat("Data path:", seurat_qs_path, "\n")
    cat("Output directory:", output_dir, "\n")
    cat(rep("=", 80), "\n\n")
    
    if (!file.exists(seurat_qs_path)) {
        cat("ERROR: Data file not found:", seurat_qs_path, "\n")
        return(NULL)
    }
    
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Test parameters based on guide.md
    # For stroke PBMC data:
    # - g3 is target variable ("2", "1", "NA")
    # - hos_no is patient variable (8-digit number)
    # - anno3.scvi is cluster variable (~22 types)
    # - GEM or set can be batch variable
    test_params <- list(
        seurat_qs_path = seurat_qs_path,
        patient_var = "hos_no",
        cluster_var = "anno3.scvi",
        target_var = "g3",
        batch_var = "GEM",  # Using GEM as batch variable
        graph_reduction = "integrated.scvi",
        layout_reduction = "umap.scvi",
        k = 30,
        d = 30,
        prop = 0.1,
        alpha = 0.1,
        target_include = c("1", "2"),  # Exclude "NA" for binary comparison
        target_levels = c("1", "2"),   # Set factor order: "1" as reference, "2" as test
        output_dir = output_dir,
        prefix = paste0("milo_", dataset_name),
        save = TRUE,
        plotting = TRUE,
        verbose = TRUE
    )
    
    # Run test
    result <- NULL
    error_occurred <- FALSE
    
    tryCatch({
        result <- do.call(run_milo_pipeline, test_params)
        
        cat("\n=== Test Results ===\n")
        cat("✓ Pipeline completed successfully!\n")
        cat("Milo object class:", class(result$milo)[1], "\n")
        cat("DA results:", nrow(result$da_results), "neighborhoods\n")
        cat("Plots:", if(is.null(result$plots)) "None" else length(result$plots), "\n")
        
        # Check cache files
        cat("\n=== Cache Files ===\n")
        cache_files <- list(
            nhoods = file.path(output_dir, paste0("milo_", dataset_name, "_01_nhoods_built.qs")),
            distances = file.path(output_dir, paste0("milo_", dataset_name, "_02_distances_calculated.qs")),
            da_milo = file.path(output_dir, paste0("milo_", dataset_name, "_03_tested.qs")),
            da_results = file.path(output_dir, paste0("milo_", dataset_name, "_03_da_results.qs")),
            plots = file.path(output_dir, paste0("milo_", dataset_name, "_04_plots.qs"))
        )
        
        for (name in names(cache_files)) {
            path <- cache_files[[name]]
            exists <- file.exists(path)
            size <- if(exists) file.info(path)$size else 0
            cat(sprintf("  %s: %s (%.2f MB)\n", name, if(exists) "✓" else "✗", size/1024/1024))
        }
        
        # Summary of DA results
        if (!is.null(result$da_results) && nrow(result$da_results) > 0) {
            cat("\n=== DA Results Summary ===\n")
            cat("Total neighborhoods:", nrow(result$da_results), "\n")
            sig_count <- sum(result$da_results$SpatialFDR < 0.1, na.rm = TRUE)
            cat("Significant neighborhoods (SpatialFDR < 0.1):", sig_count, "\n")
            cat("Mean logFC:", mean(result$da_results$logFC, na.rm = TRUE), "\n")
            cat("SD logFC:", sd(result$da_results$logFC, na.rm = TRUE), "\n")
        }
        
    }, error = function(e) {
        cat("\n=== Test Failed ===\n")
        cat("Error:", conditionMessage(e), "\n")
        cat("\nStack trace:\n")
        traceback()
        error_occurred <<- TRUE
    })
    
    if (error_occurred) {
        return(NULL)
    }
    
    return(result)
}

# Main execution
cat("=== Milo Pipeline Test Suite ===\n")
cat("Testing on is5s and is5 datasets\n\n")

# Test on is5s (downsampled data - faster for testing)
is5s_path <- "/data/user3/sobj/IS6_sex_added_0.1x_251110.qs"
is5s_output <- "/data/user3/sobj/milo_is5s"
result_is5s <- test_milo_on_dataset("is5s", is5s_path, is5s_output)

# Test on is5 (full data)
is5_path <- "/data/user3/sobj/IS6_sex_added_251110.qs"
is5_output <- "/data/user3/sobj/milo_is5"
result_is5 <- test_milo_on_dataset("is5", is5_path, is5_output)

# Final summary
cat("\n", rep("=", 80), "\n", sep = "")
cat("=== Final Summary ===\n")
cat("is5s test:", if(is.null(result_is5s)) "FAILED" else "PASSED", "\n")
cat("is5 test:", if(is.null(result_is5)) "FAILED" else "PASSED", "\n")
cat(rep("=", 80), "\n")

if (is.null(result_is5s) || is.null(result_is5)) {
    quit(status = 1)
}

cat("\n=== All Tests Passed ===\n")

