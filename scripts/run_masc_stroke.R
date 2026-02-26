#!/usr/bin/env Rscript
# MASC Analysis for Stroke Data
# Tests cluster abundance differences by g3 (target variable)

# Load required libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
})

# Source MASC functions
masc_r_path <- "/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/myR/R/masc.R"
if (!file.exists(masc_r_path)) stop("MASC.R not found")
source(masc_r_path)

# Configuration
DATA_PATH_1 <- "/data/user3/sobj/is2_IS_3_clustered.qs"
DATA_PATH_2 <- "/data/user3/sobj/is2_IS_6_subsets_Monocytes.qs"
OUTPUT_DIR <- "/data/user3/sobj/masc/stroke"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== MASC Analysis for Stroke Data ===\n\n")

# Function to prepare data and run MASC
run_masc_for_cluster <- function(seurat_obj, cluster_var, data_name, output_suffix) {
    cat(sprintf("\n--- Analyzing %s with cluster variable: %s ---\n", data_name, cluster_var))
    
    # Check if cluster variable exists
    if (!cluster_var %in% colnames(seurat_obj@meta.data)) {
        cat(sprintf("WARNING: Cluster variable '%s' not found. Skipping.\n", cluster_var))
        return(NULL)
    }
    
    # Check required variables
    if (!"g3" %in% colnames(seurat_obj@meta.data)) {
        cat("WARNING: Target variable 'g3' not found.\n")
        return(NULL)
    }
    
    # Prepare Model Variables
    # Fixed effects: GEM only (as per request)
    fixed_effects <- NULL
    if ("GEM" %in% colnames(seurat_obj@meta.data)) {
        fixed_effects <- c("GEM")
        seurat_obj@meta.data$GEM <- as.factor(seurat_obj@meta.data$GEM)
    }
    
    # Random effects: hos_no
    random_effects <- "hos_no"
    seurat_obj@meta.data$hos_no <- as.character(seurat_obj@meta.data$hos_no) # Ensure character/factor
    
    # Ensure g3 is factor
    seurat_obj@meta.data$g3 <- as.factor(seurat_obj@meta.data$g3)
    
    cat(sprintf("Fixed effects: %s\n", paste(fixed_effects, collapse = ", ")))
    cat(sprintf("Random effects: %s\n", paste(random_effects, collapse = ", ")))
    
    # Check cluster counts
    cluster_counts <- table(seurat_obj@meta.data[[cluster_var]], seurat_obj@meta.data$g3)
    cat("\nCluster counts by g3:\n")
    print(cluster_counts)
    
    # Filter clusters
    min_cells <- 10
    valid_clusters <- rownames(cluster_counts)[rowSums(cluster_counts) >= min_cells]
    
    if (length(valid_clusters) < 2) {
        cat("WARNING: Not enough valid clusters.\n")
        return(NULL)
    }
    
    cat(sprintf("\nUsing %d clusters (>= %d cells)\n", length(valid_clusters), min_cells))
    
    # Subset
    seurat_obj_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[cluster_var]] %in% valid_clusters])
    
    # Run MASC pipeline
    cat("\nRunning MASC pipeline...\n")
    results <- tryCatch({
        res <- run_masc_pipeline(
            seurat_obj = seurat_obj_subset,
            cluster_var = cluster_var,
            contrast_var = "g3",
            random_effects = random_effects,
            fixed_effects = fixed_effects,
            save = TRUE,
            output_dir = file.path(OUTPUT_DIR, data_name),
            prefix = paste0("masc_", output_suffix),
            force_run = TRUE,
            plotting = TRUE,
            save_models = FALSE,
            adjust_pvalue = TRUE,
            verbose = TRUE
        )
        res
    }, error = function(e) {
        cat(sprintf("\nERROR: %s\n", e$message))
        return(NULL)
    })
    
    if (!is.null(results)) {
        cat("\n=== MASC Results ===\n")
        print(results$masc_results)
        
        # Save summary text
        summary_path <- file.path(OUTPUT_DIR, data_name, paste0("masc_", output_suffix, "_summary.txt"))
        sink(summary_path)
        print(results$masc_results)
        sink()
        cat(sprintf("\nResults saved to: %s\n", summary_path))
    }
    
    return(results)
}

# === Analysis 1: Full data with anno3big ===
cat("\n\n===== Analysis 1: Full Data - anno3big =====\n")
if (file.exists(DATA_PATH_1)) {
    seurat_obj_1 <- qs::qread(DATA_PATH_1)
    run_masc_for_cluster(seurat_obj_1, "anno3big", "full_data", "anno3big")
    
    # === Analysis 2: Full data with anno3 ===
    cat("\n\n===== Analysis 2: Full Data - anno3 =====\n")
    run_masc_for_cluster(seurat_obj_1, "anno3", "full_data", "anno3")
} else {
    cat(sprintf("WARNING: Data file not found: %s\n", DATA_PATH_1))
}

# === Analysis 3: Monocyte subset with anno.mo ===
cat("\n\n===== Analysis 3: Monocyte Subset - anno.mo =====\n")
if (file.exists(DATA_PATH_2)) {
    seurat_obj_2 <- qs::qread(DATA_PATH_2)
    run_masc_for_cluster(seurat_obj_2, "anno.mo", "monocyte_subset", "anno_mo")
} else {
    cat(sprintf("WARNING: Data file not found: %s\n", DATA_PATH_2))
}

cat("\n\n=== Analysis Complete ===\n")
