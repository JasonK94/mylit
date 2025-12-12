#!/usr/bin/env Rscript
# Test MultiNicheNet analysis

# Set environment
# Set environment
# Note: We rely on renv being active in the project root via .Rprofile
# if (file.exists("renv/activate.R")) {
#   source("renv/activate.R")
# }

library(Seurat)
library(dplyr)
library(qs)
library(multinichenetr)

# Source functions using robust path
get_script_dir <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- "--file="
    match <- grep(file_arg, cmd_args)
    if (length(match) > 0) {
        return(dirname(normalizePath(sub(file_arg, "", cmd_args[match]))))
    } else {
        return(getwd())
    }
}

script_dir <- get_script_dir()
# Assuming script is in myR/tests/
# Wrapper is in myR/R/cci_multinichenet_wrapper.R
wrapper_path <- file.path(script_dir, "../R/cci_multinichenet_wrapper.R")

if (file.exists(wrapper_path)) {
    source(wrapper_path)
} else {
    # Fallback for interactive run or if script_dir is not correct
    # Try absolute path
    abs_path <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci_multinichenet_wrapper.R"
    if (file.exists(abs_path)) {
        source(abs_path)
    } else {
        stop("Could not find cci_multinichenet_wrapper.R")
    }
}

cat("=== MultiNicheNet Test ===\n\n")

# 1. Load Data
# Use the downsampled data if available, otherwise load full and downsample
ds_path <- "/data/user3/sobj/IS5_ds10_per_patient_cluster.qs"
if (file.exists(ds_path)) {
    cat("Loading downsampled data from:", ds_path, "\n")
    sobj <- qs::qread(ds_path)
} else {
    cat("Downsampled data not found. Loading full data...\n")
    path <- "/data/user3/sobj/IS6_sex_added_251110.qs"
    if (!file.exists(path)) stop("Data not found:", path)
    sobj <- qs::qread(path)

    # Downsample for testing
    cat("Downsampling...\n")
    sobj <- subset(sobj, downsample = 100)
}

# 2. Prepare Metadata
cat("Preparing metadata...\n")
meta_cols <- colnames(sobj@meta.data)

# Identify columns
sample_col <- if ("hos_no" %in% meta_cols) "hos_no" else "sample"
group_col <- "g3"
celltype_col <- "anno3.scvi"

if (!all(c(sample_col, group_col, celltype_col) %in% meta_cols)) {
    stop("Missing metadata columns.")
}

# Ensure group is valid
sobj <- sobj[, !is.na(sobj@meta.data[[group_col]])]
# Ensure sample is valid
sobj <- sobj[, !is.na(sobj@meta.data[[sample_col]])]

# Create a unique sample_id (patient + group) if needed, but multinichenetr handles sample_id as patient usually
# If sample_id is patient, and patients are in one group, it's fine.
# If patients are paired, it's also fine.

# 3. Run Analysis
cat("Running MultiNicheNet analysis...\n")
results <- run_multinichenet_analysis(
    sobj = sobj,
    sample_id = sample_col,
    group_id = group_col,
    celltype_id = celltype_col,
    contrasts_oi = "X2-X1", # make.names converts 1->X1, 2->X2
    # multinichenetr::get_DE_info uses contrasts. If group_id is factor, contrasts should be valid.
    # Let's check levels
    min_cells = 5,
    species = "human",
    output_dir = "test_multinichenet_output",
    verbose = TRUE,
    # Restrict to specific cell types for faster testing
    receivers_oi = unique(sobj@meta.data[[celltype_col]])[1],
    senders_oi = unique(sobj@meta.data[[celltype_col]])[1:3]
)

cat("Analysis complete.\n")

# 4. Test Plotting
cat("Testing Circos plot...\n")
# Plot for group "2"
plot_multinichenet_circos(results, group_oi = "X2", output_file = "test_multinichenet_output/circos_g2.pdf")

cat("Done.\n")
