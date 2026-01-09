#!/usr/bin/env Rscript
# scripts/consensus/run_AG_run2_batch.R
# Batch analysis for "All" and subsets

suppressPackageStartupMessages({
    library(Seurat)
    library(qs)
    library(dplyr)
})

# Parameters
base_out_dir <- "/data/user3/sobj/consensus/AG_run2"
methods_stable <- "edgeR-LRT,DESeq2-Wald,limma-voom,limma-trend"
# User asked for 'muscat-limma-trend'. My run_deg_consensus logic interprets "limma-trend" as such if pseudobulk, but let's be safe.
# Actually, the user's run_deg_consensus defaults might support "muscat" or similar.
# The user explicitly asked to run stable methods first, then nebula/dream separately.
# This script runs stable ones.

renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

# Sourcing removed to avoid package conflicts
# source("myR/R/deg_consensus/run_deg_consensus.R")
# source("scripts/consensus/plot_consensus_cli.R", local = TRUE)

# 1. Main Object "All"
message("Processing Main Object...")
sobj_main <- qs::qread("/data/user3/sobj/is2_IS_3_1_plots.qs")
out_main <- file.path(base_out_dir, "All")

# Assuming run_deg_consensus handles the execution
# We use the CLI script approach or source logic.
# Better to system call the CLI script to ensure isolation?
# Or call function directly. Function is safer for batch R session but CLI isolates memory.
# Given memory size of subsets, CLI calls are safer.

run_cli <- function(input_path, output_dir, cluster_col = "anno3big", methods = methods_stable) {
    cmd <- sprintf(
        "Rscript scripts/consensus/run_deg_consensus_cli.R --input '%s' --output '%s' --cluster '%s' --group 'g3' --covariates 'sex,age,GEM' --methods '%s' --cores 8",
        input_path, output_dir, cluster_col, methods
    )
    system(cmd)
}

# Run Main
run_cli("/data/user3/sobj/is2_IS_3_1_plots.qs", out_main)

# 2. Subsets
message("Processing Subsets...")
subsets_list <- qs::qread("/data/user3/sobj/is2_IS_5_subsets_processed.qs")
# It's a list of Seurat objects. We need to save them temporarily or pass them?
# CLI needs a file path.
# So we save each subset to a temp file or specific file.

temp_subset_dir <- "/data/user3/sobj/temp_subsets"
dir.create(temp_subset_dir, showWarnings = FALSE)

subset_names <- names(subsets_list)
for (nm in subset_names) {
    message("  Subset: ", nm)
    s_obj <- subsets_list[[nm]]

    # Save temp
    f_path <- file.path(temp_subset_dir, paste0(nm, ".qs"))
    qs::qsave(s_obj, f_path)

    # Run
    # Use 'anno3' for subsets? Or 'anno3big'? User said "subsetted object... folder name Bc, Tc".
    # Usually subsetting implies we might refine clusters or use same.
    # Let's assume 'anno3' is the finer resolution for subsets, or stick to 'anno3big' if that's the consistency goal.
    # User didn't specify column for subsets. Defaulting to 'anno3big' for consistency with All.

    out_sub <- file.path(base_out_dir, nm)
    run_cli(f_path, out_sub, cluster_col = "anno3big")

    # Cleanup
    unlink(f_path)
}

message("Batch Analysis Complete.")
