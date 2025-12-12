#!/usr/bin/env Rscript
# MultiNicheNet Analysis for IS2/IS3 Dataset - Permissive/Relaxed Mode
# Designed to recover interactions for cell types with subtle changes

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(qs)
})

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║  MultiNicheNet Analysis (Permissive Mode) - IS2/IS3       ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# Configuration with Relaxed Thresholds
CONFIG <- list(
    data_path = "/data/user3/sobj/is2_IS_3_1_plots.qs",
    sample_id = "hos_no",
    group_id = "g3",
    celltype_id = "anno3",
    contrasts = "X2-X1",
    output_dir = "multinichenet_IS2_IS3_permissive", # Separate output dir

    # RELAXED PARAMETERS
    min_cells = 5, # Lowered from 10 to include more samples
    logfc_threshold = 0.05, # Lowered from 0.10 to catch subtle shifts
    fraction_cutoff = 0.05, # Keep at 5%
    p_val_threshold = 0.20, # Relaxed from 0.05/0.10 to 0.20 (Raw p-val)
    p_val_adj = FALSE, # Strictly use raw p-value

    cores = 12, # Use more cores if available
    species = "human",
    verbose = TRUE
)

cat("Relaxed Configuration:\n")
cat("  p_val_threshold:   ", CONFIG$p_val_threshold, " (Raw P)\n")
cat("  logfc_threshold:   ", CONFIG$logfc_threshold, "\n")
cat("  min_cells:         ", CONFIG$min_cells, "\n")
cat("  output_dir:        ", CONFIG$output_dir, "\n\n")

# Source wrapper
wrapper_path <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci_multinichenet_wrapper.R"
source(wrapper_path)

# Load data
sobj <- qs::qread(CONFIG$data_path)

# Filter NAs
sobj_clean <- sobj[, !is.na(sobj@meta.data[[CONFIG$group_id]])]
sobj_clean <- sobj_clean[, !is.na(sobj_clean@meta.data[[CONFIG$sample_id]])]

# Run Analysis
run_multinichenet_analysis(
    sobj = sobj_clean,
    sample_id = CONFIG$sample_id,
    group_id = CONFIG$group_id,
    celltype_id = CONFIG$celltype_id,
    contrasts_oi = CONFIG$contrasts,
    senders_oi = NULL,
    receivers_oi = NULL,
    min_cells = CONFIG$min_cells,
    species = CONFIG$species,
    output_dir = CONFIG$output_dir,
    verbose = CONFIG$verbose,
    cores = CONFIG$cores,

    # Overriding defaults in wrapper if necessary,
    # but the wrapper function signature needs to support passing these specific threshold args
    # if they are hardcoded in the wrapper's call to multinichenetr.
    # Let's check wrapper content.
    # The wrapper in previous turn (Step 54) has:
    # fraction_cutoff = 0.05,
    # logFC_threshold = 0.10,
    # p_val_threshold = 0.05,
    # p_val_adj = FALSE
    # Hardcoded in the call to multinichenetr::multi_nichenet_analysis inside the wrapper!

    # Relaxed thresholds
    min_pct = CONFIG$fraction_cutoff,
    logfc_thresh = CONFIG$logfc_threshold,
    p_val_thresh = CONFIG$p_val_threshold,
    p_val_adj = CONFIG$p_val_adj
)
