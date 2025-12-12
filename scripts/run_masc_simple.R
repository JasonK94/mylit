#!/usr/bin/env Rscript
# Simple MASC Analysis for Stroke Data (Minimal Variables)

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
DATA_PATH <- "/data/user3/sobj/is2_IS_3_clustered.qs"
OUTPUT_DIR <- "/data/user3/sobj/masc/stroke_simple"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Simple MASC Analysis ===\n")

# Load Data
cat(sprintf("Loading data from %s...\n", DATA_PATH))
seurat_obj <- qs::qread(DATA_PATH)

# --- Preprocessing ---
# 1. Target Variable: g3
if (!"g3" %in% colnames(seurat_obj@meta.data)) stop("g3 not found")
seurat_obj@meta.data$g3 <- as.factor(seurat_obj@meta.data$g3)
cat("Target variable (g3) levels:", paste(levels(seurat_obj@meta.data$g3), collapse=", "), "\n")

# 2. Random Effect: hos_no (Ensure character/factor)
if (!"hos_no" %in% colnames(seurat_obj@meta.data)) stop("hos_no not found")
# Explicitly convert to character first to avoid numeric interpretation
seurat_obj@meta.data$hos_no <- as.character(seurat_obj@meta.data$hos_no)
cat("Random effect (hos_no) unique values:", length(unique(seurat_obj@meta.data$hos_no)), "\n")

# 3. Cluster Variable: anno3big
cluster_var <- "anno3big"
if (!cluster_var %in% colnames(seurat_obj@meta.data)) stop("Cluster variable not found")
cat("Cluster variable:", cluster_var, "\n")

# 4. Optional Batch: GEM
use_batch <- "GEM" %in% colnames(seurat_obj@meta.data)
if (use_batch) {
    seurat_obj@meta.data$GEM <- as.factor(seurat_obj@meta.data$GEM)
    cat("Batch variable (GEM) used.\n")
}

# Define Model Variables
# Only using GEM as fixed effect if available, and hos_no as random effect
fixed_effects <- if(use_batch) "GEM" else NULL
random_effects <- "hos_no"

cat(sprintf("Fixed effects: %s\n", paste(fixed_effects, collapse=", ")))
cat(sprintf("Random effects: %s\n", paste(random_effects, collapse=", ")))

cat("\nRunning MASC pipeline...\n")
results <- run_masc_pipeline(
    seurat_obj = seurat_obj,
    cluster_var = cluster_var,
    contrast_var = "g3",
    random_effects = random_effects,
    fixed_effects = fixed_effects,
    save = TRUE,
    output_dir = OUTPUT_DIR,
    prefix = "masc_simple",
    force_run = TRUE,
    plotting = TRUE,
    save_models = FALSE,
    adjust_pvalue = TRUE,
    verbose = TRUE
)

cat("\n=== Final Results ===\n")
print(results$masc_results)

# Save summary
sink(file.path(OUTPUT_DIR, "masc_simple_summary.txt"))
print(results$masc_results)
sink()


