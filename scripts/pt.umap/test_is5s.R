#!/usr/bin/env Rscript

# Test script for is5s (downsampled stroke PBMC data)
# Data location: /data/user3/sobj/IS6_sex_added_0.1x_251110.qs

# Get repo root and normalize BEFORE any directory changes
# Save to a different variable name to avoid being overwritten by start.R
pt_umap_repo_root <- Sys.getenv("MYR_REPO_ROOT", "/home/user3/data_user3/git_repo/_wt/pt.umap")
pt_umap_repo_root <- normalizePath(pt_umap_repo_root, mustWork = TRUE)
myr_path <- paste0(pt_umap_repo_root, "/myR")
if (!dir.exists(myr_path)) stop("myR directory not found at ", myr_path)

st_dir <- file.path(pt_umap_repo_root, "st")
st_dir <- normalizePath(st_dir, mustWork = TRUE)
setwd(st_dir)

sample_args <- commandArgs(trailingOnly = TRUE)
color_by    <- if (length(sample_args) >= 1) sample_args[[1]] else "g3"
label_flag  <- if (length(sample_args) >= 2) as.logical(sample_args[[2]]) else FALSE

message("▶ Running patient-level reduction test for is5s (downsampled data)")
message("  • Data: /data/user3/sobj/IS6_sex_added_0.1x_251110.qs")
message("  • Colour by: ", color_by)
message("  • Labels: ", label_flag)

source("start.R")
setwd(st_dir)

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools is required.")
# myr_path was saved before start.R could overwrite repo_root
devtools::load_all(myr_path, quiet = TRUE)

library(qs); library(Seurat); library(dplyr); library(ggplot2)

# Load is5s data (downsampled)
sobj_path <- "/data/user3/sobj/IS6_sex_added_0.1x_251110.qs"
if (!file.exists(sobj_path)) stop("Seurat object not found at ", sobj_path)

message("  • Loading Seurat object from ", sobj_path)
seurat_obj <- qs::qread(sobj_path)
message("  • Loaded cells: ", ncol(seurat_obj))
message("  • Patients (hos_no): ", length(unique(seurat_obj@meta.data$hos_no)))
message("  • Clusters (anno3.scvi): ", length(unique(seurat_obj@meta.data$anno3.scvi)))

# Compute markers
message("  • Computing markers...")
markers_df <- FindAllMarkers(
  seurat_obj,
  group.by = "anno3.scvi",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)
message("  • Markers computed: ", nrow(markers_df))

# Run patient dimensionality reduction
message("  • Running patient dimensionality reduction pipeline...")
res <- patient_dimensionality_reduction(
  seurat_obj = seurat_obj,
  markers_df = markers_df,
  reduction = "integrated.scvi",
  verbose   = TRUE
)

message("  • Combined feature matrix dims: ", paste(dim(res$combined_features), collapse = " x "))
message("  • Embedding dims: ", paste(dim(res$embedding), collapse = " x "))
message("  • Number of patients in embedding: ", nrow(res$embedding))

# Save results
output_dir <- file.path(st_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_prefix <- "patient_dim_reduction_is5s"
qs::qsave(res, file.path(output_dir, paste0(output_prefix, ".qs")))
write.csv(res$plot_df, file.path(output_dir, paste0(output_prefix, "_", color_by, ".csv")), row.names = FALSE)

# Plot
p <- plot_patient_umap(res$plot_df, color_by = color_by, label = label_flag)
ggplot2::ggsave(
  file.path(output_dir, paste0("patient_umap_is5s_", color_by, ".png")), 
  p, 
  width = 6, 
  height = 5, 
  dpi = 300
)

message("✅ Test run complete for is5s.")
message("  • Results saved to: ", output_dir)

