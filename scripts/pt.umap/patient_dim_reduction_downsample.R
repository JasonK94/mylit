#!/usr/bin/env Rscript

repo_root <- normalizePath(Sys.getenv("MYR_REPO_ROOT", "/home/user3/data_user3/git_repo/_wt/pt.umap"), mustWork = TRUE)
st_dir <- normalizePath(file.path(repo_root, "st"), mustWork = TRUE)
setwd(st_dir)

sample_args <- commandArgs(trailingOnly = TRUE)
sample_size <- if (length(sample_args) >= 1) as.integer(sample_args[[1]]) else 4000
color_by    <- if (length(sample_args) >= 2) sample_args[[2]] else "g3"
label_flag  <- if (length(sample_args) >= 3) as.logical(sample_args[[3]]) else FALSE

message("▶ Running patient-level reduction test in ", st_dir)
message("  • Sample size: ", sample_size)
message("  • Colour by: ", color_by)
message("  • Labels: ", label_flag)

source("start.R")
setwd(st_dir)

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools is required.")
devtools::load_all(file.path(repo_root, "myR"), quiet = TRUE)

library(qs); library(Seurat); library(dplyr); library(ggplot2)

sobj_path <- "sobj/IS5_g3NA_removal_251110.qs"
if (!file.exists(sobj_path)) stop("Seurat object not found at ", sobj_path)

seurat_obj <- qs::qread(sobj_path)
message("  • Loaded cells: ", ncol(seurat_obj))

set.seed(314)
if (sample_size < ncol(seurat_obj)) {
  sel_cells <- sample(colnames(seurat_obj), sample_size)
  seurat_obj <- subset(seurat_obj, cells = sel_cells)
  message("  • Downsampled cells: ", ncol(seurat_obj))
}

markers_df <- FindAllMarkers(
  seurat_obj,
  group.by = "anno3.scvi",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)
message("  • Markers computed: ", nrow(markers_df))

res <- patient_dimensionality_reduction(
  seurat_obj = seurat_obj,
  markers_df = markers_df,
  reduction = "integrated.scvi",
  verbose   = TRUE
)

message("  • Combined feature matrix dims: ", paste(dim(res$combined_features), collapse = " x "))
message("  • Embedding dims: ", paste(dim(res$embedding), collapse = " x "))

output_dir <- file.path(st_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

qs::qsave(res, file.path(output_dir, sprintf("patient_dim_reduction_downsample_%s.qs", sample_size)))
write.csv(res$plot_df, file.path(output_dir, sprintf("patient_dim_reduction_downsample_%s_%s.csv", sample_size, color_by)), row.names = FALSE)

p <- plot_patient_umap(res$plot_df, color_by = color_by, label = label_flag)
ggplot2::ggsave(file.path(output_dir, sprintf("patient_umap_%s_%s.png", color_by, sample_size)), p, width = 6, height = 5, dpi = 300)

message("✅ Test run complete.")