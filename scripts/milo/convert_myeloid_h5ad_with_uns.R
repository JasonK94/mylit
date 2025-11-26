#!/usr/bin/env Rscript
# Experimental script to convert h5ad file while retaining uns/obsp/var structures

library(reticulate)
library(Seurat)
library(qs)

main_repo_r <- "/home/jaecheon/kjc2/git_repo/mylit/myR/R/h5ad_to_seurat.R"
if (file.exists(main_repo_r)) {
  source(main_repo_r)
} else {
  stop("Could not find h5ad_to_seurat.R at: ", main_repo_r)
}

h5ad_path <- "/data/ARPAH/250305_TNBC_26/Xenium_analysis/transfer_to_python/TNBC_only/myeloid_v1_TNBC.h5ad"
output_dir <- "/data/kjc2/git_repo/_wt/h5ad2sobj"
output_name <- "myeloid_v1_TNBC_full"

if (!file.exists(h5ad_path)) {
  stop("h5ad file not found: ", h5ad_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("=", rep("=", 79))
message("Experimental conversion retaining uns/obsp/var")
message("Input: ", h5ad_path)
message("Output directory: ", output_dir)
message("=", rep("=", 79))

full_seurat <- load_h5ad_to_seurat_qs(
  h5ad_path = h5ad_path,
  condaenv = "/home/jaecheon/miniconda3/envs/scenvi",
  save_path = output_dir,
  save_name = output_name,
  use_layer = "counts",
  project_name = "TNBC_Myeloid_full",
  required_meta_cols = c("batch", "Core", "intra_T", "TIL_2", "Annot3"),
  required_obsm_keys = NULL,
  store_uns = TRUE,
  store_obsp = TRUE,
  store_var = TRUE
)

message("=", rep("=", 79))
message("Full conversion complete: ", file.path(output_dir, paste0(output_name, ".qs")))
message("Object dims: ", nrow(full_seurat), " genes x ", ncol(full_seurat), " cells")
message("Reductions: ", paste(Seurat::Reductions(full_seurat), collapse = ", "))
message("Misc entries: ", paste(names(full_seurat@misc), collapse = ", "))
message("=", rep("=", 79))