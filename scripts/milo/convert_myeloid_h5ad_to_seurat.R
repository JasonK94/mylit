#!/usr/bin/env Rscript
# Script to convert myeloid_v1_TNBC.h5ad file to Seurat object and save as .qs file

# Load required libraries
library(reticulate)
library(Seurat)
library(qs)

# Source the conversion function from main repo
# Use absolute path
main_repo_r <- "/data/kjc2/git_repo/mylit/myR/R/h5ad_to_seurat.R"
if (file.exists(main_repo_r)) {
  source(main_repo_r)
} else {
  stop("Could not find h5ad_to_seurat.R function at: ", main_repo_r)
}

# File paths
h5ad_path <- "/data/ARPAH/250305_TNBC_26/Xenium_analysis/transfer_to_python/TNBC_only/myeloid_v1_TNBC.h5ad"
# Save in current worktree directory
output_dir <- "/data/kjc2/git_repo/_wt/h5ad2sobj"
output_name <- "myeloid_v1_TNBC"

# Check if file exists
if (!file.exists(h5ad_path)) {
  stop("h5ad file not found: ", h5ad_path)
}

# Convert h5ad to Seurat and save as .qs
message("=", rep("=", 79))
message("Converting h5ad to Seurat format...")
message("Input: ", h5ad_path)
message("Output directory: ", output_dir)
message("=", rep("=", 79))

seurat_obj <- load_h5ad_to_seurat_qs(
  h5ad_path = h5ad_path,
  condaenv = "/home/jaecheon/miniconda3/envs/scenvi",
  save_path = output_dir,
  save_name = output_name,
  use_layer = "counts",  # Use counts layer if available
  project_name = "TNBC_Myeloid",
  # 필요한 메타데이터 컬럼만 선택 (파일 크기 감소)
  required_meta_cols = c("batch", "Core", "intra_T", "TIL_2", "Annot3"),
  # 필요한 reduction만 선택 (X_scVI만 사용)
  required_obsm_keys = c("X_scVI", "X_umap")  # integrated.scvi와 umap만
)

message("\n", "=", rep("=", 79))
message("Conversion complete!")
message("Seurat object saved to: ", file.path(output_dir, paste0(output_name, ".qs")))

# Verify the conversion
message("\n", "=", rep("=", 79))
message("Verifying Seurat object...")
message("Dimensions: ", nrow(seurat_obj), " genes x ", ncol(seurat_obj), " cells")
message("Available reductions: ", paste(Seurat::Reductions(seurat_obj), collapse = ", "))

# Check required metadata columns
required_cols <- c("batch", "Core", "intra_T", "TIL_2", "Annot3")
message("\nMetadata columns check:")
for (col in required_cols) {
  if (col %in% colnames(seurat_obj@meta.data)) {
    message("  ✓ ", col, ": Found")
    message("    Unique values: ", paste(unique(seurat_obj@meta.data[[col]]), collapse = ", "))
  } else {
    message("  ✗ ", col, ": MISSING")
  }
}

message("\n", "=", rep("=", 79))
message("All done!")

