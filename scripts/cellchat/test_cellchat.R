# scripts/cellchat/test_cellchat.R

# Set library path to shared renv
renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) {
  .libPaths(c(renv_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(qs)
  library(future)
})

# Source wrapper
wrapper_path <- "myR/R/cci_cellchat_wrapper.R"
if (!file.exists(wrapper_path)) {
  wrapper_path <- "/home/user3/data_user3/git_repo/_wt/cellchat/myR/R/cci_cellchat_wrapper.R"
}
source(wrapper_path)

# Data path
data_path <- "/data/user3/sobj/IS6_sex_added_0.1x_251110.qs"
group_col <- "anno3.scvi"

# Run analysis using the wrapper with path input
# This tests:
# 1. Path input support
# 2. Logging (check logs/cc/runX)
# 3. Checkpointing (will create output_dir/checkpoints)
# 4. qs usage

output_dir <- "docs/cellchat/test_run_cli"
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)

message("Running CellChat analysis via wrapper...")

cellchat <- run_cellchat_analysis(
  input_data = data_path,
  group.by = group_col,
  species = "human",
  output_dir = output_dir,
  n_cores = 4, # Use fewer cores for testing
  verbose = TRUE
)

message("Test completed successfully!")
