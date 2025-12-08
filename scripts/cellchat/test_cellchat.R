
# scripts/cellchat/test_cellchat.R

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(qs)
  library(future)
})

# Source wrapper
# Assuming running from project root
if (file.exists("myR/R/cci_cellchat_wrapper.R")) {
  source("myR/R/cci_cellchat_wrapper.R")
} else {
  stop("Wrapper file not found. Please run from project root.")
}

# Load data
data_path <- "/data/user3/sobj/IS6_sex_added_0.1x_251110.qs"
if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

message("Loading data from ", data_path)
sobj <- qs::qread(data_path)

# Subset for quick testing
message("Subsetting to 500 cells for quick testing...")
set.seed(123)
sobj <- sobj[, sample(colnames(sobj), 500)]

# Check metadata
message("Metadata columns: ", paste(colnames(sobj@meta.data), collapse = ", "))
# User said cluster info is in $anno3.scvi
group_col <- "anno3.scvi"

if (!group_col %in% colnames(sobj@meta.data)) {
  stop("Group column '", group_col, "' not found in metadata.")
}

# Run analysis
output_dir <- "docs/cellchat/test_run"
message("Running CellChat analysis...")

# Use a subset of DB for faster testing if needed, or full DB
# db.use = "Secreted Signaling" 
# For now, let's try full DB but maybe limit cores to avoid OOM if machine is small
# User said "renv installed", so environment should be fine.

cellchat <- run_cellchat_analysis(
  sobj = sobj,
  group.by = group_col,
  species = "human", 
  output_dir = output_dir,
  n_cores = 4,
  verbose = TRUE
)

message("Test completed successfully!")
