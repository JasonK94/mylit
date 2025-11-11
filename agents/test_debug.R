# Test script to debug the functions
# Load required packages and test functions

# Source the functions
source("/data/kjc1/mylit/myR/R/test_cursor.R")

# Check if data file exists and load it
data_file <- "/data/kjc1/mylit/projects/stroke/final_downsampled_251105.qs"
if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

# Try loading the data
if (requireNamespace("qs", quietly = TRUE)) {
  data_test <- qs::qread(data_file)
} else {
  # Try RDS format
  data_test <- readRDS(sub("\\.qs$", ".rds", data_file))
}

# Check if it's a Seurat object
if (!inherits(data_test, "Seurat")) {
  stop("data_test is not a Seurat object")
}

# Test the data
message("Data loaded successfully")
message("Number of cells: ", ncol(data_test))
message("Number of genes: ", nrow(data_test))
message("Columns in meta.data: ", paste(head(colnames(data_test@meta.data), 10), collapse=", "))

# Modify type column
if ("type" %in% colnames(data_test@meta.data)) {
  data_test$type <- make.names(data_test$type)
  message("Type column modified with make.names")
  message("Unique types: ", paste(unique(data_test$type), collapse=", "))
} else {
  stop("'type' column not found in meta.data")
}

# Check required columns
required_cols <- c("seurat_clusters", "hos_no", "type")
missing_cols <- setdiff(required_cols, colnames(data_test@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse=", "))
}

# Test runMUSCAT
message("\n=== Testing runMUSCAT ===")
tryCatch({
  musc2 <- runMUSCAT(
    data_test, 
    cluster_id = "seurat_clusters", 
    sample_id = "hos_no",
    group_id = "type",
    formula_str = "~ 0 + group", 
    contrast = "groupICH-groupIS", 
    method = "edgeR"
  )
  message("SUCCESS: runMUSCAT completed")
  print(str(musc2))
}, error = function(e) {
  message("ERROR in runMUSCAT: ", e$message)
  traceback()
})

