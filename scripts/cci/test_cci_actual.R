# Actual test with data loading and function execution
# Run this in R session: source("scripts/cci/test_cci_actual.R")

cat("=== CCI Tool Actual Test ===\n\n")

# Load required packages
cat("1. Loading packages...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(qs)
})
cat("  ✓ Packages loaded\n\n")

# Load functions
cat("2. Loading CCI functions...\n")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
cat("  ✓ CCI functions loaded\n\n")

# Load run_nichenet_analysis
cat("3. Loading run_nichenet_analysis...\n")
if (file.exists("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")) {
  source("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")
  cat("  ✓ run_nichenet_analysis loaded\n\n")
} else {
  stop("CCI.R not found!")
}

# Load data
cat("4. Loading test data...\n")
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
cat("  ✓ Data loaded: ", ncol(sobj_test), " cells, ", nrow(sobj_test), " genes\n\n")

# Check metadata
cat("5. Checking metadata...\n")
colnames_meta <- colnames(sobj_test@meta.data)
cat("  - Metadata columns: ", length(colnames_meta), "\n")

if (!"anno3.scvi" %in% colnames_meta) {
  stop("anno3.scvi column not found!")
}
if (!"g3" %in% colnames_meta) {
  stop("g3 column not found!")
}

clusters <- unique(sobj_test@meta.data$anno3.scvi)
clusters <- clusters[!is.na(clusters)]
cat("  - Clusters: ", length(clusters), "\n")
cat("  - Cluster IDs: ", paste(head(clusters, 5), collapse = ", "), "\n")

g3_vals <- unique(sobj_test@meta.data$g3)
g3_vals <- g3_vals[!is.na(g3_vals)]
cat("  - g3 values: ", paste(g3_vals, collapse = ", "), "\n\n")

# Prepare DEG list
cat("6. Preparing DEG list...\n")
# Find a receiver cluster with sufficient cells in both conditions
receiver_cluster_test <- NULL
for (clust in clusters) {
  cells_in_clust <- sobj_test@meta.data[["anno3.scvi"]] == clust
  g3_in_clust <- sobj_test@meta.data[["g3"]][cells_in_clust]
  g3_in_clust <- g3_in_clust[!is.na(g3_in_clust)]
  unique_g3 <- unique(g3_in_clust)
  if (length(unique_g3) >= 2) {
    counts_per_cond <- table(g3_in_clust)
    if (min(counts_per_cond) >= 10) {  # At least 10 cells per condition
      receiver_cluster_test <- as.character(clust)
      cat("  - Selected receiver cluster: ", receiver_cluster_test, "\n")
      cat("    - Cells per condition: ", paste(names(counts_per_cond), "=", counts_per_cond, collapse = ", "), "\n")
      break
    }
  }
}

if (is.null(receiver_cluster_test)) {
  receiver_cluster_test <- as.character(clusters[1])
  cat("  - Using first cluster as receiver: ", receiver_cluster_test, "\n")
  cat("    (Warning: May not have sufficient cells in both conditions)\n")
}

# Use actual genes from the data that are likely to be in NicheNet
# Get genes that are expressed in the receiver
receiver_cells <- sobj_test@meta.data[["anno3.scvi"]] == receiver_cluster_test
receiver_expr <- Seurat::GetAssayData(sobj_test, assay = "RNA", slot = "counts")[, receiver_cells]
expressed_genes <- rownames(receiver_expr)[Matrix::rowSums(receiver_expr > 0) > 0]
actual_genes <- head(expressed_genes, 20)  # Use first 20 expressed genes

deg_df_example <- data.frame(
  gene = actual_genes,
  cluster = rep(receiver_cluster_test, length(actual_genes)),
  avg_log2FC = c(rep(1.5, 10), rep(-1.2, length(actual_genes) - 10)),
  p_val_adj = c(0.001, 0.0001, 0.01, 0.05, 0.1, rep(0.01, length(actual_genes) - 5)),
  stringsAsFactors = FALSE
)
cat("  ✓ DEG list prepared with ", nrow(deg_df_example), " genes\n\n")

# Test input validation
cat("7. Testing input validation...\n")
tryCatch({
  validation <- validate_cci_inputs(
    sobj = sobj_test,
    cluster_col = "anno3.scvi",
    deg_df = deg_df_example,
    receiver_cluster = receiver_cluster_test,
    sender_clusters = NULL
  )
  cat("  ✓ Input validation passed\n\n")
}, error = function(e) {
  cat("  ✗ Input validation failed: ", e$message, "\n\n")
  stop("Validation failed")
})

# Test DEG extraction
cat("8. Testing DEG extraction...\n")
tryCatch({
  receiver_degs <- extract_receiver_degs(
    deg_df_example,
    receiver_cluster_test,
    p_val_adj_cutoff = 0.05,
    logfc_cutoff = 0.25
  )
  cat("  ✓ DEG extraction passed: ", nrow(receiver_degs), " DEGs\n\n")
}, error = function(e) {
  cat("  ✗ DEG extraction failed: ", e$message, "\n\n")
  stop("DEG extraction failed")
})

# Test sender identification
cat("9. Testing sender identification...\n")
tryCatch({
  sender_clusters <- identify_sender_clusters(
    sobj_test,
    "anno3.scvi",
    receiver_cluster_test,
    NULL
  )
  cat("  ✓ Sender identification passed: ", length(sender_clusters), " senders\n\n")
}, error = function(e) {
  cat("  ✗ Sender identification failed: ", e$message, "\n\n")
  stop("Sender identification failed")
})

# Check if we can proceed with full analysis
cat("10. Checking if conditions are suitable for full analysis...\n")
receiver_cells <- sobj_test@meta.data$anno3.scvi == receiver_cluster_test
receiver_g3 <- sobj_test@meta.data$g3[receiver_cells]
receiver_g3 <- receiver_g3[!is.na(receiver_g3)]
unique_g3 <- unique(receiver_g3)

if (length(unique_g3) >= 2) {
  condition_oi <- as.character(unique_g3[1])
  condition_ref <- as.character(unique_g3[2])
  cat("  - Condition OI: ", condition_oi, "\n")
  cat("  - Condition Ref: ", condition_ref, "\n")
  cat("  ✓ Conditions suitable for analysis\n\n")
  
  cat("11. Running full CCI analysis (this may take several minutes)...\n")
  cat("    Note: This will download NicheNet data if not already cached\n")
  cat("    Note: Network issues may cause download failures - this is normal\n\n")
  
  # Check if NicheNet data already exists
  nichenet_dir <- file.path(getwd(), "nichenet_data_human")
  if (dir.exists(nichenet_dir)) {
    files <- list.files(nichenet_dir, pattern = "\\.rds$")
    cat("    Found ", length(files), " NicheNet data files in cache\n\n")
  }
  
  # Use existing NicheNet data directory
  nichenet_data_dir <- "/data/user3/git_repo/human"
  cat("    Using NicheNet data from: ", nichenet_data_dir, "\n\n")
  
  # Check if SCT assay needs PrepSCTFindMarkers
  cat("    Checking Seurat assay setup...\n")
  if ("SCT" %in% names(sobj_test@assays)) {
    cat("    - SCT assay found, may need PrepSCTFindMarkers\n")
    # Try to run PrepSCTFindMarkers if needed
    tryCatch({
      sobj_test <- Seurat::PrepSCTFindMarkers(sobj_test)
      cat("    - PrepSCTFindMarkers completed\n")
    }, error = function(e) {
      cat("    - PrepSCTFindMarkers warning: ", e$message, "\n")
      cat("    - Will try with RNA assay instead\n")
    })
  }
  
  tryCatch({
    results <- run_cci_analysis(
      sobj = sobj_test,
      cluster_col = "anno3.scvi",
      deg_df = deg_df_example,
      receiver_cluster = receiver_cluster_test,
      condition_col = "g3",
      condition_oi = condition_oi,
      condition_ref = condition_ref,
      species = "human",
      assay_name = "RNA",  # Use RNA assay to avoid SCT issues
      verbose = TRUE,
      auto_save = TRUE,
      top_n_ligands = 10,  # Reduced for faster testing
      run_circos = FALSE,  # Skip circos for faster testing
      nichenet_data_dir = nichenet_data_dir,  # Use existing data directory
      p_val_adj_cutoff = 1.1,  # Very lenient to ensure DEGs pass filter (p_val_adj can be > 1.0)
      logfc_cutoff = 0.05  # Very lenient logfc threshold
    )
    
    cat("\n  ✓ Full CCI analysis completed!\n")
    cat("  - Top ligands: ", length(results$nichenet_results$best_upstream_ligands), "\n")
    cat("  - Sender clusters: ", length(results$sender_clusters), "\n")
    if (!is.null(results$saved_path)) {
      cat("  - Results saved to: ", results$saved_path, "\n")
    }
    
    # Show top ligands
    if (length(results$nichenet_results$best_upstream_ligands) > 0) {
      cat("  - Top 5 ligands: ", paste(head(results$nichenet_results$best_upstream_ligands, 5), collapse = ", "), "\n")
    }
    
  }, error = function(e) {
    cat("  ✗ Full analysis failed: ", e$message, "\n")
    if (grepl("download", e$message, ignore.case = TRUE)) {
      cat("    This is a network/download issue. NicheNet data download may have failed.\n")
      cat("    You can try again later or manually download the data.\n")
    } else {
      cat("    Check the error message above for details.\n")
    }
  })
  
} else {
  cat("  ✗ Not enough conditions in receiver cluster for analysis\n")
  cat("    Unique g3 values in receiver: ", paste(unique_g3, collapse = ", "), "\n")
  cat("    Need at least 2 different conditions\n\n")
}

cat("\n=== Test Complete ===\n")

