#' Pseudobulk Differential Expression Analysis
#'
#' This module provides functions for performing pseudo-bulk differential
#' gene expression analysis using edgeR.
#'
#' @name pseudobulk_deg
NULL

#' Prepare Pseudobulk Data for edgeR
#'
#' Aggregates single-cell counts to pseudo-bulk level for differential expression analysis.
#'
#' @param object Seurat object
#' @param sample_col Metadata column identifying samples/replicates
#' @param group_col Metadata column for grouping (e.g., condition, treatment)
#' @param cluster_col Optional cluster column for per-cluster analysis (default: NULL)
#' @param target_cluster If cluster_col is specified, analyze only this cluster (default: NULL)
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#' @param min_cells Minimum cells per sample to include (default: 10)
#' @param min_counts Minimum total counts per gene (default: 10)
#'
#' @return List containing:
#'   \item{counts}{Aggregated count matrix}
#'   \item{metadata}{Sample-level metadata}
#'   \item{cluster}{Cluster being analyzed (if applicable)}
#'
#' @examples
#' \dontrun{
#' # Overall pseudobulk
#' pb_data <- prepare_pseudobulk_edgeR(seurat_obj, "sample", "condition")
#' 
#' # Per-cluster pseudobulk
#' pb_data <- prepare_pseudobulk_edgeR(
#'   seurat_obj, 
#'   "sample", 
#'   "condition",
#'   cluster_col = "seurat_clusters",
#'   target_cluster = "0"
#' )
#' }
#'
#' @export
prepare_pseudobulk_edgeR <- function(object, 
                                     sample_col, 
                                     group_col,
                                     cluster_col = NULL,
                                     target_cluster = NULL,
                                     assay = "RNA",
                                     slot = "counts",
                                     min_cells = 10,
                                     min_counts = 10) {
  
  # Validate inputs
  required_cols <- c(sample_col, group_col)
  if (!is.null(cluster_col)) {
    required_cols <- c(required_cols, cluster_col)
  }
  
  missing_cols <- setdiff(required_cols, colnames(object@meta.data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  # Subset to target cluster if specified
  if (!is.null(cluster_col) && !is.null(target_cluster)) {
    cells_to_use <- colnames(object)[object@meta.data[[cluster_col]] == target_cluster]
    
    if (length(cells_to_use) == 0) {
      stop("No cells found for cluster: ", target_cluster)
    }
    
    object <- object[, cells_to_use]
    message("Analyzing cluster ", target_cluster, " (", length(cells_to_use), " cells)")
  }
  
  # Get counts
  counts <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  
  # Get metadata
  metadata <- object@meta.data[, required_cols, drop = FALSE]
  
  # Create sample identifiers (unique per sample-cluster combination if applicable)
  if (!is.null(cluster_col) && is.null(target_cluster)) {
    metadata$pb_sample <- paste0(metadata[[sample_col]], "_", metadata[[cluster_col]])
  } else {
    metadata$pb_sample <- as.character(metadata[[sample_col]])
  }
  
  # Check for minimum cells per sample
  cells_per_sample <- table(metadata$pb_sample)
  valid_samples <- names(cells_per_sample)[cells_per_sample >= min_cells]
  
  if (length(valid_samples) == 0) {
    stop("No samples have >= ", min_cells, " cells")
  }
  
  if (length(valid_samples) < length(cells_per_sample)) {
    n_removed <- length(cells_per_sample) - length(valid_samples)
    message("Removing ", n_removed, " samples with < ", min_cells, " cells")
    
    keep_cells <- metadata$pb_sample %in% valid_samples
    counts <- counts[, keep_cells, drop = FALSE]
    metadata <- metadata[keep_cells, , drop = FALSE]
  }
  
  # Aggregate counts by sample
  unique_samples <- unique(metadata$pb_sample)
  pb_counts <- matrix(
    0, 
    nrow = nrow(counts), 
    ncol = length(unique_samples),
    dimnames = list(rownames(counts), unique_samples)
  )
  
  for (samp in unique_samples) {
    cells_in_sample <- metadata$pb_sample == samp
    pb_counts[, samp] <- Matrix::rowSums(counts[, cells_in_sample, drop = FALSE])
  }
  
  # Filter lowly expressed genes
  keep_genes <- Matrix::rowSums(pb_counts) >= min_counts
  pb_counts <- pb_counts[keep_genes, , drop = FALSE]
  
  message("Pseudobulk matrix: ", nrow(pb_counts), " genes x ", ncol(pb_counts), " samples")
  
  # Create sample-level metadata
  sample_metadata <- metadata[!duplicated(metadata$pb_sample), , drop = FALSE]
  rownames(sample_metadata) <- sample_metadata$pb_sample
  sample_metadata <- sample_metadata[colnames(pb_counts), , drop = FALSE]
  
  return(list(
    counts = pb_counts,
    metadata = sample_metadata,
    cluster = target_cluster
  ))
}

#' Run Pseudobulk Differential Expression with edgeR
#'
#' Performs edgeR-based pseudo-bulk differential gene expression analysis.
#'
#' @param object Seurat object or prepared pseudobulk list from prepare_pseudobulk_edgeR
#' @param sample_col Metadata column identifying samples (required if object is Seurat)
#' @param group_col Metadata column for comparison (required if object is Seurat)
#' @param comparison Comparison to make (e.g., c("Treatment", "Control"))
#' @param cluster_col Optional cluster column (default: NULL)
#' @param target_cluster Specific cluster to analyze (default: NULL for overall)
#' @param mode Analysis mode: "overall", "per_cluster", "specific_cluster" (default: "overall")
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "counts")
#' @param min_cells Minimum cells per sample (default: 10)
#' @param min_counts Minimum total counts per gene (default: 10)
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param logfc_threshold Log fold change threshold (default: 0)
#'
#' @return Data frame with differential expression results, or list of data frames if mode="per_cluster"
#'
#' @examples
#' \dontrun{
#' # Overall analysis
#' deg_results <- run_pseudobulk_deg(
#'   seurat_obj, 
#'   sample_col = "sample",
#'   group_col = "condition",
#'   comparison = c("Treatment", "Control")
#' )
#' 
#' # Per-cluster analysis
#' deg_list <- run_pseudobulk_deg(
#'   seurat_obj,
#'   sample_col = "sample",
#'   group_col = "condition",
#'   comparison = c("Treatment", "Control"),
#'   cluster_col = "seurat_clusters",
#'   mode = "per_cluster"
#' )
#' }
#'
#' @export
run_pseudobulk_deg <- function(object,
                               sample_col = NULL,
                               group_col = NULL,
                               comparison = NULL,
                               cluster_col = NULL,
                               target_cluster = NULL,
                               mode = c("overall", "per_cluster", "specific_cluster"),
                               assay = "RNA",
                               slot = "counts",
                               min_cells = 10,
                               min_counts = 10,
                               fdr_threshold = 0.05,
                               logfc_threshold = 0) {
  
  mode <- match.arg(mode)
  
  # Check if object is already prepared pseudobulk data
  is_prepared <- is.list(object) && all(c("counts", "metadata") %in% names(object))
  
  if (!is_prepared && (is.null(sample_col) || is.null(group_col))) {
    stop("sample_col and group_col are required when object is a Seurat object")
  }
  
  if (is.null(comparison) || length(comparison) != 2) {
    stop("comparison must be a character vector of length 2 (e.g., c('Treatment', 'Control'))")
  }
  
  # Perform analysis based on mode
  if (mode == "overall") {
    if (is_prepared) {
      pb_data <- object
    } else {
      pb_data <- prepare_pseudobulk_edgeR(
        object, sample_col, group_col,
        cluster_col = NULL,
        target_cluster = NULL,
        assay = assay,
        slot = slot,
        min_cells = min_cells,
        min_counts = min_counts
      )
    }
    
    results <- .run_edger_analysis(
      pb_data$counts, 
      pb_data$metadata, 
      group_col,
      comparison,
      fdr_threshold,
      logfc_threshold
    )
    
    return(results)
    
  } else if (mode == "per_cluster") {
    if (is.null(cluster_col)) {
      stop("cluster_col is required for per_cluster mode")
    }
    
    if (is_prepared) {
      stop("object must be a Seurat object for per_cluster mode")
    }
    
    clusters <- unique(object@meta.data[[cluster_col]])
    message("Running pseudobulk DEG for ", length(clusters), " clusters")
    
    results_list <- list()
    
    for (clust in clusters) {
      message("\n--- Cluster ", clust, " ---")
      
      tryCatch({
        pb_data <- prepare_pseudobulk_edgeR(
          object, sample_col, group_col,
          cluster_col = cluster_col,
          target_cluster = clust,
          assay = assay,
          slot = slot,
          min_cells = min_cells,
          min_counts = min_counts
        )
        
        results_list[[as.character(clust)]] <- .run_edger_analysis(
          pb_data$counts, 
          pb_data$metadata, 
          group_col,
          comparison,
          fdr_threshold,
          logfc_threshold
        )
      }, error = function(e) {
        message("Failed for cluster ", clust, ": ", e$message)
      })
    }
    
    return(results_list)
    
  } else if (mode == "specific_cluster") {
    if (is.null(cluster_col) || is.null(target_cluster)) {
      stop("cluster_col and target_cluster are required for specific_cluster mode")
    }
    
    if (is_prepared) {
      pb_data <- object
    } else {
      pb_data <- prepare_pseudobulk_edgeR(
        object, sample_col, group_col,
        cluster_col = cluster_col,
        target_cluster = target_cluster,
        assay = assay,
        slot = slot,
        min_cells = min_cells,
        min_counts = min_counts
      )
    }
    
    results <- .run_edger_analysis(
      pb_data$counts, 
      pb_data$metadata, 
      group_col,
      comparison,
      fdr_threshold,
      logfc_threshold
    )
    
    return(results)
  }
}

#' Internal edgeR Analysis Function
#'
#' @keywords internal
.run_edger_analysis <- function(counts, 
                                metadata, 
                                group_col,
                                comparison,
                                fdr_threshold,
                                logfc_threshold) {
  
  # Check if groups are present
  groups <- metadata[[group_col]]
  
  if (!all(comparison %in% groups)) {
    stop("Comparison groups not found in ", group_col, ": ", 
         paste(setdiff(comparison, groups), collapse = ", "))
  }
  
  # Create DGEList
  dge <- edgeR::DGEList(counts = counts, group = groups)
  
  # Normalization
  dge <- edgeR::calcNormFactors(dge)
  
  # Design matrix
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(factor(groups))
  
  # Estimate dispersion
  dge <- edgeR::estimateDisp(dge, design)
  
  # Fit model
  fit <- edgeR::glmQLFit(dge, design)
  
  # Make contrast
  contrast_formula <- paste0(comparison[1], "-", comparison[2])
  contrast <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
  
  # Test
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  
  # Get results
  results <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
  
  # Add gene column
  results$gene <- rownames(results)
  
  # Filter by thresholds
  results$significant <- (results$FDR < fdr_threshold) & 
                         (abs(results$logFC) > logfc_threshold)
  
  # Reorder columns
  col_order <- c("gene", "logFC", "logCPM", "F", "PValue", "FDR", "significant")
  results <- results[, col_order]
  
  n_sig <- sum(results$significant)
  n_up <- sum(results$significant & results$logFC > 0)
  n_down <- sum(results$significant & results$logFC < 0)
  
  message("DEG results: ", n_sig, " significant genes (", n_up, " up, ", n_down, " down)")
  
  return(results)
}

