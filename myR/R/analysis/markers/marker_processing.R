#' Marker Gene Processing Functions
#'
#' This module provides functions for processing, filtering, and formatting
#' marker gene results from Seurat's FindMarkers/FindAllMarkers.
#'
#' @name marker_processing
NULL

#' Trim Marker Results by Direction and Significance
#'
#' Filters marker gene results based on log fold change direction and 
#' adjusted p-value threshold.
#'
#' @param markers Data frame from FindMarkers/FindAllMarkers
#' @param direction Filter direction: "up" (positive logFC), "down" (negative logFC), 
#'   or "both" (default: "both")
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param logfc_col Name of log fold change column (default: "avg_log2FC")
#' @param padj_col Name of adjusted p-value column (default: "p_val_adj")
#'
#' @return Filtered data frame
#'
#' @examples
#' \dontrun{
#' # Keep only upregulated markers
#' up_markers <- marker_trim(all_markers, direction = "up")
#' 
#' # Keep both directions with stricter p-value
#' sig_markers <- marker_trim(all_markers, padj_threshold = 0.01)
#' }
#'
#' @export
marker_trim <- function(markers, 
                       direction = c("both", "up", "down"),
                       padj_threshold = 0.05,
                       logfc_col = "avg_log2FC",
                       padj_col = "p_val_adj") {
  
  direction <- match.arg(direction)
  
  # Check required columns
  if (!logfc_col %in% colnames(markers)) {
    stop("Column '", logfc_col, "' not found in markers")
  }
  if (!padj_col %in% colnames(markers)) {
    stop("Column '", padj_col, "' not found in markers")
  }
  
  # Filter by significance
  markers <- markers[markers[[padj_col]] < padj_threshold, , drop = FALSE]
  
  # Filter by direction
  if (direction == "up") {
    markers <- markers[markers[[logfc_col]] > 0, , drop = FALSE]
  } else if (direction == "down") {
    markers <- markers[markers[[logfc_col]] < 0, , drop = FALSE]
  }
  
  return(markers)
}

#' Filter Out Unwanted Gene Types
#'
#' Removes ribosomal, mitochondrial, hemoglobin, and other unwanted genes
#' from marker results.
#'
#' @param markers Data frame from FindMarkers/FindAllMarkers
#' @param remove_ribo Remove ribosomal genes (^RP[SL]) (default: TRUE)
#' @param remove_mito Remove mitochondrial genes (^MT-) (default: TRUE)
#' @param remove_hb Remove hemoglobin genes (^HB[AB]) (default: TRUE)
#' @param remove_ig Remove immunoglobulin genes (^IG[HKL]) (default: FALSE)
#' @param remove_tr Remove T-cell receptor genes (^TR[ABGD]) (default: FALSE)
#' @param custom_patterns Additional regex patterns to remove (default: NULL)
#' @param gene_col Name of gene column (default: "gene" or rownames if not present)
#'
#' @return Filtered data frame
#'
#' @examples
#' \dontrun{
#' # Remove ribosomal and mitochondrial genes only
#' clean_markers <- marker_filter(all_markers)
#' 
#' # Also remove immunoglobulin genes
#' clean_markers <- marker_filter(all_markers, remove_ig = TRUE)
#' 
#' # Custom pattern removal
#' clean_markers <- marker_filter(all_markers, custom_patterns = c("^LINC", "^MIR"))
#' }
#'
#' @export
marker_filter <- function(markers, 
                         remove_ribo = TRUE,
                         remove_mito = TRUE,
                         remove_hb = TRUE,
                         remove_ig = FALSE,
                         remove_tr = FALSE,
                         custom_patterns = NULL,
                         gene_col = "gene") {
  
  # Get gene names
  if (gene_col %in% colnames(markers)) {
    genes <- markers[[gene_col]]
  } else if (!is.null(rownames(markers))) {
    genes <- rownames(markers)
  } else {
    stop("Cannot identify gene column. Specify gene_col or ensure rownames are set.")
  }
  
  # Build filter patterns
  patterns <- character(0)
  
  if (remove_ribo) patterns <- c(patterns, "^RP[SL]")
  if (remove_mito) patterns <- c(patterns, "^MT-")
  if (remove_hb) patterns <- c(patterns, "^HB[AB]")
  if (remove_ig) patterns <- c(patterns, "^IG[HKL]")
  if (remove_tr) patterns <- c(patterns, "^TR[ABGD]")
  
  if (!is.null(custom_patterns)) {
    patterns <- c(patterns, custom_patterns)
  }
  
  if (length(patterns) == 0) {
    return(markers)
  }
  
  # Combine patterns
  combined_pattern <- paste(patterns, collapse = "|")
  
  # Filter
  keep <- !grepl(combined_pattern, genes, ignore.case = FALSE)
  
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message("Removed ", n_removed, " genes matching unwanted patterns")
  }
  
  return(markers[keep, , drop = FALSE])
}

#' Remove Intercept Terms from Linear Model Results
#'
#' Filters out "(Intercept)" terms from linear mixed model results.
#'
#' @param results Data frame with model results
#' @param term_col Name of column containing term names (default: "term")
#'
#' @return Filtered data frame
#'
#' @export
lrf <- function(results, term_col = "term") {
  if (!term_col %in% colnames(results)) {
    warning("Column '", term_col, "' not found. Returning unfiltered results.")
    return(results)
  }
  
  results[results[[term_col]] != "(Intercept)", , drop = FALSE]
}

#' Convert All Markers to Named List
#'
#' Converts FindAllMarkers output to a named list, with one element per cluster.
#'
#' @param all_markers Data frame from FindAllMarkers
#' @param cluster_col Name of cluster column (default: "cluster")
#' @param gene_col Name of gene column (default: "gene" or rownames)
#'
#' @return Named list of marker data frames, one per cluster
#'
#' @examples
#' \dontrun{
#' all_markers <- FindAllMarkers(seurat_obj)
#' marker_list <- all_markers_to_list(all_markers)
#' # Access cluster 0 markers
#' cluster0_markers <- marker_list[["0"]]
#' }
#'
#' @export
all_markers_to_list <- function(all_markers, 
                               cluster_col = "cluster",
                               gene_col = "gene") {
  
  if (!cluster_col %in% colnames(all_markers)) {
    stop("Column '", cluster_col, "' not found in all_markers")
  }
  
  # Ensure gene column exists
  if (!gene_col %in% colnames(all_markers) && !is.null(rownames(all_markers))) {
    all_markers[[gene_col]] <- rownames(all_markers)
  }
  
  # Split by cluster
  marker_list <- split(all_markers, all_markers[[cluster_col]])
  
  return(marker_list)
}

#' Print Top Marker Genes for All Clusters
#'
#' Prints top N marker genes for each cluster from FindAllMarkers output.
#'
#' @param all_markers Data frame from FindAllMarkers
#' @param n Number of top markers to print per cluster (default: 10)
#' @param order_by Column to order by (default: "avg_log2FC")
#' @param decreasing Sort in decreasing order (default: TRUE)
#' @param cluster_col Name of cluster column (default: "cluster")
#' @param gene_col Name of gene column (default: "gene" or rownames)
#'
#' @return Invisibly returns the input data frame
#'
#' @export
marker_print_all <- function(all_markers, 
                            n = 10,
                            order_by = "avg_log2FC",
                            decreasing = TRUE,
                            cluster_col = "cluster",
                            gene_col = "gene") {
  
  if (!cluster_col %in% colnames(all_markers)) {
    stop("Column '", cluster_col, "' not found")
  }
  
  # Ensure gene column
  if (!gene_col %in% colnames(all_markers) && !is.null(rownames(all_markers))) {
    all_markers[[gene_col]] <- rownames(all_markers)
  }
  
  clusters <- unique(all_markers[[cluster_col]])
  
  for (clust in clusters) {
    cat("\n=== Cluster", clust, "===\n")
    cluster_markers <- all_markers[all_markers[[cluster_col]] == clust, ]
    marker_print(cluster_markers, n = n, order_by = order_by, 
                decreasing = decreasing, gene_col = gene_col)
  }
  
  invisible(all_markers)
}

#' Print Top Marker Genes
#'
#' Prints top N marker genes from a single marker data frame.
#'
#' @param markers Data frame from FindMarkers
#' @param n Number of top markers to print (default: 10)
#' @param order_by Column to order by (default: "avg_log2FC")
#' @param decreasing Sort in decreasing order (default: TRUE)
#' @param gene_col Name of gene column (default: "gene" or rownames)
#' @param columns Columns to display (default: c("avg_log2FC", "p_val_adj"))
#'
#' @return Invisibly returns the input data frame
#'
#' @export
marker_print <- function(markers, 
                        n = 10,
                        order_by = "avg_log2FC",
                        decreasing = TRUE,
                        gene_col = "gene",
                        columns = c("avg_log2FC", "p_val_adj")) {
  
  if (nrow(markers) == 0) {
    cat("No markers found.\n")
    return(invisible(markers))
  }
  
  # Ensure gene column
  if (!gene_col %in% colnames(markers) && !is.null(rownames(markers))) {
    markers[[gene_col]] <- rownames(markers)
  }
  
  # Order
  if (order_by %in% colnames(markers)) {
    markers <- markers[order(markers[[order_by]], decreasing = decreasing), ]
  }
  
  # Select top N
  top_n <- head(markers, n)
  
  # Select columns to display
  display_cols <- c(gene_col, intersect(columns, colnames(top_n)))
  
  print(top_n[, display_cols, drop = FALSE], row.names = FALSE)
  
  invisible(markers)
}

#' Synthesize Multiple Marker Results (Meta-Analysis)
#'
#' Combines p-values and log fold changes from multiple FindMarkers results
#' using Fisher's method for p-values and weighted mean for log fold changes.
#'
#' @param marker_list Named list of marker data frames
#' @param gene_col Name of gene column (default: "gene" or rownames)
#' @param logfc_col Name of log fold change column (default: "avg_log2FC")
#' @param pval_col Name of p-value column (default: "p_val")
#' @param weights Optional weights for each dataset (default: equal weights)
#'
#' @return Data frame with synthesized results
#'
#' @export
synthesize_markers <- function(marker_list, 
                              gene_col = "gene",
                              logfc_col = "avg_log2FC",
                              pval_col = "p_val",
                              weights = NULL) {
  
  if (length(marker_list) == 0) {
    stop("marker_list is empty")
  }
  
  # Equal weights if not provided
  if (is.null(weights)) {
    weights <- rep(1, length(marker_list))
  }
  
  if (length(weights) != length(marker_list)) {
    stop("Length of weights must match length of marker_list")
  }
  
  # Normalize weights
  weights <- weights / sum(weights)
  
  # Get all unique genes
  all_genes <- unique(unlist(lapply(marker_list, function(m) {
    if (gene_col %in% colnames(m)) {
      m[[gene_col]]
    } else {
      rownames(m)
    }
  })))
  
  # Initialize results
  results <- data.frame(
    gene = all_genes,
    combined_logfc = NA_real_,
    combined_pval = NA_real_,
    n_datasets = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # Combine for each gene
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]
    
    logfcs <- numeric(0)
    pvals <- numeric(0)
    gene_weights <- numeric(0)
    
    for (j in seq_along(marker_list)) {
      markers <- marker_list[[j]]
      
      # Get gene row
      if (gene_col %in% colnames(markers)) {
        gene_row <- markers[markers[[gene_col]] == gene, ]
      } else {
        gene_row <- markers[rownames(markers) == gene, ]
      }
      
      if (nrow(gene_row) > 0) {
        logfcs <- c(logfcs, gene_row[[logfc_col]][1])
        pvals <- c(pvals, gene_row[[pval_col]][1])
        gene_weights <- c(gene_weights, weights[j])
      }
    }
    
    if (length(logfcs) > 0) {
      # Weighted mean of log fold changes
      results$combined_logfc[i] <- weighted.mean(logfcs, gene_weights)
      
      # Fisher's method for combining p-values
      chi_sq <- -2 * sum(log(pvals))
      results$combined_pval[i] <- pchisq(chi_sq, df = 2 * length(pvals), lower.tail = FALSE)
      
      results$n_datasets[i] <- length(logfcs)
    }
  }
  
  # Adjust p-values
  results$combined_padj <- p.adjust(results$combined_pval, method = "BH")
  
  # Order by combined p-value
  results <- results[order(results$combined_pval), ]
  
  # Remove genes with no data
  results <- results[!is.na(results$combined_pval), ]
  
  return(results)
}

#' Synthesize Ranks from Multiple Marker Results
#'
#' Combines ranks from multiple FindMarkers results using geometric mean of ranks.
#'
#' @param marker_list Named list of marker data frames
#' @param gene_col Name of gene column (default: "gene" or rownames)
#' @param rank_by Column to rank by (default: "avg_log2FC")
#' @param decreasing Rank in decreasing order (default: TRUE)
#'
#' @return Data frame with synthesized ranks
#'
#' @export
synthesize_ranks <- function(marker_list, 
                            gene_col = "gene",
                            rank_by = "avg_log2FC",
                            decreasing = TRUE) {
  
  if (length(marker_list) == 0) {
    stop("marker_list is empty")
  }
  
  # Get all unique genes
  all_genes <- unique(unlist(lapply(marker_list, function(m) {
    if (gene_col %in% colnames(m)) {
      m[[gene_col]]
    } else {
      rownames(m)
    }
  })))
  
  # Initialize results
  results <- data.frame(
    gene = all_genes,
    geometric_mean_rank = NA_real_,
    n_datasets = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # Get ranks for each gene
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]
    ranks <- numeric(0)
    
    for (j in seq_along(marker_list)) {
      markers <- marker_list[[j]]
      
      # Rank markers
      if (rank_by %in% colnames(markers)) {
        markers$rank <- rank(-markers[[rank_by]] * ifelse(decreasing, 1, -1), 
                            ties.method = "average")
      } else {
        next
      }
      
      # Get gene rank
      if (gene_col %in% colnames(markers)) {
        gene_rank <- markers$rank[markers[[gene_col]] == gene]
      } else {
        gene_rank <- markers$rank[rownames(markers) == gene]
      }
      
      if (length(gene_rank) > 0) {
        ranks <- c(ranks, gene_rank[1])
      }
    }
    
    if (length(ranks) > 0) {
      # Geometric mean of ranks
      results$geometric_mean_rank[i] <- exp(mean(log(ranks)))
      results$n_datasets[i] <- length(ranks)
    }
  }
  
  # Order by geometric mean rank
  results <- results[order(results$geometric_mean_rank), ]
  
  # Remove genes with no data
  results <- results[!is.na(results$geometric_mean_rank), ]
  
  return(results)
}

