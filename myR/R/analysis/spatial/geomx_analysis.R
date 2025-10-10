#' GeoMx Digital Spatial Profiling Analysis
#'
#' This module provides functions for analyzing NanoString GeoMx Digital Spatial
#' Profiling data, including Q3 normalization, quality control, differential
#' expression, and visualization.
#'
#' @name geomx_analysis
NULL

#' Prepare GeoMx Data
#'
#' Prepares raw GeoMx count data for analysis by organizing into a standard format.
#'
#' @param counts Count matrix (genes x ROIs)
#' @param metadata Data frame with ROI-level metadata
#' @param feature_info Optional data frame with gene/probe information
#'
#' @return List containing:
#'   \item{counts}{Count matrix}
#'   \item{metadata}{Sample metadata}
#'   \item{features}{Feature information}
#'
#' @export
prepare_geomx_data <- function(counts, metadata, feature_info = NULL) {
  
  # Validate inputs
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts must be a matrix or data frame")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  
  # Convert to matrix if needed
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  
  # Check dimensions match
  if (ncol(counts) != nrow(metadata)) {
    stop("Number of columns in counts (", ncol(counts), 
         ") must match number of rows in metadata (", nrow(metadata), ")")
  }
  
  # Align column names
  if (!is.null(colnames(counts)) && !is.null(rownames(metadata))) {
    if (!all(colnames(counts) == rownames(metadata))) {
      warning("Column names of counts and row names of metadata don't match. Attempting to align...")
      common_samples <- intersect(colnames(counts), rownames(metadata))
      counts <- counts[, common_samples, drop = FALSE]
      metadata <- metadata[common_samples, , drop = FALSE]
    }
  } else {
    colnames(counts) <- rownames(metadata)
  }
  
  # Prepare feature info if not provided
  if (is.null(feature_info)) {
    feature_info <- data.frame(
      gene = rownames(counts),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    counts = counts,
    metadata = metadata,
    features = feature_info
  ))
}

#' Q3 Normalization for GeoMx Data
#'
#' Performs Q3 (third quartile) normalization on GeoMx count data.
#'
#' @param geomx_data List from prepare_geomx_data or count matrix
#' @param scale_factor Scaling factor after normalization (default: 1000)
#'
#' @return Normalized data in same format as input
#'
#' @examples
#' \dontrun{
#' geomx_data <- prepare_geomx_data(counts, metadata)
#' geomx_norm <- q3_normalize(geomx_data)
#' }
#'
#' @export
q3_normalize <- function(geomx_data, scale_factor = 1000) {
  
  # Extract counts if list format
  if (is.list(geomx_data) && "counts" %in% names(geomx_data)) {
    counts <- geomx_data$counts
    is_list_input <- TRUE
  } else {
    counts <- geomx_data
    is_list_input <- FALSE
  }
  
  # Calculate Q3 for each sample
  q3_values <- apply(counts, 2, function(x) quantile(x, probs = 0.75))
  
  # Calculate scaling factors
  scaling_factors <- q3_values / median(q3_values)
  
  # Normalize
  norm_counts <- sweep(counts, 2, scaling_factors, FUN = "/")
  
  # Apply scale factor
  norm_counts <- norm_counts * scale_factor
  
  message("Q3 normalization complete. Median Q3 value: ", round(median(q3_values), 2))
  
  # Return in same format as input
  if (is_list_input) {
    geomx_data$counts <- norm_counts
    geomx_data$q3_factors <- scaling_factors
    return(geomx_data)
  } else {
    return(norm_counts)
  }
}

#' Perform Quality Control on GeoMx Data
#'
#' Performs basic QC checks and filtering on GeoMx data.
#'
#' @param geomx_data List from prepare_geomx_data
#' @param min_counts Minimum total counts per ROI (default: 1000)
#' @param min_genes Minimum number of detected genes per ROI (default: 100)
#' @param min_rois Minimum number of ROIs expressing a gene (default: 2)
#'
#' @return Filtered GeoMx data list with QC metrics added to metadata
#'
#' @export
perform_qc <- function(geomx_data, 
                       min_counts = 1000,
                       min_genes = 100,
                       min_rois = 2) {
  
  counts <- geomx_data$counts
  metadata <- geomx_data$metadata
  
  # Calculate QC metrics
  metadata$total_counts <- colSums(counts)
  metadata$n_genes_detected <- colSums(counts > 0)
  
  # Identify ROIs passing QC
  pass_qc <- (metadata$total_counts >= min_counts) & 
             (metadata$n_genes_detected >= min_genes)
  
  n_fail <- sum(!pass_qc)
  if (n_fail > 0) {
    message("Removing ", n_fail, " ROIs failing QC")
  }
  
  # Filter ROIs
  counts <- counts[, pass_qc, drop = FALSE]
  metadata <- metadata[pass_qc, , drop = FALSE]
  
  # Filter genes
  genes_detected <- rowSums(counts > 0) >= min_rois
  n_genes_removed <- sum(!genes_detected)
  
  if (n_genes_removed > 0) {
    message("Removing ", n_genes_removed, " genes detected in < ", min_rois, " ROIs")
  }
  
  counts <- counts[genes_detected, , drop = FALSE]
  
  if (!is.null(geomx_data$features)) {
    geomx_data$features <- geomx_data$features[genes_detected, , drop = FALSE]
  }
  
  geomx_data$counts <- counts
  geomx_data$metadata <- metadata
  
  message("After QC: ", nrow(counts), " genes x ", ncol(counts), " ROIs")
  
  return(geomx_data)
}

#' Find Differentially Expressed Genes in GeoMx Data
#'
#' Performs differential expression analysis on GeoMx data using limma or Wilcoxon test.
#'
#' @param geomx_data List from prepare_geomx_data (normalized)
#' @param group_col Metadata column for grouping
#' @param comparison Comparison to make (e.g., c("Treatment", "Control"))
#' @param method Method to use: "limma" or "wilcox" (default: "limma")
#' @param covariates Optional covariate columns to include in limma model
#' @param log_transform Apply log2(x+1) transformation (default: TRUE for limma)
#'
#' @return Data frame with differential expression results
#'
#' @export
find_deg_geomx <- function(geomx_data,
                           group_col,
                           comparison,
                           method = c("limma", "wilcox"),
                           covariates = NULL,
                           log_transform = NULL) {
  
  method <- match.arg(method)
  
  # Set default for log_transform
  if (is.null(log_transform)) {
    log_transform <- (method == "limma")
  }
  
  counts <- geomx_data$counts
  metadata <- geomx_data$metadata
  
  # Check group column
  if (!group_col %in% colnames(metadata)) {
    stop("group_col '", group_col, "' not found in metadata")
  }
  
  # Check comparison groups
  groups <- metadata[[group_col]]
  if (!all(comparison %in% groups)) {
    stop("Comparison groups not found: ", 
         paste(setdiff(comparison, groups), collapse = ", "))
  }
  
  # Log transform if requested
  if (log_transform) {
    expr_data <- log2(counts + 1)
  } else {
    expr_data <- counts
  }
  
  if (method == "limma") {
    results <- .geomx_limma(expr_data, metadata, group_col, comparison, covariates)
  } else {
    results <- .geomx_wilcox(expr_data, metadata, group_col, comparison)
  }
  
  return(results)
}

#' Internal limma Analysis
#' @keywords internal
.geomx_limma <- function(expr_data, metadata, group_col, comparison, covariates) {
  
  # Create design matrix
  groups <- factor(metadata[[group_col]])
  
  if (is.null(covariates)) {
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- levels(groups)
  } else {
    # Include covariates
    covar_formula <- paste(covariates, collapse = " + ")
    formula_str <- paste("~ 0 + groups +", covar_formula)
    design <- model.matrix(as.formula(formula_str), data = metadata)
    colnames(design) <- gsub("groups", "", colnames(design))
  }
  
  # Fit model
  fit <- limma::lmFit(expr_data, design)
  
  # Make contrast
  contrast_formula <- paste0(comparison[1], "-", comparison[2])
  contrast <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
  
  # Test
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)
  
  # Get results
  results <- limma::topTable(fit2, number = Inf, sort.by = "P")
  
  # Standardize column names
  results <- data.frame(
    gene = rownames(results),
    logFC = results$logFC,
    AveExpr = results$AveExpr,
    t = results$t,
    pvalue = results$P.Value,
    padj = results$adj.P.Val,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' Internal Wilcoxon Test
#' @keywords internal
.geomx_wilcox <- function(expr_data, metadata, group_col, comparison) {
  
  groups <- metadata[[group_col]]
  
  group1_samples <- groups == comparison[1]
  group2_samples <- groups == comparison[2]
  
  results <- data.frame(
    gene = rownames(expr_data),
    logFC = numeric(nrow(expr_data)),
    mean_group1 = numeric(nrow(expr_data)),
    mean_group2 = numeric(nrow(expr_data)),
    pvalue = numeric(nrow(expr_data)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(expr_data))) {
    group1_vals <- expr_data[i, group1_samples]
    group2_vals <- expr_data[i, group2_samples]
    
    results$mean_group1[i] <- mean(group1_vals)
    results$mean_group2[i] <- mean(group2_vals)
    results$logFC[i] <- results$mean_group1[i] - results$mean_group2[i]
    
    # Wilcoxon test
    test_result <- wilcox.test(group1_vals, group2_vals)
    results$pvalue[i] <- test_result$p.value
  }
  
  # Adjust p-values
  results$padj <- p.adjust(results$pvalue, method = "BH")
  
  # Sort by p-value
  results <- results[order(results$pvalue), ]
  
  return(results)
}

#' Plot Volcano Plot for GeoMx DEG Results
#'
#' Creates a volcano plot from differential expression results.
#'
#' @param deg_results Data frame from find_deg_geomx
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param logfc_threshold Log fold change threshold (default: 1)
#' @param label_top Number of top genes to label (default: 10)
#' @param title Plot title
#'
#' @return ggplot object
#'
#' @export
plot_deg_volcano <- function(deg_results,
                             padj_threshold = 0.05,
                             logfc_threshold = 1,
                             label_top = 10,
                             title = "Volcano Plot") {
  
  # Add significance category
  deg_results$category <- "NS"
  deg_results$category[deg_results$padj < padj_threshold & 
                       deg_results$logFC > logfc_threshold] <- "Up"
  deg_results$category[deg_results$padj < padj_threshold & 
                       deg_results$logFC < -logfc_threshold] <- "Down"
  
  # Get top genes to label
  top_genes <- head(deg_results[deg_results$category != "NS", ], label_top)
  
  # Create plot
  p <- ggplot2::ggplot(deg_results, ggplot2::aes(x = logFC, y = -log10(padj), color = category)) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    ggplot2::geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
                       linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(padj_threshold), 
                       linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Regulation"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  # Add labels if ggrepel is available
  if (nrow(top_genes) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = 20
    )
  }
  
  return(p)
}

#' Plot Heatmap for GeoMx DEG Results
#'
#' Creates a heatmap of top differentially expressed genes.
#'
#' @param geomx_data GeoMx data list (normalized)
#' @param deg_results Data frame from find_deg_geomx
#' @param group_col Metadata column for grouping
#' @param n_genes Number of top genes to plot (default: 50)
#' @param scale_rows Z-score scale rows (default: TRUE)
#'
#' @return Heatmap (ComplexHeatmap or pheatmap object)
#'
#' @export
plot_deg_heatmap <- function(geomx_data,
                             deg_results,
                             group_col,
                             n_genes = 50,
                             scale_rows = TRUE) {
  
  # Get top genes
  top_genes <- head(deg_results$gene, n_genes)
  
  # Extract expression
  expr_data <- geomx_data$counts[top_genes, , drop = FALSE]
  
  # Log transform if not already
  if (max(expr_data) > 100) {
    expr_data <- log2(expr_data + 1)
  }
  
  # Z-score if requested
  if (scale_rows) {
    expr_data <- t(scale(t(expr_data)))
  }
  
  # Get annotation
  annotation_col <- data.frame(
    Group = geomx_data$metadata[[group_col]],
    row.names = colnames(expr_data)
  )
  
  # Create heatmap
  pheatmap::pheatmap(
    expr_data,
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = (n_genes <= 50),
    show_colnames = TRUE,
    main = paste("Top", n_genes, "DEGs")
  )
}

#' Run Complete GeoMx Analysis Pipeline
#'
#' Wrapper function to run the complete GeoMx analysis workflow.
#'
#' @param counts Count matrix
#' @param metadata Sample metadata
#' @param group_col Column for differential expression
#' @param comparison Two-element vector for comparison
#' @param normalize Perform Q3 normalization (default: TRUE)
#' @param qc Perform QC filtering (default: TRUE)
#' @param method DE method: "limma" or "wilcox" (default: "limma")
#' @param plot Create plots (default: TRUE)
#'
#' @return List containing processed data, DEG results, and plots
#'
#' @export
run_geomx_analysis <- function(counts,
                               metadata,
                               group_col,
                               comparison,
                               normalize = TRUE,
                               qc = TRUE,
                               method = "limma",
                               plot = TRUE) {
  
  message("=== GeoMx Analysis Pipeline ===")
  
  # Prepare data
  message("\n1. Preparing data...")
  geomx_data <- prepare_geomx_data(counts, metadata)
  
  # Normalize
  if (normalize) {
    message("\n2. Q3 normalization...")
    geomx_data <- q3_normalize(geomx_data)
  }
  
  # QC
  if (qc) {
    message("\n3. Quality control...")
    geomx_data <- perform_qc(geomx_data)
  }
  
  # Differential expression
  message("\n4. Differential expression (", method, ")...")
  deg_results <- find_deg_geomx(
    geomx_data,
    group_col = group_col,
    comparison = comparison,
    method = method
  )
  
  n_sig <- sum(deg_results$padj < 0.05)
  message("Found ", n_sig, " significant DEGs (padj < 0.05)")
  
  # Create plots
  plots <- list()
  if (plot) {
    message("\n5. Creating plots...")
    
    plots$volcano <- plot_deg_volcano(
      deg_results,
      title = paste(comparison[1], "vs", comparison[2])
    )
    
    plots$heatmap <- plot_deg_heatmap(
      geomx_data,
      deg_results,
      group_col = group_col
    )
  }
  
  message("\n=== Analysis Complete ===")
  
  return(list(
    data = geomx_data,
    deg_results = deg_results,
    plots = plots
  ))
}

