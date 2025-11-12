#' Heatmap Plot Functions
#'
#' Standardized heatmap functions for gene/gene set expression visualization
#' with flexible grouping and normalization.
#'
#' @name plots_heatmap
NULL

#' Heatmap: Gene Expression
#'
#' Creates a heatmap of gene expression across groups (clusters, samples, etc.).
#' Supports z-score normalization and flexible grouping.
#'
#' @param data Seurat object or data.frame
#' @param features Character vector. Genes/features to plot
#' @param group.by Character. Grouping variable (e.g., "seurat_clusters", "cell_type")
#' @param split.by Character. Optional secondary grouping variable
#' @param aggregate Logical. Whether to aggregate cells (default: TRUE)
#' @param aggregate_by Character vector. Variables to aggregate by (default: group.by)
#' @param agg_fun Character. Aggregation function: "mean" (default), "sum", "median"
#' @param normalize Logical. Apply z-score normalization (default: TRUE)
#' @param normalize_by Character. Normalize by: "row" (genes), "column" (groups), "both" (default: "row")
#' @param assay Character. Assay to use (Seurat only)
#' @param layer Character. Layer to use: "counts", "data", "scale.data" (default: "data")
#' @param metadata_df Optional data.frame. Additional metadata to join
#' @param match_id Character. Column name in metadata_df for matching
#' @param title Character. Plot title
#' @param x_label Character. X-axis label (default: group.by)
#' @param y_label Character. Y-axis label (default: "Genes")
#' @param cluster_rows Logical. Cluster rows (default: FALSE)
#' @param cluster_cols Logical. Cluster columns (default: FALSE)
#' @param show_rownames Logical. Show row names (default: TRUE)
#' @param show_colnames Logical. Show column names (default: TRUE)
#' @param color_palette Character vector. Color palette (default: blue-white-red)
#' @param return_data Logical. Return data instead of plot (default: FALSE)
#'
#' @return ggplot object or data.frame (if return_data=TRUE)
#'
#' @examples
#' \dontrun{
#' # Basic gene heatmap
#' p <- plot_heatmap_genes(sobj, features = c("CD3D", "CD4", "CD8A"), 
#'                         group.by = "seurat_clusters")
#'
#' # With split.by
#' p <- plot_heatmap_genes(sobj, features = genes, group.by = "cell_type",
#'                          split.by = "condition")
#'
#' # Return data
#' heatmap_data <- plot_heatmap_genes(sobj, features = genes, 
#'                                     group.by = "cluster", return_data = TRUE)
#' }
#'
#' @export
plot_heatmap_genes <- function(data,
                               features,
                               group.by,
                               split.by = NULL,
                               aggregate = TRUE,
                               aggregate_by = NULL,
                               agg_fun = "mean",
                               normalize = TRUE,
                               normalize_by = c("row", "column", "both"),
                               assay = NULL,
                               layer = "data",
                               metadata_df = NULL,
                               match_id = NULL,
                               title = NULL,
                               x_label = NULL,
                               y_label = NULL,
                               cluster_rows = FALSE,
                               cluster_cols = FALSE,
                               show_rownames = TRUE,
                               show_colnames = TRUE,
                               color_palette = NULL,
                               return_data = FALSE) {
  
  normalize_by <- match.arg(normalize_by)
  
  # Helper function for NULL coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  if (is.null(features) || length(features) == 0) {
    stop("features must be specified")
  }
  
  # Determine aggregation variables
  if (is.null(aggregate_by)) {
    aggregate_by <- unique(c(group.by, split.by))
    aggregate_by <- aggregate_by[!is.null(aggregate_by)]
  }
  
  # Prepare data
  plot_df <- .prepare_data_with_aggregation(
    data = data,
    features = features,
    group.by = group.by,
    split.by = split.by,
    aggregate = aggregate,
    aggregate_by = aggregate_by,
    agg_fun = agg_fun,
    assay = assay,
    layer = layer,
    metadata_df = metadata_df,
    match_id = match_id
  )
  
  # Reshape for heatmap
  heatmap_matrix <- .reshape_for_heatmap(
    plot_df = plot_df,
    features = features,
    row_var = group.by,
    col_var = NULL,  # Features become columns
    aggregate = FALSE,  # Already aggregated
    agg_fun = agg_fun
  )
  
  # Normalize if requested
  if (normalize) {
    if (normalize_by %in% c("row", "both")) {
      # Normalize across groups (columns) for each gene
      heatmap_matrix <- t(scale(t(heatmap_matrix)))
    }
    if (normalize_by %in% c("column", "both")) {
      # Normalize across genes (rows) for each group
      heatmap_matrix <- scale(heatmap_matrix)
    }
  }
  
  # Handle single group case
  if (ncol(heatmap_matrix) == 1) {
    if (normalize_by == "row") {
      # Normalize across genes
      heatmap_matrix <- scale(heatmap_matrix)
    } else {
      # Set to 0 for neutral
      heatmap_matrix[] <- 0
    }
  }
  
  # Cluster if requested
  if (cluster_rows && nrow(heatmap_matrix) > 1) {
    row_order <- hclust(dist(heatmap_matrix))$order
    heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]
  }
  
  if (cluster_cols && ncol(heatmap_matrix) > 1) {
    col_order <- hclust(dist(t(heatmap_matrix)))$order
    heatmap_matrix <- heatmap_matrix[, col_order, drop = FALSE]
  }
  
  # Sort row names if numeric
  row_names <- rownames(heatmap_matrix)
  numeric_test <- suppressWarnings(as.numeric(row_names))
  if (!all(is.na(numeric_test)) && sum(is.na(numeric_test)) == 0) {
    row_order <- order(numeric_test)
    heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]
  }
  
  # Return data if requested
  if (return_data) {
    return(as.data.frame(heatmap_matrix))
  }
  
  # Convert to long format for ggplot
  heatmap_df <- as.data.frame(heatmap_matrix) %>%
    tibble::rownames_to_column(var = group.by) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group.by),
      names_to = "Gene",
      values_to = "Zscore"
    )
  
  # Set factor levels to preserve order
  heatmap_df[[group.by]] <- factor(heatmap_df[[group.by]], levels = rownames(heatmap_matrix))
  heatmap_df$Gene <- factor(heatmap_df$Gene, levels = colnames(heatmap_matrix))
  
  # Create plot
  p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(
    x = .data[[group.by]],
    y = .data[["Gene"]],
    fill = .data[["Zscore"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = if(is.null(color_palette)) "blue" else color_palette[1],
      mid = if(is.null(color_palette)) "white" else color_palette[2],
      high = if(is.null(color_palette)) "red" else color_palette[3],
      midpoint = 0,
      name = if(normalize) "Z-score" else "Expression"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = if (is.null(title)) "Normalized Gene Expression" else title,
      x = if (is.null(x_label)) group.by else x_label,
      y = if (is.null(y_label)) "Genes" else y_label
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(face = "bold", size = 14),
      axis.title.y = ggplot2::element_text(face = "bold", size = 14),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  if (!show_rownames) {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }
  if (!show_colnames) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
  
  return(p)
}

#' Heatmap: Gene Sets (Metadata-based)
#'
#' Creates a heatmap of gene set expression. Gene sets should be pre-computed
#' as metadata columns (e.g., using AddModuleScore or similar).
#'
#' @param data Seurat object or data.frame
#' @param features Character vector. Gene set metadata column names
#' @param group.by Character. Grouping variable
#' @param split.by Character. Optional secondary grouping variable
#' @param aggregate Logical. Whether to aggregate cells (default: TRUE)
#' @param aggregate_by Character vector. Variables to aggregate by
#' @param agg_fun Character. Aggregation function (default: "mean")
#' @param normalize Logical. Apply z-score normalization (default: TRUE)
#' @param assay Character. Assay to use (Seurat only, not used for metadata)
#' @param layer Character. Layer to use (Seurat only, not used for metadata)
#' @param metadata_df Optional data.frame. Additional metadata
#' @param match_id Character. Column name in metadata_df for matching
#' @param title Character. Plot title
#' @param x_label Character. X-axis label
#' @param y_label Character. Y-axis label
#' @param assign_celltype Logical. Assign cell type based on highest gene set (default: FALSE)
#' @param return_data Logical. Return data instead of plot (default: FALSE)
#'
#' @return ggplot object or data.frame (if return_data=TRUE)
#'
#' @examples
#' \dontrun{
#' # Gene sets as metadata columns
#' sobj$T_cell_score <- AddModuleScore(...)
#' sobj$B_cell_score <- AddModuleScore(...)
#' 
#' p <- plot_heatmap_genesets(sobj, 
#'                             features = c("T_cell_score", "B_cell_score"),
#'                             group.by = "seurat_clusters")
#' }
#'
#' @export
plot_heatmap_genesets <- function(data,
                                   features,
                                   group.by,
                                   split.by = NULL,
                                   aggregate = TRUE,
                                   aggregate_by = NULL,
                                   agg_fun = "mean",
                                   normalize = TRUE,
                                   assay = NULL,
                                   layer = "data",
                                   metadata_df = NULL,
                                   match_id = NULL,
                                   title = NULL,
                                   x_label = NULL,
                                   y_label = NULL,
                                   assign_celltype = FALSE,
                                   return_data = FALSE) {
  
  `%||%` <- rlang::`%||%`
  
  if (is.null(features) || length(features) == 0) {
    stop("features (gene set column names) must be specified")
  }
  
  # Determine aggregation variables
  if (is.null(aggregate_by)) {
    aggregate_by <- unique(c(group.by, split.by))
    aggregate_by <- aggregate_by[!is.null(aggregate_by)]
  }
  
  # Prepare data (features are metadata columns)
  plot_df <- .prepare_data_with_aggregation(
    data = data,
    features = features,
    group.by = group.by,
    split.by = split.by,
    aggregate = aggregate,
    aggregate_by = aggregate_by,
    agg_fun = agg_fun,
    assay = assay,
    layer = layer,
    metadata_df = metadata_df,
    match_id = match_id
  )
  
  # Reshape for heatmap
  heatmap_matrix <- .reshape_for_heatmap(
    plot_df = plot_df,
    features = features,
    row_var = group.by,
    col_var = NULL,
    aggregate = FALSE,
    agg_fun = agg_fun
  )
  
  # Normalize
  if (normalize && nrow(heatmap_matrix) > 1) {
    heatmap_matrix <- t(scale(t(heatmap_matrix)))
  } else if (normalize && nrow(heatmap_matrix) == 1) {
    heatmap_matrix[] <- 0
  }
  
  # Assign cell type if requested
  if (assign_celltype && ncol(heatmap_matrix) > 0) {
    assigned <- apply(heatmap_matrix, 1, function(x) {
      if (all(is.na(x))) return(NA)
      colnames(heatmap_matrix)[which.max(x)]
    })
    # Add as attribute or column (for return_data)
    attr(heatmap_matrix, "assigned_celltype") <- assigned
  }
  
  # Sort row names if numeric
  row_names <- rownames(heatmap_matrix)
  numeric_test <- suppressWarnings(as.numeric(row_names))
  if (!all(is.na(numeric_test)) && sum(is.na(numeric_test)) == 0) {
    row_order <- order(numeric_test)
    heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]
  }
  
  # Return data if requested
  if (return_data) {
    result_df <- as.data.frame(heatmap_matrix)
    if (assign_celltype) {
      result_df$Assigned_CellType <- attr(heatmap_matrix, "assigned_celltype")
    }
    return(result_df)
  }
  
  # Convert to long format
  heatmap_df <- as.data.frame(heatmap_matrix) %>%
    tibble::rownames_to_column(var = group.by) %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(group.by),
      names_to = "GeneSet",
      values_to = "Zscore"
    )
  
  if (assign_celltype) {
    assigned_df <- data.frame(
      .row = rownames(heatmap_matrix),
      Assigned_CellType = attr(heatmap_matrix, "assigned_celltype")
    )
    names(assigned_df)[1] <- group.by
    heatmap_df <- dplyr::left_join(heatmap_df, assigned_df, by = group.by)
  }
  
  # Set factor levels
  heatmap_df[[group.by]] <- factor(heatmap_df[[group.by]], levels = rownames(heatmap_matrix))
  heatmap_df$GeneSet <- factor(heatmap_df$GeneSet, levels = colnames(heatmap_matrix))
  
  # Create plot
  p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(
    x = .data[[group.by]],
    y = .data[["GeneSet"]],
    fill = .data[["Zscore"]]
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = if(normalize) "Z-score" else "Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = if (is.null(title)) "Normalized Gene Set Expression" else title,
      x = if (is.null(x_label)) group.by else x_label,
      y = if (is.null(y_label)) "Gene Set" else y_label
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(face = "bold", size = 14),
      axis.title.y = ggplot2::element_text(face = "bold", size = 14),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

