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
                               normalize_transpose = FALSE,
                               remove_na = FALSE,
                               show_group_separator = TRUE,
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
  
  # If normalize_transpose is FALSE, use column normalization (sample/group level)
  # If TRUE, use row normalization (gene level) - default behavior
  if (!normalize_transpose) {
    # User wants column normalization (sample/group level)
    if (normalize_by == "row") {
      normalize_by <- "column"
    } else if (normalize_by == "column") {
      normalize_by <- "row"  # Swap back if already column
    }
    # "both" stays "both"
  }
  # If normalize_transpose is TRUE (default), keep row normalization (gene level)
  
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
  
  # Handle split.by for faceting
  has_split <- !is.null(split.by) && split.by %in% names(plot_df)
  
  if (has_split) {
    # Split data by split.by variable
    # Include NA as a separate level if present
    split_values_all <- unique(plot_df[[split.by]])
    has_na <- any(is.na(split_values_all))
    
    # Separate NA and non-NA values
    split_values <- split_values_all[!is.na(split_values_all)]
    
    # If only NA values exist, warn but continue
    if (length(split_values) == 0 && !has_na) {
      stop("No valid split.by values found")
    }
    
    # If we have both NA and non-NA, include NA as a level
    if (has_na && length(split_values) > 0) {
      # Convert NA to "NA" string for processing
      plot_df[[split.by]] <- as.character(plot_df[[split.by]])
      plot_df[[split.by]][is.na(plot_df[[split.by]])] <- "NA"
      split_values <- c(split_values, "NA")
    } else if (has_na && length(split_values) == 0) {
      # Only NA values - this is the problematic case (6번 플롯 문제)
      warning("Only NA values found in split.by. Filtering out NA values or using group.by only.")
      # Don't use split.by, treat as no split
      has_split <- FALSE
    }
    
    # Reshape for each split value
    heatmap_list <- list()
    for (split_val in split_values) {
      # Handle "NA" string vs actual NA
      if (split_val == "NA") {
        plot_df_split <- plot_df[is.na(plot_df[[split.by]]) | plot_df[[split.by]] == "NA", , drop = FALSE]
      } else {
        plot_df_split <- plot_df[plot_df[[split.by]] == split_val & !is.na(plot_df[[split.by]]), , drop = FALSE]
      }
      
      # Skip if no data for this split
      if (nrow(plot_df_split) == 0) {
        warning("No data found for split.by value: ", split_val)
        next
      }
      
      heatmap_matrix <- .reshape_for_heatmap(
        plot_df = plot_df_split,
        features = features,
        row_var = group.by,
        col_var = NULL,  # Features become columns
        aggregate = FALSE,  # Already aggregated
        agg_fun = agg_fun
      )
      
      heatmap_list[[as.character(split_val)]] <- heatmap_matrix
    }
    
    # Combine matrices with split labels
    # Create a combined data frame with split.by column
    heatmap_df_list <- list()
    for (split_val in names(heatmap_list)) {
      hm <- heatmap_list[[split_val]]
      hm_df <- as.data.frame(hm) %>%
        tibble::rownames_to_column(var = group.by) %>%
        tidyr::pivot_longer(
          cols = -dplyr::all_of(group.by),
          names_to = "Gene",
          values_to = "Zscore"
        )
      hm_df[[split.by]] <- split_val
      heatmap_df_list[[split_val]] <- hm_df
    }
    
    heatmap_df <- dplyr::bind_rows(heatmap_df_list)
    
    # Normalize across all splits (gene level normalization)
    if (normalize) {
      if (normalize_by %in% c("row", "both")) {
        # Normalize by gene across all groups and splits
        heatmap_df <- heatmap_df %>%
          dplyr::group_by(.data[["Gene"]]) %>%
          dplyr::mutate(
            Zscore = scale(.data[["Zscore"]])[, 1]
          ) %>%
          dplyr::ungroup()
      }
      if (normalize_by %in% c("column", "both")) {
        # Normalize by group within each split
        heatmap_df <- heatmap_df %>%
          dplyr::group_by(.data[[split.by]], .data[[group.by]]) %>%
          dplyr::mutate(
            Zscore = scale(.data[["Zscore"]])[, 1]
          ) %>%
          dplyr::ungroup()
      }
    }
    
    # Set factor levels
    all_row_names <- unique(heatmap_df[[group.by]])
    all_col_names <- unique(heatmap_df$Gene)
    heatmap_df[[group.by]] <- factor(heatmap_df[[group.by]], levels = all_row_names)
    heatmap_df$Gene <- factor(heatmap_df$Gene, levels = all_col_names)
    heatmap_df[[split.by]] <- factor(heatmap_df[[split.by]], levels = split_values)
    
    # Remove NA values if requested
    if (remove_na) {
      na_before <- nrow(heatmap_df)
      heatmap_df <- heatmap_df[!is.na(heatmap_df$Zscore), ]
      na_removed <- na_before - nrow(heatmap_df)
      if (na_removed > 0) {
        warning("Removed ", na_removed, " rows with NA Zscore values")
      }
    }
    
  } else {
    # No split.by - original logic
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
    
    # Convert to long format for ggplot
    heatmap_df <- as.data.frame(heatmap_matrix) %>%
      tibble::rownames_to_column(var = group.by) %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(group.by),
        names_to = "Gene",
        values_to = "Zscore"
      )
    
    # Set factor levels to preserve order (remove duplicates)
    unique_row_names <- unique(rownames(heatmap_matrix))
    unique_col_names <- unique(colnames(heatmap_matrix))
    heatmap_df[[group.by]] <- factor(heatmap_df[[group.by]], levels = unique_row_names)
    heatmap_df$Gene <- factor(heatmap_df$Gene, levels = unique_col_names)
    
    # Remove NA values if requested
    if (remove_na) {
      na_before <- nrow(heatmap_df)
      heatmap_df <- heatmap_df[!is.na(heatmap_df$Zscore), ]
      na_removed <- na_before - nrow(heatmap_df)
      if (na_removed > 0) {
        warning("Removed ", na_removed, " rows with NA Zscore values")
      }
    }
  }
  
  # Cluster if requested (only if no split)
  if (!has_split) {
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
  }
  
  # Return data if requested
  if (return_data) {
    if (has_split) {
      return(heatmap_df)
    } else {
      return(as.data.frame(heatmap_matrix))
    }
  }
  
  # Create plot
  if (has_split) {
    # Faceted plot with split.by
    p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(
      x = .data[[group.by]],
      y = .data[["Gene"]],
      fill = .data[["Zscore"]]
    )) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(ggplot2::vars(.data[[split.by]]), scales = "free_x", nrow = 1) +
      # Add vertical line between facets (g3=1 and g3=2)
      ggplot2::geom_vline(
        data = data.frame(x = 0.5, split_val = split_values),
        ggplot2::aes(xintercept = x),
        linetype = "solid",
        color = "black",
        linewidth = 1,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_gradient2(
        low = if(is.null(color_palette)) "blue" else color_palette[1],
        mid = if(is.null(color_palette)) "white" else color_palette[2],
        high = if(is.null(color_palette)) "red" else color_palette[3],
        midpoint = 0,
        name = if(normalize) "Z-score" else "Expression",
        na.value = "grey90"
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
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
        strip.text = ggplot2::element_text(face = "bold", size = 12),
        # Add border between panels to separate g3=1 and g3=2
        panel.spacing = ggplot2::unit(0.05, "lines"),
        strip.background = ggplot2::element_rect(fill = "grey90", color = "black", linewidth = 1.5)
      )
    
    # Log facet information for user reference
    if (nrow(samples_per_split) >= 2 && show_group_separator) {
      # Get number of samples in first facet
      n_first_facet <- samples_per_split$n_samples[1]
      split_levels <- levels(heatmap_df[[split.by]])
      if (length(split_levels) >= 2) {
        message("Heatmap facets: ", split_levels[1], " has ", n_first_facet, " samples, ", 
                split_levels[2], " has ", samples_per_split$n_samples[2], " samples")
      }
      # Note: Visual separation is provided by strip.background border and panel.spacing
      # In facet_wrap with free_x scales, direct vertical lines between facets are complex
    }
  } else {
    # Single plot without split
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
    
    # Add vertical separator lines between groups if requested
    if (show_group_separator && nlevels(heatmap_df[[group.by]]) > 1) {
      # Calculate positions for vertical lines between groups
      group_levels <- levels(heatmap_df[[group.by]])
      n_groups <- length(group_levels)
      
      if (n_groups > 1) {
        # Get x positions (0.5 between each group)
        # Convert factor positions to numeric
        x_positions <- seq(0.5, n_groups - 0.5, by = 1)
        x_positions <- x_positions[-c(1, length(x_positions))]  # Remove first and last
        
        if (length(x_positions) > 0) {
          p <- p + ggplot2::geom_vline(
            xintercept = x_positions,
            linetype = "dashed",
            color = "black",
            linewidth = 0.5,
            alpha = 0.5
          )
        }
      }
    }
  }
  
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

