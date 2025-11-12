#' Cell Aggregation Utilities
#'
#' Functions for aggregating cells by various grouping variables
#' (cluster, sample, patient, etc.) for pseudobulk analysis.
#'
#' @name utils_aggregation
NULL

#' Aggregate Cells by Grouping Variables
#'
#' Aggregates cell-level data by specified grouping variables.
#' Supports multiple levels of aggregation (e.g., patient -> group).
#'
#' @param plot_df data.frame with cell-level data
#' @param features Character vector of features to aggregate
#' @param group_vars Character vector of grouping variables, in order of aggregation.
#'   Example: c("patient", "cell_type") aggregates by patient first, then by cell_type.
#' @param fun Character. Aggregation function: "mean" (default), "sum", "median", "count"
#' @param min_cells Integer. Minimum number of cells per group (default: 1)
#'
#' @return data.frame with aggregated data
#'
#' @examples
#' \dontrun{
#' # Aggregate by patient
#' agg_df <- .aggregate_cells(plot_df, features = c("CD3D", "CD4"), group_vars = "patient")
#'
#' # Multi-level aggregation: patient -> cell_type
#' agg_df <- .aggregate_cells(plot_df, features = "CD3D", group_vars = c("patient", "cell_type"))
#' }
#'
#' @keywords internal
.aggregate_cells <- function(plot_df,
                            features,
                            group_vars,
                            fun = "mean",
                            min_cells = 1) {
  
  if (!is.data.frame(plot_df)) {
    stop("plot_df must be a data.frame")
  }
  
  # Check required columns
  missing_cols <- setdiff(c(features, group_vars), names(plot_df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Select aggregation function
  agg_fun <- switch(
    fun,
    "mean" = function(x) mean(x, na.rm = TRUE),
    "sum" = function(x) sum(x, na.rm = TRUE),
    "median" = function(x) median(x, na.rm = TRUE),
    "count" = function(x) length(x),
    stop("Unknown aggregation function: ", fun)
  )
  
  # Aggregate
  agg_df <- plot_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(features),
        agg_fun,
        .names = "{.col}"
      ),
      .n_cells = dplyr::n(),
      .groups = "drop"
    )
  
  # Filter by minimum cells
  if (min_cells > 1) {
    agg_df <- agg_df %>%
      dplyr::filter(.n_cells >= min_cells)
  }
  
  return(agg_df)
}

#' Prepare Data with Flexible Aggregation
#'
#' Main function for preparing plot data with optional aggregation.
#' Supports Seurat/df input, metadata joining, and multi-level aggregation.
#'
#' @param data Seurat object or data.frame
#' @param features Character vector of features
#' @param group.by Character. Grouping variable for faceting
#' @param split.by Character. Splitting variable for x-axis/color
#' @param sample_col Character. Sample/patient ID column (for aggregation)
#' @param aggregate Logical. Whether to aggregate cells (default: FALSE)
#' @param aggregate_by Character vector. Variables to aggregate by.
#'   If NULL and sample_col is provided, uses sample_col.
#'   If NULL and aggregate=TRUE, uses group.by and split.by.
#' @param agg_fun Character. Aggregation function (default: "mean")
#' @param min_cells Integer. Minimum cells per group for aggregation
#' @param assay Character. Assay to use (Seurat only)
#' @param layer Character. Layer to use (Seurat only)
#' @param metadata_df Optional data.frame. Additional metadata
#' @param match_id Character. Column for matching metadata
#'
#' @return data.frame ready for plotting
#'
#' @examples
#' \dontrun{
#' # Cell-level data
#' plot_df <- .prepare_data_with_aggregation(
#'   data = sobj,
#'   features = "CD3D",
#'   group.by = "cell_type",
#'   split.by = "condition"
#' )
#'
#' # Aggregated by patient
#' plot_df <- .prepare_data_with_aggregation(
#'   data = sobj,
#'   features = "CD3D",
#'   group.by = "cell_type",
#'   split.by = "condition",
#'   sample_col = "patient",
#'   aggregate = TRUE
#' )
#' }
#'
#' @keywords internal
.prepare_data_with_aggregation <- function(data,
                                           features,
                                           group.by = NULL,
                                           split.by = NULL,
                                           sample_col = NULL,
                                           aggregate = FALSE,
                                           aggregate_by = NULL,
                                           agg_fun = "mean",
                                           min_cells = 1,
                                           assay = NULL,
                                           layer = "data",
                                           metadata_df = NULL,
                                           match_id = NULL) {
  
  # Prepare base data
  plot_df <- .prepare_plot_data(
    data = data,
    features = features,
    group.by = group.by,
    split.by = split.by,
    sample_col = sample_col,
    assay = assay,
    layer = layer,
    metadata_df = metadata_df,
    match_id = match_id
  )
  
  # Determine aggregation variables
  if (aggregate) {
    if (is.null(aggregate_by)) {
      # Default: use sample_col if provided, otherwise group.by + split.by
      if (!is.null(sample_col)) {
        aggregate_by <- sample_col
      } else {
        aggregate_by <- unique(c(group.by, split.by))
        aggregate_by <- aggregate_by[!is.null(aggregate_by)]
      }
    }
    
    if (length(aggregate_by) == 0) {
      warning("No aggregation variables specified, returning cell-level data")
      return(plot_df)
    }
    
    # Aggregate
    plot_df <- .aggregate_cells(
      plot_df = plot_df,
      features = features,
      group_vars = aggregate_by,
      fun = agg_fun,
      min_cells = min_cells
    )
  }
  
  return(plot_df)
}

#' Reshape Data for Heatmap
#'
#' Reshapes data from long to wide format for heatmap plotting.
#' Handles multiple grouping variables and feature aggregation.
#'
#' @param plot_df data.frame with plot data
#' @param features Character vector of features (genes)
#' @param row_var Character. Variable for heatmap rows (e.g., "cluster", "patient")
#' @param col_var Character. Variable for heatmap columns (e.g., "gene", "cell_type")
#' @param value_var Character. Variable containing values (default: uses first feature)
#' @param aggregate Logical. Whether data needs aggregation (default: FALSE)
#' @param agg_fun Character. Aggregation function if aggregate=TRUE
#'
#' @return data.frame in wide format (rows x columns matrix)
#'
#' @keywords internal
.reshape_for_heatmap <- function(plot_df,
                                  features,
                                  row_var,
                                  col_var = NULL,
                                  value_var = NULL,
                                  aggregate = FALSE,
                                  agg_fun = "mean") {
  
  # If col_var is NULL, features become columns
  if (is.null(col_var)) {
    # Check if features are already in columns (already in wide format)
    features_in_cols <- all(features %in% names(plot_df))
    
    if (features_in_cols) {
      # Already in wide format, just select relevant columns
      heatmap_df <- plot_df %>%
        dplyr::select(dplyr::all_of(c(row_var, features)))
    } else {
      # Need to pivot to long first, then aggregate if needed, then pivot to wide
      plot_df_long <- plot_df %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(features),
          names_to = ".feature",
          values_to = ".value"
        )
      
      # Aggregate if needed
      if (aggregate) {
        group_vars <- unique(c(row_var, ".feature"))
        group_vars <- group_vars[!is.null(group_vars)]
        
        plot_df_long <- .aggregate_cells(
          plot_df = plot_df_long,
          features = ".value",
          group_vars = group_vars,
          fun = agg_fun
        )
      }
      
      # Reshape to wide
      heatmap_df <- plot_df_long %>%
        tidyr::pivot_wider(
          id_cols = dplyr::all_of(row_var),
          names_from = ".feature",
          values_from = ".value",
          values_fill = NA
        )
    }
  } else {
    # col_var is specified
    if (is.null(value_var)) {
      value_var <- features[1]
    }
    
    # Aggregate if needed
    if (aggregate) {
      group_vars <- unique(c(row_var, col_var))
      group_vars <- group_vars[!is.null(group_vars)]
      
      plot_df <- .aggregate_cells(
        plot_df = plot_df,
        features = value_var,
        group_vars = group_vars,
        fun = agg_fun
      )
    }
    
    # Reshape to wide
    heatmap_df <- plot_df %>%
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(row_var),
        names_from = dplyr::all_of(col_var),
        values_from = dplyr::all_of(value_var),
        values_fill = NA
      )
  }
  
  # Convert to matrix format
  row_names <- heatmap_df[[row_var]]
  heatmap_df[[row_var]] <- NULL
  heatmap_matrix <- as.matrix(heatmap_df)
  rownames(heatmap_matrix) <- row_names
  
  return(heatmap_matrix)
}

