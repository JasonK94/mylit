#' Boxplot Functions
#'
#' Standardized boxplot functions for cell-level and pseudobulk visualization
#' with flexible grouping and statistical testing.
#'
#' @name plots_box
NULL

#' Boxplot: Cell-level or Aggregated
#'
#' Creates boxplots of features across groups. Supports both cell-level and
#' aggregated (pseudobulk) visualization.
#'
#' @param data Seurat object or data.frame
#' @param features Character vector. Features to plot (genes or metadata)
#' @param group.by Character. Grouping variable for faceting
#' @param split.by Character. Splitting variable for x-axis
#' @param idents Character vector. Specific identities to include from group.by
#' @param aggregate Logical. Whether to aggregate cells (default: FALSE)
#' @param sample_col Character. Sample/patient ID column (required if aggregate=TRUE)
#' @param aggregate_by Character vector. Variables to aggregate by
#' @param agg_fun Character. Aggregation function: "mean" (default), "sum", "median"
#' @param assay Character. Assay to use (Seurat only)
#' @param layer Character. Layer to use: "counts", "data", "scale.data" (default: "data")
#' @param metadata_df Optional data.frame. Additional metadata to join
#' @param match_id Character. Column name in metadata_df for matching
#' @param ncol Integer. Number of columns for multi-feature plots (default: 3)
#' @param pt.size Numeric. Point size for jitter (0 to hide, default: 0.2 for cell, 1.0 for aggregated)
#' @param violin Logical. Overlay violin plot (default: FALSE)
#' @param add_stats Logical. Add statistical comparisons (default: TRUE)
#' @param ... Additional arguments passed to stat_compare_means
#'
#' @return ggplot object (patchwork for multiple features)
#'
#' @examples
#' \dontrun{
#' # Cell-level boxplot
#' p <- plot_box(sobj, features = c("CD3D", "CD4"), 
#'               group.by = "cell_type", split.by = "condition")
#'
#' # Aggregated boxplot
#' p <- plot_box(sobj, features = "CD3D", group.by = "cell_type",
#'                split.by = "condition", sample_col = "patient", aggregate = TRUE)
#' }
#'
#' @export
plot_box <- function(data,
                     features,
                     group.by,
                     split.by,
                     idents = NULL,
                     aggregate = FALSE,
                     sample_col = NULL,
                     aggregate_by = NULL,
                     agg_fun = "mean",
                     assay = NULL,
                     layer = "data",
                     metadata_df = NULL,
                     match_id = NULL,
                     ncol = 3,
                     pt.size = NULL,
                     violin = FALSE,
                     add_stats = TRUE,
                     ...) {
  
  if (is.null(features) || length(features) == 0) {
    stop("features must be specified")
  }
  
  if (is.null(split.by)) {
    stop("split.by must be specified")
  }
  
  `%||%` <- rlang::`%||%`
  
  # Set default point size
  if (is.null(pt.size)) {
    pt.size <- if (aggregate) 1.0 else 0.2
  }
  
  # Determine aggregation variables
  if (aggregate) {
    if (is.null(sample_col) && is.null(aggregate_by)) {
      stop("sample_col or aggregate_by must be specified when aggregate=TRUE")
    }
    if (is.null(aggregate_by)) {
      aggregate_by <- unique(c(sample_col, group.by, split.by))
      aggregate_by <- aggregate_by[!is.null(aggregate_by)]
    } else {
      # Ensure group.by and split.by are included in aggregate_by
      # They are needed for plotting even if user didn't specify them
      if (!is.null(group.by) && !group.by %in% aggregate_by) {
        aggregate_by <- c(aggregate_by, group.by)
      }
      if (!is.null(split.by) && !split.by %in% aggregate_by) {
        aggregate_by <- c(aggregate_by, split.by)
      }
    }
  }
  
  # Prepare data
  plot_df <- .prepare_data_with_aggregation(
    data = data,
    features = features,
    group.by = group.by,
    split.by = split.by,
    sample_col = sample_col,
    aggregate = aggregate,
    aggregate_by = aggregate_by,
    agg_fun = agg_fun,
    assay = assay,
    layer = layer,
    metadata_df = metadata_df,
    match_id = match_id
  )
  
  # Handle group.by = NULL case
  group.by.internal <- group.by
  if (is.null(group.by)) {
    group.by.internal <- ".internal_placeholder_group"
    plot_df[[group.by.internal]] <- "Overall"
  } else {
    # Check if group.by column exists after aggregation
    if (!group.by.internal %in% names(plot_df)) {
      stop("Column '", group.by.internal, "' not found after aggregation. ",
           "Please include '", group.by.internal, "' in aggregate_by parameter.")
    }
  }
  
  # Check if split.by column exists
  if (!split.by %in% names(plot_df)) {
    stop("Column '", split.by, "' not found after aggregation. ",
         "Please include '", split.by, "' in aggregate_by parameter.")
  }
  
  # Filter by idents if specified
  if (!is.null(group.by) && !is.null(idents)) {
    original_col_type <- class(plot_df[[group.by.internal]])
    plot_df[[group.by.internal]] <- as.character(plot_df[[group.by.internal]])
    plot_df <- plot_df[plot_df[[group.by.internal]] %in% idents, ]
    
    if (nrow(plot_df) == 0) {
      warning("No data remaining after idents filtering")
      return(invisible(NULL))
    }
    
    if ("factor" %in% original_col_type) {
      # Preserve original factor levels
      if (inherits(data, "Seurat")) {
        original_levels <- levels(data@meta.data[[group.by]])
      } else {
        original_levels <- levels(data[[group.by]])
      }
      valid_subset_levels <- intersect(original_levels, idents)
      plot_df[[group.by.internal]] <- factor(plot_df[[group.by.internal]], 
                                             levels = valid_subset_levels)
    } else {
      plot_df[[group.by.internal]] <- as.factor(plot_df[[group.by.internal]])
    }
  } else {
    # Convert to factors (only if not already converted in idents filtering)
    plot_df[[group.by.internal]] <- as.factor(plot_df[[group.by.internal]])
  }
  
  plot_df[[split.by]] <- as.factor(plot_df[[split.by]])
  
  # Create plots for each feature
  all_plots <- list()
  
  for (feature_name in features) {
    if (!feature_name %in% names(plot_df)) {
      warning("Feature '", feature_name, "' not found, skipping")
      next
    }
    
    # Create base plot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data[[split.by]],
      y = .data[[feature_name]],
      fill = .data[[split.by]]
    ))
    
    # Add violin if requested
    if (violin) {
      p <- p + ggplot2::geom_violin(
        alpha = 0.5,
        scale = "width",
        trim = FALSE,
        na.rm = TRUE
      )
    }
    
    # Add boxplot
    p <- p + ggplot2::geom_boxplot(
      na.rm = TRUE,
      outlier.shape = NA,
      alpha = 0.7
    )
    
    # Add jitter points
    if (pt.size > 0) {
      p <- p + ggplot2::geom_jitter(
        color = "black",
        size = pt.size,
        alpha = if (aggregate) 0.7 else 0.3,
        width = 0.2,
        height = 0,
        na.rm = TRUE
      )
    }
    
    # Add statistical comparisons
    if (add_stats) {
      if (requireNamespace("ggpubr", quietly = TRUE)) {
        p <- p + ggpubr::stat_compare_means(
          label = "p.format",
          ...
        )
      } else {
        warning("ggpubr not available, skipping statistical tests")
      }
    }
    
    # Add faceting and theme
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", group.by.internal))) +
      ggplot2::labs(
        title = feature_name,
        x = split.by,
        y = if (aggregate) paste("Average", feature_name) else feature_name
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    all_plots[[feature_name]] <- p
  }
  
  if (length(all_plots) == 0) {
    message("No plots created")
    return(invisible(NULL))
  }
  
  # Combine plots if multiple features
  if (length(all_plots) == 1) {
    return(all_plots[[1]])
  } else {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(all_plots, ncol = ncol))
    } else {
      warning("patchwork not available, returning first plot only")
      return(all_plots[[1]])
    }
  }
}

