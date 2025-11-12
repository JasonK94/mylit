#' Volcano Plot Functions
#'
#' Standardized volcano plot functions for visualizing effect sizes and p-values.
#'
#' @name plots_volcano
NULL

#' Volcano Plot
#'
#' Creates a volcano plot from statistical results (e.g., differential expression,
#' linear model results). Supports flexible column naming and filtering.
#'
#' @param data data.frame with statistical results
#' @param x_col Character. Column name for effect size/estimate (default: "estimate")
#' @param y_col Character. Column name for p-value (default: "p_value")
#' @param gene_col Character. Column name for gene/feature names (default: "gene")
#' @param filter_col Character. Column name to filter by pattern (default: "term")
#' @param filter_pattern Character. Regex pattern to filter rows
#' @param effect_threshold Numeric. Threshold for "large effect" (default: 0.5)
#' @param p_threshold Numeric. P-value significance threshold (default: 0.05)
#' @param p_adj_col Character. Column name for adjusted p-value (optional)
#' @param use_p_adj Logical. Use adjusted p-value if available (default: TRUE)
#' @param label_top Integer. Number of top genes to label (default: 10)
#' @param label_col Character. Column name for labels (default: gene_col)
#' @param resize_x Logical. Auto-resize x-axis limits (default: TRUE)
#' @param resize_y Logical. Auto-resize y-axis limits (default: FALSE)
#' @param title Character. Plot title (default: "Volcano Plot")
#' @param xlab Character. X-axis label (default: "Effect Size")
#' @param ylab Character. Y-axis label (default: "-log10(p-value)")
#' @param colors Named character vector. Custom colors for categories
#' @param point_size Numeric. Point size (default: 2)
#' @param point_alpha Numeric. Point transparency (default: 0.6)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#' # Basic volcano plot
#' p <- plot_volcano(lmm_results, x_col = "estimate", y_col = "p_value")
#'
#' # With filtering
#' p <- plot_volcano(lmm_results, filter_pattern = "drug.*timepoint", label_top = 20)
#'
#' # Using adjusted p-values
#' p <- plot_volcano(lmm_results, p_adj_col = "p_adj", use_p_adj = TRUE)
#' }
#'
#' @export
plot_volcano <- function(data,
                         x_col = "estimate",
                         y_col = "p_value",
                         gene_col = "gene",
                         filter_col = "term",
                         filter_pattern = NULL,
                         effect_threshold = 0.5,
                         p_threshold = 0.05,
                         p_adj_col = NULL,
                         use_p_adj = TRUE,
                         label_top = 10,
                         label_col = NULL,
                         resize_x = TRUE,
                         resize_y = FALSE,
                         title = "Volcano Plot",
                         xlab = NULL,
                         ylab = NULL,
                         colors = NULL,
                         point_size = 2,
                         point_alpha = 0.6) {
  
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  # Helper function for NULL coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  label_col <- if (is.null(label_col)) gene_col else label_col
  
  # Filter data if pattern provided
  plot_data <- data
  if (!is.null(filter_pattern) && !is.null(filter_col)) {
    if (!filter_col %in% names(plot_data)) {
      warning("filter_col '", filter_col, "' not found, skipping filter")
    } else {
      plot_data <- plot_data %>%
        dplyr::filter(grepl(filter_pattern, .data[[filter_col]], perl = TRUE))
    }
  }
  
  if (nrow(plot_data) == 0) {
    stop("No data remaining after filtering")
  }
  
  # Check required columns
  required_cols <- c(x_col, y_col)
  if (label_top > 0 && !label_col %in% names(plot_data)) {
    warning("label_col '", label_col, "' not found, skipping labels")
    label_top <- 0
  }
  
  missing_cols <- setdiff(required_cols, names(plot_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Determine which p-value to use
  p_col <- y_col
  if (use_p_adj && !is.null(p_adj_col) && p_adj_col %in% names(plot_data)) {
    p_col <- p_adj_col
  }
  
  # Prepare plotting variables
  plot_data <- plot_data %>%
    dplyr::mutate(
      effect_size = .data[[x_col]],
      p_value_plot = .data[[p_col]],
      log_p = -log10(.data[[p_col]]),
      category = dplyr::case_when(
        abs(.data[["effect_size"]]) > effect_threshold & 
          .data[["p_value_plot"]] < p_threshold ~ "Significant",
        abs(.data[["effect_size"]]) > effect_threshold ~ "Large effect only",
        .data[["p_value_plot"]] < p_threshold ~ "Significant only",
        TRUE ~ "Not significant"
      ),
      category = factor(
        category,
        levels = c("Significant", "Large effect only", 
                   "Significant only", "Not significant")
      )
    )
  
  # Set default colors
  default_colors <- c(
    "Significant" = "#d62728",
    "Large effect only" = "#ff7f0e",
    "Significant only" = "#1f77b4",
    "Not significant" = "gray70"
  )
  
  plot_colors <- if (is.null(colors)) default_colors else colors
  
  # Identify top genes to label
  top_genes <- NULL
  if (label_top > 0 && label_col %in% names(plot_data)) {
    top_genes <- plot_data %>%
      dplyr::filter(category == "Significant") %>%
      dplyr::arrange(.data[["p_value_plot"]]) %>%
      head(label_top)
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data[["effect_size"]],
    y = .data[["log_p"]],
    color = category
  )) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::geom_hline(
      yintercept = -log10(p_threshold),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-effect_threshold, effect_threshold),
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::scale_color_manual(values = plot_colors, name = "Category") +
    ggplot2::labs(
      title = title,
      x = if (is.null(xlab)) "Effect Size (Estimate)" else xlab,
      y = if (is.null(ylab)) "-log10(p-value)" else ylab
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  # Add gene labels
  if (!is.null(top_genes) && nrow(top_genes) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = top_genes,
        ggplot2::aes(label = .data[[label_col]]),
        size = 3,
        max.overlaps = 20,
        color = "black"
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = top_genes,
        ggplot2::aes(label = .data[[label_col]]),
        size = 3,
        color = "black",
        hjust = 0,
        vjust = 0
      )
    }
  }
  
  # Resize axes if requested
  if (resize_x) {
    xlim_val <- quantile(abs(plot_data$effect_size), 0.95, na.rm = TRUE)
    p <- p + ggplot2::coord_cartesian(xlim = c(-xlim_val, xlim_val))
  }
  
  if (resize_y) {
    ylim_val <- quantile(plot_data$log_p, 0.95, na.rm = TRUE)
    p <- p + ggplot2::coord_cartesian(ylim = c(0, ylim_val))
  }
  
  return(p)
}

