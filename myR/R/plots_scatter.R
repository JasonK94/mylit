#' Scatter Plot Functions
#'
#' Standardized scatter plot functions for expression vs covariate visualization
#' with flexible aggregation and styling options.
#'
#' @name plots_scatter
NULL

#' Scatter Plot: Expression vs Covariate
#'
#' Creates scatter plots of gene expression vs numeric covariates with optional
#' aggregation, coloring, and smoothing. Supports both Seurat objects and dataframes.
#'
#' @param data Seurat object or data.frame
#' @param feature Character. Gene/feature name to plot
#' @param x_var Character. Numeric covariate column name
#' @param group.by Character. Column to aggregate by (default: "sample_no")
#' @param split.by Character. Optional categorical variable to split/color by
#' @param color_by Character. Optional column to color points by
#' @param aggregate Logical. Whether to aggregate cells by group.by (default: TRUE)
#' @param transpose Logical. Swap x and y axes (default: FALSE)
#' @param fitted_line Character. Smoothing method: "linear", "loess", "lasso", or NULL
#' @param palette Character vector or palette name for colors
#' @param transparency Logical. Map numeric color_by to alpha (default: TRUE)
#' @param transparency_desc Logical. Inverse alpha mapping (default: FALSE)
#' @param assay Character. Assay to use (Seurat only)
#' @param layer Character. Layer to use: "counts", "data", "scale.data" (default: "data")
#' @param metadata_df Optional data.frame. Additional metadata to join
#' @param match_id Character. Column name in metadata_df for matching
#' @param point_size Numeric. Point size (default: 3)
#' @param title Character. Plot title (auto-generated if NULL)
#' @param xlab Character. X-axis label (auto-generated if NULL)
#' @param ylab Character. Y-axis label (auto-generated if NULL)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#' # Basic scatter plot
#' p <- plot_scatter(sobj, feature = "CD3D", x_var = "nih_change", group.by = "sample")
#'
#' # With split.by
#' p <- plot_scatter(sobj, feature = "CD3D", x_var = "nih_change", 
#'                   split.by = "condition", fitted_line = "linear")
#'
#' # With color_by
#' p <- plot_scatter(sobj, feature = "CD3D", x_var = "nih_change",
#'                   color_by = "age", palette = "viridis")
#' }
#'
#' @export
plot_scatter <- function(data,
                         feature,
                         x_var,
                         group.by = "sample_no",
                         split.by = NULL,
                         color_by = NULL,
                         aggregate = TRUE,
                         aggregate_by = NULL,
                         transpose = FALSE,
                         fitted_line = c("linear", "loess", "lasso"),
                         palette = NULL,
                         transparency = TRUE,
                         transparency_desc = FALSE,
                         assay = NULL,
                         layer = "data",
                         metadata_df = NULL,
                         match_id = NULL,
                         point_size = 3,
                         title = NULL,
                         xlab = NULL,
                         ylab = NULL) {
  
  # Handle fitted_line (can be NULL)
  if (!is.null(fitted_line)) {
    fitted_line <- match.arg(fitted_line, choices = c("linear", "loess", "lasso"))
  }
  stopifnot(is.character(feature), length(feature) == 1)
  stopifnot(is.character(x_var), length(x_var) == 1)
  
  # Helper function for NULL coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Prepare data
  if (aggregate) {
    # Determine aggregation variables
    if (is.null(aggregate_by)) {
      # Include split.by in aggregation if provided
      aggregate_by <- unique(c(group.by, split.by))
      aggregate_by <- aggregate_by[!is.null(aggregate_by)]
    }
    
    plot_df <- .prepare_data_with_aggregation(
      data = data,
      features = c(feature, x_var),
      group.by = group.by,
      split.by = split.by,
      sample_col = group.by,
      aggregate = TRUE,
      aggregate_by = aggregate_by,
      agg_fun = "mean",
      assay = assay,
      layer = layer,
      metadata_df = metadata_df,
      match_id = match_id
    )
    
    # Rename aggregated columns
    if (feature %in% names(plot_df)) {
      plot_df$avg_expr <- plot_df[[feature]]
    } else {
      stop("Feature '", feature, "' not found after aggregation")
    }
    
    if (x_var %in% names(plot_df)) {
      plot_df$x_val <- plot_df[[x_var]]
    } else {
      stop("Covariate '", x_var, "' not found after aggregation")
    }
    
    # Handle color_by and split.by in aggregated data
    if (!is.null(color_by) && color_by %in% names(plot_df)) {
      plot_df$colour <- plot_df[[color_by]]
    }
    if (!is.null(split.by) && split.by %in% names(plot_df)) {
      plot_df$split_col <- plot_df[[split.by]]
    }
    
  } else {
    # Cell-level data
    plot_df <- .prepare_plot_data(
      data = data,
      features = c(feature, x_var),
      group.by = group.by,
      split.by = split.by,
      assay = assay,
      layer = layer,
      metadata_df = metadata_df,
      match_id = match_id
    )
    
    # Aggregate by group.by for plotting
    is_color_numeric <- if (!is.null(color_by)) {
      is.numeric(plot_df[[color_by]])
    } else {
      FALSE
    }
    
    agg_df <- plot_df %>%
      dplyr::group_by(.data[[group.by]]) %>%
      dplyr::summarise(
        avg_expr = mean(.data[[feature]], na.rm = TRUE),
        x_val = mean(.data[[x_var]], na.rm = TRUE),
        colour = if (!is.null(color_by)) {
          if (is_color_numeric) {
            mean(.data[[color_by]], na.rm = TRUE)
          } else {
            dplyr::first(.data[[color_by]], na_rm = TRUE)
          }
        } else {
          NA
        },
        split_col = if (!is.null(split.by)) {
          dplyr::first(.data[[split.by]], na_rm = TRUE)
        } else {
          NA
        },
        .groups = "drop"
      )
    
    plot_df <- agg_df
  }
  
  # Set axis columns
  x_col <- if (transpose) "avg_expr" else "x_val"
  y_col <- if (transpose) "x_val" else "avg_expr"
  
  # Create base plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
  
  # Add points with coloring
  if (!is.null(split.by)) {
    if (!is.null(color_by)) {
      warning("`color_by` is ignored when `split.by` is provided")
    }
    
    if (!"split_col" %in% names(plot_df)) {
      stop("split.by column '", split.by, "' not found in aggregated data")
    }
    
    # Convert to factor, handling numeric values
    if (is.numeric(plot_df$split_col)) {
      plot_df$split_col <- as.factor(as.character(plot_df$split_col))
    } else {
      plot_df$split_col <- as.factor(plot_df$split_col)
    }
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = .data[["split_col"]]),
      size = point_size
    )
    
    # Set palette
    n_levels <- length(unique(plot_df$split_col))
    if (is.null(palette)) {
      # Use Set1 for small number of levels, otherwise use a larger palette
      if (n_levels <= 9) {
        pal <- RColorBrewer::brewer.pal(max(3, n_levels), "Set1")[1:n_levels]
      } else {
        pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_levels)
      }
    } else {
      pal <- palette
      if (length(pal) < n_levels) {
        pal <- grDevices::colorRampPalette(pal)(n_levels)
      }
    }
    p <- p + ggplot2::scale_colour_manual(values = pal, name = split.by)
    
  } else if (!is.null(color_by)) {
    if (!"colour" %in% names(plot_df)) {
      stop("color_by column '", color_by, "' not found in aggregated data")
    }
    
    is_color_numeric <- is.numeric(plot_df$colour)
    
    if (is_color_numeric) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(colour = colour, alpha = if(transparency) colour else NULL),
        size = point_size
      )
      
      if (transparency) {
        alpha_range <- if (transparency_desc) c(1, 0.2) else c(0.2, 1)
        p <- p + ggplot2::scale_alpha(range = alpha_range, guide = "none")
      }
      
      pal <- if (is.null(palette)) {
        viridisLite::viridis(256)
      } else if (is.character(palette) && length(palette) == 1 &&
                 palette %in% c("magma", "plasma", "inferno", "viridis", "cividis", "turbo")) {
        viridisLite::viridis(256, option = palette)
      } else {
        palette
      }
      p <- p + ggplot2::scale_colour_gradientn(colours = pal, name = color_by)
      
    } else {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(colour = colour),
        size = point_size
      )
      
      n_levels <- length(unique(plot_df$colour))
      pal <- if (is.null(palette)) RColorBrewer::brewer.pal(max(3, n_levels), "Set1") else palette
      if (length(pal) < n_levels) {
        pal <- grDevices::colorRampPalette(pal)(n_levels)
      }
      p <- p + ggplot2::scale_colour_manual(values = pal, name = color_by)
    }
    
  } else {
    p <- p + ggplot2::geom_point(size = point_size)
  }
  
  # Add smoothing line
  if (!is.null(fitted_line)) {
    method_val <- switch(
      fitted_line,
      "linear" = "lm",
      "loess" = "loess",
      "lasso" = NULL,  # Special handling
      NULL
    )
    
    if (!is.null(split.by)) {
      # Split by smoothing
      if (fitted_line == "lasso") {
        warning("Lasso smoothing not supported with split.by, using linear")
        method_val <- "lm"
      }
      
      p <- p + ggplot2::geom_smooth(
        ggplot2::aes(colour = .data[["split_col"]]),
        method = method_val,
        se = TRUE,
        show.legend = FALSE
      )
      
    } else {
      # Single smoothing line
      if (fitted_line == "linear") {
        p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black")
        
        # Add statistics annotation
        fit <- tryCatch(
          stats::lm(plot_df[[y_col]] ~ plot_df[[x_col]]),
          error = function(e) NULL
        )
        if (!is.null(fit) && nrow(summary(fit)$coefficients) >= 2) {
          coef <- round(stats::coef(fit), 3)
          pval <- signif(summary(fit)$coefficients[2, 4], 3)
          annot <- paste0("y = ", coef[1], " + ", coef[2], " * x\np = ", pval)
          p <- p + ggplot2::annotate(
            "text",
            x = min(plot_df[[x_col]], na.rm = TRUE),
            y = max(plot_df[[y_col]], na.rm = TRUE),
            label = annot,
            hjust = 0,
            vjust = 1,
            size = 4
          )
        }
        
      } else if (fitted_line == "loess") {
        p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "black")
        
      } else if (fitted_line == "lasso") {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
          warning("glmnet not installed; falling back to linear fit")
          p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black")
        } else {
          xmat <- as.matrix(plot_df[[x_col]])
          fit <- glmnet::cv.glmnet(xmat, plot_df[[y_col]], alpha = 1)
          preds <- as.numeric(glmnet::predict.glmnet(
            fit$glmnet.fit,
            newx = xmat,
            s = fit$lambda.min
          ))
          pred_df <- plot_df %>% dplyr::mutate(pred = preds)
          p <- p + ggplot2::geom_line(
            data = pred_df[order(pred_df[[x_col]]), ],
            ggplot2::aes(x = .data[[x_col]], y = pred),
            colour = "red",
            linewidth = 1
          )
        }
      }
    }
  }
  
  # Set labels
  x_label <- if (is.null(xlab)) {
    if (transpose) paste("Average", feature, "expression") else x_var
  } else {
    xlab
  }
  
  y_label <- if (is.null(ylab)) {
    if (transpose) x_var else paste("Average", feature, "expression")
  } else {
    ylab
  }
  
  plot_title <- if (is.null(title)) {
    if (transpose) paste("Average", feature, "expression vs", x_var) else paste(x_var, "vs Average", feature, "expression")
  } else {
    title
  }
  
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = plot_title
    )
  
  return(p)
}

