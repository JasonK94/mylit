#' Basic Plotting Functions
#'
#' This module provides basic plotting utilities for single vectors and simple distributions.
#'
#' @name basic_plots
NULL

#' Create Histogram for Numeric Vector
#'
#' Generates a histogram for a numeric vector or a column in a data frame.
#'
#' @param x Numeric vector or data frame
#' @param col If x is data frame, column name or index
#' @param bins Number of bins (default: 30)
#' @param title Plot title
#' @param xlab X-axis label
#' @param color Fill color (default: "steelblue")
#'
#' @return ggplot2 object
#'
#' @export
mybar <- function(x, col = NULL, bins = 30, title = NULL, xlab = NULL, color = "steelblue") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  # Extract data
  if (is.data.frame(x)) {
    if (is.null(col)) {
      stop("col must be specified when x is a data frame")
    }
    data_vec <- x[[col]]
    if (is.null(xlab)) xlab <- as.character(col)
  } else {
    data_vec <- x
    if (is.null(xlab)) xlab <- "Value"
  }
  
  if (!is.numeric(data_vec)) {
    stop("Data must be numeric")
  }
  
  # Create plot
  p <- ggplot2::ggplot(data.frame(x = data_vec), ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(bins = bins, fill = color, color = "black", alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% "Histogram",
      x = xlab,
      y = "Count"
    )
  
  return(p)
}

#' Create Density Plot for Numeric Vector
#'
#' Generates a density plot for a numeric vector or a column in a data frame.
#'
#' @param x Numeric vector or data frame
#' @param col If x is data frame, column name or index
#' @param title Plot title
#' @param xlab X-axis label
#' @param color Fill color (default: "steelblue")
#' @param alpha Transparency (default: 0.5)
#'
#' @return ggplot2 object
#'
#' @export
mydensity <- function(x, col = NULL, title = NULL, xlab = NULL, 
                      color = "steelblue", alpha = 0.5) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  # Extract data
  if (is.data.frame(x)) {
    if (is.null(col)) {
      stop("col must be specified when x is a data frame")
    }
    data_vec <- x[[col]]
    if (is.null(xlab)) xlab <- as.character(col)
  } else {
    data_vec <- x
    if (is.null(xlab)) xlab <- "Value"
  }
  
  if (!is.numeric(data_vec)) {
    stop("Data must be numeric")
  }
  
  # Create plot
  p <- ggplot2::ggplot(data.frame(x = data_vec), ggplot2::aes(x = x)) +
    ggplot2::geom_density(fill = color, alpha = alpha, color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% "Density Plot",
      x = xlab,
      y = "Density"
    )
  
  return(p)
}

#' Create Single Vector Line Plot
#'
#' Creates a line plot for a numeric vector with optional cutoff line.
#'
#' @param vec Numeric vector
#' @param cutoff Optional horizontal cutoff line value
#' @param show_intersection Show intersection point with cutoff (default: TRUE)
#' @param title Plot title
#' @param xlab X-axis label (default: "Index")
#' @param ylab Y-axis label (default: "Value")
#' @param color Line color (default: "steelblue")
#'
#' @return ggplot2 object
#'
#' @export
myline <- function(vec, cutoff = NULL, show_intersection = TRUE,
                   title = NULL, xlab = "Index", ylab = "Value",
                   color = "steelblue") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  if (!is.numeric(vec)) {
    stop("vec must be numeric")
  }
  
  # Create data frame
  plot_data <- data.frame(
    index = seq_along(vec),
    value = vec
  )
  
  # Base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = index, y = value)) +
    ggplot2::geom_line(color = color, size = 1) +
    ggplot2::geom_point(color = color, size = 2, alpha = 0.6) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% "Line Plot",
      x = xlab,
      y = ylab
    )
  
  # Add cutoff line
  if (!is.null(cutoff)) {
    p <- p + ggplot2::geom_hline(yintercept = cutoff, 
                                  linetype = "dashed", 
                                  color = "red", 
                                  size = 1)
    
    # Find intersection
    if (show_intersection && any(vec > cutoff)) {
      first_above <- which(vec > cutoff)[1]
      if (!is.na(first_above)) {
        p <- p + ggplot2::annotate(
          "point",
          x = first_above,
          y = vec[first_above],
          color = "red",
          size = 4
        ) +
        ggplot2::annotate(
          "text",
          x = first_above,
          y = vec[first_above],
          label = paste0("n=", first_above),
          vjust = -1,
          color = "red"
        )
      }
    }
  }
  
  return(p)
}

#' Create Multiple Vector Line Plot with Dual Y-Axis
#'
#' Creates a line plot for multiple numeric vectors with optional dual y-axis.
#'
#' @param vec1 First numeric vector
#' @param vec2 Second numeric vector (optional)
#' @param cutoff1 Cutoff for vec1
#' @param cutoff2 Cutoff for vec2
#' @param title Plot title
#' @param ylab1 Y-axis label for vec1
#' @param ylab2 Y-axis label for vec2
#' @param color1 Color for vec1 (default: "steelblue")
#' @param color2 Color for vec2 (default: "darkred")
#'
#' @return ggplot2 object
#'
#' @export
mylines <- function(vec1, vec2 = NULL, cutoff1 = NULL, cutoff2 = NULL,
                    title = NULL, ylab1 = "Vector 1", ylab2 = "Vector 2",
                    color1 = "steelblue", color2 = "darkred") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  if (!is.numeric(vec1)) {
    stop("vec1 must be numeric")
  }
  
  # Create data frame
  plot_data <- data.frame(
    index = seq_along(vec1),
    vec1 = vec1
  )
  
  if (!is.null(vec2)) {
    if (length(vec2) != length(vec1)) {
      stop("vec1 and vec2 must have the same length")
    }
    plot_data$vec2 <- vec2
    
    # Use dual y-axis
    scaling_factor <- max(vec1, na.rm = TRUE) / max(vec2, na.rm = TRUE)
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = index)) +
      ggplot2::geom_line(ggplot2::aes(y = vec1, color = "Vector 1"), size = 1) +
      ggplot2::geom_line(ggplot2::aes(y = vec2 * scaling_factor, color = "Vector 2"), size = 1) +
      ggplot2::scale_y_continuous(
        name = ylab1,
        sec.axis = ggplot2::sec_axis(~./scaling_factor, name = ylab2)
      ) +
      ggplot2::scale_color_manual(values = c("Vector 1" = color1, "Vector 2" = color2)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title %||% "Multi-Line Plot",
        x = "Index",
        color = "Variable"
      )
    
    # Add cutoffs
    if (!is.null(cutoff1)) {
      p <- p + ggplot2::geom_hline(yintercept = cutoff1, 
                                    linetype = "dashed", 
                                    color = color1)
    }
    if (!is.null(cutoff2)) {
      p <- p + ggplot2::geom_hline(yintercept = cutoff2 * scaling_factor, 
                                    linetype = "dashed", 
                                    color = color2)
    }
  } else {
    # Single vector
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = index, y = vec1)) +
      ggplot2::geom_line(color = color1, size = 1) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title %||% "Line Plot",
        x = "Index",
        y = ylab1
      )
    
    if (!is.null(cutoff1)) {
      p <- p + ggplot2::geom_hline(yintercept = cutoff1, 
                                    linetype = "dashed", 
                                    color = "red")
    }
  }
  
  return(p)
}

#' Create Cumulative Distribution Function Plot
#'
#' Computes and plots the CDF for probability, logit, or ratio values.
#'
#' @param values Numeric vector of values
#' @param type Type of values: "prob", "logit", or "ratio" (default: "prob")
#' @param title Plot title
#' @param xlab X-axis label
#' @param color Line color (default: "steelblue")
#'
#' @return ggplot2 object
#'
#' @export
cdf <- function(values, type = c("prob", "logit", "ratio"), 
                title = NULL, xlab = NULL, color = "steelblue") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  type <- match.arg(type)
  
  # Transform if needed
  transformed_values <- switch(
    type,
    "prob" = values,
    "logit" = log(values / (1 - values)),
    "ratio" = log2(values)
  )
  
  # Compute CDF
  sorted_values <- sort(transformed_values)
  cdf_values <- seq_along(sorted_values) / length(sorted_values)
  
  plot_data <- data.frame(
    value = sorted_values,
    cdf = cdf_values
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = cdf)) +
    ggplot2::geom_line(color = color, size = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste("CDF -", toupper(type)),
      x = xlab %||% paste(toupper(type), "Value"),
      y = "Cumulative Probability"
    )
  
  return(p)
}

#' Create Multiple CDF Plots
#'
#' Computes and plots CDFs across multiple or grouped datasets.
#'
#' @param data_list List of numeric vectors or data frame with grouping column
#' @param group_col If data is data.frame, column name for grouping
#' @param value_col If data is data.frame, column name for values
#' @param type Type of values: "prob", "logit", or "ratio" (default: "prob")
#' @param title Plot title
#' @param xlab X-axis label
#'
#' @return ggplot2 object
#'
#' @export
cdf_multi <- function(data_list, group_col = NULL, value_col = NULL,
                      type = c("prob", "logit", "ratio"), 
                      title = NULL, xlab = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  type <- match.arg(type)
  
  # Prepare data
  if (is.data.frame(data_list)) {
    if (is.null(group_col) || is.null(value_col)) {
      stop("group_col and value_col must be specified when data_list is a data frame")
    }
    
    groups <- unique(data_list[[group_col]])
    plot_data_list <- list()
    
    for (grp in groups) {
      values <- data_list[data_list[[group_col]] == grp, value_col]
      
      # Transform
      transformed <- switch(
        type,
        "prob" = values,
        "logit" = log(values / (1 - values)),
        "ratio" = log2(values)
      )
      
      sorted_values <- sort(transformed)
      cdf_values <- seq_along(sorted_values) / length(sorted_values)
      
      plot_data_list[[as.character(grp)]] <- data.frame(
        value = sorted_values,
        cdf = cdf_values,
        group = as.character(grp)
      )
    }
    
    plot_data <- do.call(rbind, plot_data_list)
    
  } else if (is.list(data_list)) {
    plot_data_list <- list()
    
    for (i in seq_along(data_list)) {
      grp_name <- names(data_list)[i] %||% paste0("Group", i)
      values <- data_list[[i]]
      
      # Transform
      transformed <- switch(
        type,
        "prob" = values,
        "logit" = log(values / (1 - values)),
        "ratio" = log2(values)
      )
      
      sorted_values <- sort(transformed)
      cdf_values <- seq_along(sorted_values) / length(sorted_values)
      
      plot_data_list[[grp_name]] <- data.frame(
        value = sorted_values,
        cdf = cdf_values,
        group = grp_name
      )
    }
    
    plot_data <- do.call(rbind, plot_data_list)
    
  } else {
    stop("data_list must be a list or data frame")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value, y = cdf, color = group)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste("CDF Comparison -", toupper(type)),
      x = xlab %||% paste(toupper(type), "Value"),
      y = "Cumulative Probability",
      color = "Group"
    )
  
  return(p)
}

