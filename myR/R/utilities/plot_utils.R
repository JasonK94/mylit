#' Plot Utility Functions
#'
#' This module provides utilities for plot creation, saving, and styling.
#'
#' @name plot_utils
NULL

#' Save Plot with Conflict Resolution
#'
#' Saves a plot to a file, automatically handling filename conflicts by 
#' appending a number to the filename if it already exists.
#'
#' @param plot_object ggplot2 object or other plot object
#' @param base_filename Base filename (with extension)
#' @param output_dir Directory to save plot to
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution in dots per inch
#' @param device Device to use for saving (auto-detected from extension if NULL)
#'
#' @return The full path to the saved file
#'
#' @examples
#' \dontrun{
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
#' save_plot_with_conflict_resolution(
#'   p, 
#'   "my_plot.png", 
#'   "plots", 
#'   width = 7, 
#'   height = 5, 
#'   dpi = 300
#' )
#' }
#'
#' @export
save_plot_with_conflict_resolution <- function(plot_object, 
                                               base_filename, 
                                               output_dir, 
                                               width, 
                                               height, 
                                               dpi,
                                               device = NULL) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Output directory created: ", output_dir)
  }
  
  # Parse filename
  filename_no_ext <- tools::file_path_sans_ext(base_filename)
  ext <- tools::file_ext(base_filename)
  if (ext == "") ext <- "png" # Default extension
  
  # Find available filename
  full_path <- file.path(output_dir, paste0(filename_no_ext, ".", ext))
  counter <- 1
  while (file.exists(full_path)) {
    full_path <- file.path(output_dir, paste0(filename_no_ext, "_", counter, ".", ext))
    counter <- counter + 1
  }
  
  # Save plot
  ggplot2::ggsave(
    filename = full_path, 
    plot = plot_object, 
    width = width, 
    height = height, 
    dpi = dpi,
    device = device
  )
  
  message("Plot saved: ", full_path)
  return(full_path)
}

#' Get Default Plot Dimensions
#'
#' Returns sensible default plot dimensions based on plot type.
#'
#' @param plot_type Type of plot ("heatmap", "scatter", "boxplot", "default")
#' @param n_items Number of items (e.g., genes, clusters) - used for sizing
#'
#' @return List with width and height in inches
#'
#' @examples
#' get_default_plot_dims("heatmap", n_items = 20)
#'
#' @export
get_default_plot_dims <- function(plot_type = "default", n_items = NULL) {
  dims <- switch(
    plot_type,
    "heatmap" = list(
      width = 7 + min(n_items %||% 0, 10) * 0.3,
      height = 5 + min(n_items %||% 0, 20) * 0.2
    ),
    "scatter" = list(width = 7, height = 5),
    "boxplot" = list(
      width = 5 + min(n_items %||% 0, 10) * 0.4,
      height = 6
    ),
    "default" = list(width = 7, height = 5)
  )
  
  # Cap maximum dimensions
  dims$width <- min(dims$width, 20)
  dims$height <- min(dims$height, 20)
  
  return(dims)
}

#' Sanitize Filename
#'
#' Removes or replaces characters that are problematic in filenames.
#'
#' @param filename Filename to sanitize
#' @param replacement Character to replace problematic characters with (default: "_")
#'
#' @return Sanitized filename
#'
#' @examples
#' sanitize_filename("my/file:name?.txt")
#' # Returns: "my_file_name_.txt"
#'
#' @export
sanitize_filename <- function(filename, replacement = "_") {
  # Remove or replace problematic characters
  sanitized <- gsub("[^a-zA-Z0-9_.-]", replacement, filename)
  # Remove multiple consecutive replacements
  sanitized <- gsub(paste0(replacement, "+"), replacement, sanitized)
  return(sanitized)
}

#' Create Plot Grid Layout
#'
#' Helper to calculate optimal grid layout dimensions.
#'
#' @param n_plots Number of plots to arrange
#' @param ncol Desired number of columns (if NULL, calculated automatically)
#' @param nrow Desired number of rows (if NULL, calculated automatically)
#'
#' @return List with ncol and nrow
#'
#' @examples
#' create_plot_grid(7)  # Returns list(ncol=3, nrow=3)
#' create_plot_grid(7, ncol=2)  # Returns list(ncol=2, nrow=4)
#'
#' @export
create_plot_grid <- function(n_plots, ncol = NULL, nrow = NULL) {
  if (!is.null(ncol) && !is.null(nrow)) {
    # Both specified
    if (ncol * nrow < n_plots) {
      warning("Grid ", ncol, "x", nrow, " too small for ", n_plots, " plots. ",
              "Adjusting dimensions.")
      ncol <- NULL
      nrow <- NULL
    }
  }
  
  if (is.null(ncol) && is.null(nrow)) {
    # Calculate square-ish grid
    ncol <- ceiling(sqrt(n_plots))
    nrow <- ceiling(n_plots / ncol)
  } else if (is.null(ncol)) {
    # Calculate ncol from nrow
    ncol <- ceiling(n_plots / nrow)
  } else if (is.null(nrow)) {
    # Calculate nrow from ncol
    nrow <- ceiling(n_plots / ncol)
  }
  
  return(list(ncol = ncol, nrow = nrow))
}

#' Add Significance Stars
#'
#' Converts p-values to significance stars.
#'
#' @param p_values Numeric vector of p-values
#' @param thresholds Named vector of thresholds (default: c("***"=0.001, "**"=0.01, "*"=0.05))
#' @param ns_symbol Symbol for non-significant (default: "ns")
#'
#' @return Character vector of significance symbols
#'
#' @examples
#' add_significance_stars(c(0.0001, 0.02, 0.3))
#' # Returns: c("***", "*", "ns")
#'
#' @export
add_significance_stars <- function(p_values, 
                                   thresholds = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                                   ns_symbol = "ns") {
  
  stars <- character(length(p_values))
  
  for (i in seq_along(p_values)) {
    if (is.na(p_values[i])) {
      stars[i] <- NA_character_
      next
    }
    
    assigned <- FALSE
    for (star_name in names(thresholds)) {
      if (p_values[i] < thresholds[star_name]) {
        stars[i] <- star_name
        assigned <- TRUE
        break
      }
    }
    
    if (!assigned) {
      stars[i] <- ns_symbol
    }
  }
  
  return(stars)
}

#' Get Color Palette
#'
#' Returns a color palette for consistent styling across plots.
#'
#' @param name Palette name ("default", "viridis", "Set1", "custom")
#' @param n Number of colors needed
#'
#' @return Vector of color codes
#'
#' @export
get_color_palette <- function(name = "default", n = NULL) {
  palette <- switch(
    name,
    "default" = RColorBrewer::brewer.pal(min(n %||% 8, 8), "Set2"),
    "viridis" = viridisLite::viridis(n %||% 8),
    "Set1" = RColorBrewer::brewer.pal(min(n %||% 9, 9), "Set1"),
    "custom" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                 "#FF7F00", "#FFFF33", "#A65628", "#F781BF"),
    RColorBrewer::brewer.pal(min(n %||% 8, 8), "Set2")
  )
  
  # Extend palette if needed
  if (!is.null(n) && length(palette) < n) {
    palette <- grDevices::colorRampPalette(palette)(n)
  }
  
  return(palette)
}


