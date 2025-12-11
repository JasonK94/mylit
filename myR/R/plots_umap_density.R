#' UMAP Density Visualization
#'
#' Generates kernel density estimates on UMAP (or any two-dimensional reduction)
#' to highlight cell frequency patterns. Supports facetting by metadata,
#' zebra/contour style rendering, and signed density-difference maps between two
#' identity groups.
#'
#' @param data Seurat object or data.frame that already contains embedding
#'   columns.
#' @param reduction Character. Reduction slot to use when `data` is a Seurat
#'   object (default: `"umap.scvi"`).
#' @param dims Integer or character vector of length 2 specifying which
#'   dimensions to plot (default: `c(1, 2)`).
#' @param split.by Optional metadata column used to facet density plots.
#' @param group.by Optional metadata column (e.g., clusters) used for overlay
#'   contours/labels.
#' @param style Rendering style. `"zebra"` mimics flow-cytometry contour plots,
#'   `"color"` uses a continuous heatmap (default: `c("zebra", "color")`).
#' @param ident.by Metadata column used to interpret `ident.1`/`ident.2`. When
#'   `NULL`, the function falls back to `split.by`, and finally to the active
#'   Seurat identities if `data` is a Seurat object.
#' @param ident.1,ident.2 Optional vectors of identity labels used to compute a
#'   signed density difference (ident.1 minus ident.2). When supplied,
#'   `split.by` facets are ignored.
#' @param palette Character vector or palette name passed to the fill scale.
#'   Recognised keywords: `"viridis"`, `"magma"`, `"inferno"`, `"plasma"`,
#'   `"cividis"`, `"turbo"`, `"hot"`.
#' @param contour_bins Number of contour bands (default: 12).
#' @param resolution Integer grid resolution used by `stat_density_2d`
#'   (default: 200). Values below 50 are automatically raised.
#' @param bandwidth Optional numeric (length 1 or 2) bandwidth passed to
#'   `MASS::kde2d` when computing signed differences.
#' @param metadata_df Optional data.frame with extra metadata to join by cell.
#' @param match_id Column name inside `metadata_df` that matches `.cell_id`. If
#'   `NULL`, the function looks for `.cell_id` or rownames.
#' @param show_points Logical. Overlay raw points (default: `FALSE`).
#' @param point_size,point_alpha Aesthetics for optional raw points.
#' @param facet_scales Passed to `facet_wrap()` when faceting (default: `"fixed"`).
#' @param overlay Logical. Draw per-cluster overlays using `group.by`.
#' @param overlay_type `"contour"` (default) for contour lines, `"grid"` for
#'   dominant-glyph overlays, or `"both"` to combine them.
#' @param overlay_bins Integer. Number of contour levels per cluster (default: 4).
#' @param overlay_grid_res Integer grid resolution for `"grid"` overlays.
#' @param overlay_grid_min_density Numeric in `[0, 1]` indicating the minimum
#'   relative density required to draw glyphs (default: 0.1).
#' @param overlay_palette Optional vector of colors for clusters. Defaults to
#'   `RColorBrewer::Set2`.
#' @param overlay_symbols Optional vector of shape codes for `"grid"` overlays.
#' @param overlay_label Logical. Place text labels for each cluster overlay.
#' @param overlay_label_stat `"centroid"` (default) or `"median"` to control label
#'   placement.
#' @param overlay_label_size Numeric text size for overlay labels.
#' @param overlay_symbol_size Size for `"grid"` glyphs.
#' @param overlay_symbol_alpha Alpha for `"grid"` glyphs.
#' @param return_data Logical. When `TRUE`, returns a list with the ggplot
#'   object and the plotting dataframe.
#'
#' @return A ggplot object (or list with `plot` and `data` when
#'   `return_data = TRUE`).
#'
#' @examples
#' \dontrun{
#' # Single density map
#' plot_umap_density(sobj, split.by = NULL, style = "color")
#'
#' # Facet by condition with zebra contours
#' plot_umap_density(sobj, split.by = "condition", style = "zebra")
#'
#' # Signed density difference between responders and non-responders
#' plot_umap_density(
#'   sobj,
#'   split.by = "condition",
#'   ident.by = "condition",
#'   ident.1 = "responder",
#'   ident.2 = "non_responder",
#'   style = "color"
#' )
#' }
#' @export
plot_umap_density <- function(data,
                              reduction = "umap.scvi",
                              dims = c(1, 2),
                              split.by = NULL,
                              group.by = NULL,
                              style = c("zebra", "color"),
                              ident.by = NULL,
                              ident.1 = NULL,
                              ident.2 = NULL,
                              palette = "inferno",
                              contour_bins = 12,
                              resolution = 200,
                              bandwidth = NULL,
                              metadata_df = NULL,
                              match_id = NULL,
                              show_points = FALSE,
                              point_size = 0.3,
                              point_alpha = 0.35,
                              facet_scales = c("fixed", "free"),
                              overlay = FALSE,
                              overlay_type = c("contour", "grid", "both"),
                              overlay_bins = 4,
                              overlay_grid_res = 40,
                              overlay_grid_min_density = 0.1,
                              overlay_palette = NULL,
                              overlay_symbols = c(3, 4, 8, 11, 13, 15),
                              overlay_label = FALSE,
                              overlay_label_stat = c("centroid", "median"),
                              overlay_label_size = 3,
                              overlay_symbol_size = 2.5,
                              overlay_symbol_alpha = 0.7,
                              return_data = FALSE) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  style <- match.arg(style)
  facet_scales <- match.arg(facet_scales)
  overlay_type <- match.arg(overlay_type)
  overlay_label_stat <- match.arg(overlay_label_stat)
  contour_bins <- max(3, as.integer(contour_bins))
  resolution <- max(50, as.integer(resolution))
  overlay_bins <- max(1, as.integer(overlay_bins))
  overlay_grid_res <- max(10, as.integer(overlay_grid_res))
  overlay_grid_min_density <- max(0, min(overlay_grid_min_density, 1))
  
  if (!is.null(bandwidth)) {
    if (!is.numeric(bandwidth) || length(bandwidth) > 2) {
      stop("`bandwidth` must be numeric of length 1 or 2.")
    }
    if (length(bandwidth) == 1) {
      bandwidth <- rep(bandwidth, 2)
    }
  }
  
  plot_df <- .collect_umap_density_data(
    data = data,
    reduction = reduction,
    dims = dims,
    split.by = split.by,
    group.by = group.by,
    ident.by = ident.by,
    metadata_df = metadata_df,
    match_id = match_id
  )
  
  if (nrow(plot_df) < 10) {
    stop("Not enough cells to estimate densities (need >= 10 rows).")
  }
  
  # Remove rows with NA coordinates
  plot_df <- plot_df[stats::complete.cases(plot_df[, c("dim1", "dim2")]), ]
  if (nrow(plot_df) < 10) {
    stop("Not enough complete coordinate pairs after removing NA values.")
  }
  
  if (!is.null(split.by)) {
    if (!split.by %in% names(plot_df)) {
      stop("`split.by` column '", split.by, "' was not found in metadata.")
    }
    plot_df$split_col <- as.factor(plot_df[[split.by]])
  }
  
  wants_difference <- !is.null(ident.1) || !is.null(ident.2)
  if (wants_difference) {
    if (is.null(ident.1) || is.null(ident.2)) {
      stop("Both ident.1 and ident.2 must be supplied to compute a difference map.")
    }
    if (!is.null(split.by)) {
      warning("`split.by` is ignored when ident.1/ident.2 are provided.")
    }
    
    ident_col <- ident.by %||% split.by
    if (is.null(ident_col)) {
      if (inherits(data, "Seurat")) {
        active_ids <- as.character(Seurat::Idents(data))
        plot_df$.active_ident <- active_ids[match(plot_df$.cell_id, names(active_ids))]
        ident_col <- ".active_ident"
      } else {
        stop("Provide `ident.by` (or `split.by`) when using ident.1/ident.2 with data.frames.")
      }
    }
    
    if (!ident_col %in% names(plot_df)) {
      stop("`ident.by` column '", ident_col, "' was not found.")
    }
    
    diff_plot <- .plot_umap_density_difference(
      plot_df = plot_df,
      ident_col = ident_col,
      ident.1 = ident.1,
      ident.2 = ident.2,
      style = style,
      palette = palette,
      contour_bins = contour_bins,
      resolution = resolution,
      bandwidth = bandwidth
    )
    
    if (overlay) {
      if (is.null(group.by)) {
        stop("`overlay = TRUE` requires `group.by` when using ident.1/ident.2.")
      }
      overlay_layers <- .build_overlay_layers(
        plot_df = plot_df,
        cluster_col = group.by,
        palette = overlay_palette,
        overlay_bins = overlay_bins,
        overlay_label = overlay_label,
        overlay_label_stat = overlay_label_stat,
        overlay_label_size = overlay_label_size,
        overlay_type = overlay_type,
        overlay_symbols = overlay_symbols,
        overlay_symbol_size = overlay_symbol_size,
        overlay_symbol_alpha = overlay_symbol_alpha,
        overlay_grid_res = overlay_grid_res,
        overlay_grid_min_density = overlay_grid_min_density,
        bandwidth = bandwidth
      )
      for (layer in overlay_layers) {
        diff_plot <- diff_plot + layer
      }
    }
    
    if (return_data) {
      return(list(plot = diff_plot, data = plot_df))
    }
    return(diff_plot)
  }
  
  density_plot <- .plot_umap_density_single(
    plot_df = plot_df,
    style = style,
    palette = palette,
    contour_bins = contour_bins,
    resolution = resolution,
    split.by = split.by,
    show_points = show_points,
    point_size = point_size,
    point_alpha = point_alpha,
    facet_scales = facet_scales
  )
  
  if (return_data) {
    return(list(plot = density_plot, data = plot_df))
  }
  
  return(density_plot)
}

.collect_umap_density_data <- function(data,
                                       reduction,
                                       dims,
                                       split.by,
                                       group.by,
                                       ident.by,
                                       metadata_df,
                                       match_id) {
  meta_cols <- unique(c(split.by, ident.by, group.by))
  meta_cols <- meta_cols[!is.null(meta_cols)]
  axis_labels <- c("dim1", "dim2")
  
  if (inherits(data, "Seurat")) {
    if (!reduction %in% Seurat::Reductions(data)) {
      stop("Reduction '", reduction, "' not found in Seurat object.")
    }
    embed <- Seurat::Embeddings(data, reduction = reduction)
    dim_idx <- .resolve_dim_indices(embed, dims)
    embed_df <- as.data.frame(embed[, dim_idx, drop = FALSE])
    axis_labels <- colnames(embed_df)
    names(embed_df) <- c("dim1", "dim2")
    embed_df$.cell_id <- rownames(embed_df)
    
    meta_df <- data@meta.data
    meta_df$.cell_id <- rownames(meta_df)
    if (length(meta_cols) > 0) {
      keep_cols <- unique(c(".cell_id", meta_cols))
      keep_cols <- keep_cols[keep_cols %in% colnames(meta_df)]
      if (length(keep_cols) > 1) {
        embed_df <- dplyr::left_join(
          embed_df,
          meta_df[, keep_cols, drop = FALSE],
          by = ".cell_id"
        )
      }
    }
    
    if (!is.null(metadata_df)) {
      embed_df <- .join_extra_metadata(embed_df, metadata_df, match_id)
    }
    
    attr(embed_df, "axis_labels") <- axis_labels
    attr(embed_df, "reduction") <- reduction
    return(embed_df)
  }
  
  if (is.data.frame(data)) {
    df <- data
    dims <- .resolve_df_dims(df, dims)
    
    coord_df <- data.frame(
      dim1 = df[[dims[1]]],
      dim2 = df[[dims[2]]],
      stringsAsFactors = FALSE
    )
    if (".cell_id" %in% names(df)) {
      coord_df$.cell_id <- as.character(df$.cell_id)
    } else if (!is.null(rownames(df))) {
      coord_df$.cell_id <- rownames(df)
    } else {
      coord_df$.cell_id <- paste0("cell_", seq_len(nrow(df)))
    }
    
    if (length(meta_cols) > 0) {
      meta_cols <- meta_cols[meta_cols %in% names(df)]
      if (length(meta_cols) > 0) {
        coord_df <- cbind(coord_df, df[, meta_cols, drop = FALSE])
      }
    }
    
    if (!is.null(metadata_df)) {
      coord_df <- .join_extra_metadata(coord_df, metadata_df, match_id)
    }
    
    attr(coord_df, "axis_labels") <- dims
    attr(coord_df, "reduction") <- "custom"
    return(coord_df)
  }
  
  stop("`data` must be a Seurat object or data.frame.")
}

.resolve_dim_indices <- function(embedding, dims) {
  if (length(dims) != 2) {
    stop("`dims` must specify exactly two axes.")
  }
  if (is.numeric(dims)) {
    return(as.integer(dims))
  }
  idx <- match(dims, colnames(embedding))
  if (any(is.na(idx))) {
    stop("Could not match dims: ", paste(dims, collapse = ", "))
  }
  return(idx)
}

.resolve_df_dims <- function(df, dims) {
  if (length(dims) != 2) {
    stop("`dims` must specify exactly two axes.")
  }
  if (is.numeric(dims)) {
    dims <- names(df)[dims]
  }
  missing_dims <- setdiff(dims, names(df))
  if (length(missing_dims) > 0) {
    stop("Columns not found in dataframe: ", paste(missing_dims, collapse = ", "))
  }
  return(dims)
}

.join_extra_metadata <- function(plot_df, metadata_df, match_id) {
  md <- metadata_df
  if (!is.null(match_id)) {
    if (!match_id %in% names(md)) {
      stop("match_id '", match_id, "' not found in metadata_df.")
    }
    md$.cell_id <- md[[match_id]]
  } else if (".cell_id" %in% names(md)) {
    md$.cell_id <- as.character(md$.cell_id)
  } else if (!is.null(rownames(md))) {
    md$.cell_id <- rownames(md)
  } else {
    stop("metadata_df must contain `.cell_id` or you must provide `match_id`.")
  }
  
  md <- md[!duplicated(md$.cell_id), , drop = FALSE]
  plot_df <- dplyr::left_join(plot_df, md, by = ".cell_id")
  return(plot_df)
}

.plot_umap_density_single <- function(plot_df,
                                      style,
                                      palette,
                                      contour_bins,
                                      resolution,
                                      split.by,
                                      show_points,
                                      point_size,
                                      point_alpha,
                                      facet_scales) {
  axis_labels <- attr(plot_df, "axis_labels")
  x_lab <- axis_labels[1] %||% "dim1"
  y_lab <- axis_labels[2] %||% "dim2"
  reduction_name <- attr(plot_df, "reduction") %||% "reduction"
  
  base <- ggplot2::ggplot(plot_df, ggplot2::aes(dim1, dim2))
  
  if (style == "color") {
    pal <- .resolve_gradient_palette(palette, n = 256)
    base <- base +
      ggplot2::stat_density_2d_filled(
        ggplot2::aes(fill = after_stat(as.numeric(level))),
        contour_var = "density",
        geom = "polygon",
        n = resolution,
        alpha = 0.95
      ) +
      ggplot2::geom_density_2d(
        color = "white",
        bins = contour_bins,
        size = 0.25,
        n = resolution
      ) +
      ggplot2::scale_fill_gradientn(
        colours = pal,
        name = "Density",
        limits = c(0, NA)
      )
  } else {
    pal <- .resolve_zebra_palette(palette, contour_bins)
    base <- base +
      ggplot2::stat_density_2d_filled(
        ggplot2::aes(fill = after_stat(as.numeric(level))),
        contour_var = "density",
        bins = contour_bins,
        geom = "polygon",
        n = resolution,
        alpha = 0.95
      ) +
      ggplot2::geom_density_2d(
        color = "#101010",
        bins = contour_bins,
        size = 0.35,
        n = resolution
      ) +
      ggplot2::scale_fill_stepsn(
        colours = pal,
        n.breaks = contour_bins + 1,
        name = "Density"
      )
  }
  
  if (show_points) {
    base <- base + ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(dim1, dim2),
      size = point_size,
      alpha = point_alpha,
      inherit.aes = FALSE,
      color = "black"
    )
  }
  
  if (!is.null(split.by) && "split_col" %in% names(plot_df)) {
    base <- base + ggplot2::facet_wrap(~split_col, scales = facet_scales)
  }
  
  base <- base +
    ggplot2::labs(
      x = paste0(reduction_name, " ", x_lab),
      y = paste0(reduction_name, " ", y_lab),
      title = "UMAP Density"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = "#f5f5f5", colour = "grey80"),
      legend.position = "right"
    )
  
  if (overlay) {
    if (is.null(group.by)) {
      stop("`overlay = TRUE` requires `group.by` to be specified.")
    }
    overlay_layers <- .build_overlay_layers(
      plot_df = plot_df,
      cluster_col = group.by,
      palette = overlay_palette,
      overlay_bins = overlay_bins,
      overlay_label = overlay_label,
      overlay_label_stat = overlay_label_stat,
      overlay_label_size = overlay_label_size,
      overlay_type = overlay_type,
      overlay_symbols = overlay_symbols,
      overlay_symbol_size = overlay_symbol_size,
      overlay_symbol_alpha = overlay_symbol_alpha,
      overlay_grid_res = overlay_grid_res,
      overlay_grid_min_density = overlay_grid_min_density,
      bandwidth = bandwidth
    )
    for (layer in overlay_layers) {
      base <- base + layer
    }
  }
  
  return(base)
}

.plot_umap_density_difference <- function(plot_df,
                                          ident_col,
                                          ident.1,
                                          ident.2,
                                          style,
                                          palette,
                                          contour_bins,
                                          resolution,
                                          bandwidth) {
  plot_df$.ident_tmp <- as.character(plot_df[[ident_col]])
  ident.1 <- as.character(ident.1)
  ident.2 <- as.character(ident.2)
  
  group1 <- plot_df[plot_df$.ident_tmp %in% ident.1, ]
  group2 <- plot_df[plot_df$.ident_tmp %in% ident.2, ]
  if (nrow(group1) < 10 || nrow(group2) < 10) {
    stop("Each identity group must contain at least 10 cells for KDE.")
  }
  
  if (is.null(bandwidth)) {
    kde1 <- MASS::kde2d(group1$dim1, group1$dim2, n = resolution)
    kde2 <- MASS::kde2d(group2$dim1, group2$dim2, n = resolution)
  } else {
    kde1 <- MASS::kde2d(group1$dim1, group1$dim2, n = resolution, h = bandwidth)
    kde2 <- MASS::kde2d(group2$dim1, group2$dim2, n = resolution, h = bandwidth)
  }
  
  diff_df <- expand.grid(
    dim1 = kde1$x,
    dim2 = kde1$y,
    KEEP.OUT.ATTRS = FALSE
  )
  diff_df$density_diff <- as.vector(kde1$z - kde2$z)
  diff_df$density_ident1 <- as.vector(kde1$z)
  diff_df$density_ident2 <- as.vector(kde2$z)
  
  axis_labels <- attr(plot_df, "axis_labels")
  x_lab <- axis_labels[1] %||% "dim1"
  y_lab <- axis_labels[2] %||% "dim2"
  reduction_name <- attr(plot_df, "reduction") %||% "reduction"
  diff_label <- paste0(
    paste(ident.1, collapse = ", "),
    " - ",
    paste(ident.2, collapse = ", ")
  )
  
  fill_scale <- if (style == "color") {
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "#f7f7f7",
      high = "#b2182b",
      midpoint = 0,
      name = diff_label
    )
  } else {
    ggplot2::scale_fill_steps2(
      low = "#053061",
      mid = "#fdfdfd",
      high = "#67001f",
      name = diff_label,
      n.breaks = contour_bins + 2
    )
  }
  
  ggplot2::ggplot(diff_df, ggplot2::aes(dim1, dim2, fill = density_diff)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = density_diff),
      color = "white",
      bins = contour_bins,
      size = 0.3
    ) +
    fill_scale +
    ggplot2::labs(
      x = paste0(reduction_name, " ", x_lab),
      y = paste0(reduction_name, " ", y_lab),
      title = "UMAP Density Difference"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

.build_overlay_layers <- function(plot_df,
                                  cluster_col,
                                  palette = NULL,
                                  overlay_bins = 4,
                                  overlay_label = FALSE,
                                  overlay_label_stat = "centroid",
                                  overlay_label_size = 3,
                                  overlay_type = "contour",
                                  overlay_symbols = c(3, 4, 8, 11, 13, 15),
                                  overlay_symbol_size = 2.5,
                                  overlay_symbol_alpha = 0.7,
                                  overlay_grid_res = 40,
                                  overlay_grid_min_density = 0.1,
                                  bandwidth = NULL) {
  if (!cluster_col %in% names(plot_df)) {
    stop("`group.by` column '", cluster_col, "' not found for overlay.")
  }
  
  overlay_df <- plot_df[!is.na(plot_df[[cluster_col]]), ]
  if (nrow(overlay_df) == 0) {
    stop("No valid rows for overlay using '", cluster_col, "'.")
  }
  
  overlay_df$.overlay_cluster <- droplevels(as.factor(overlay_df[[cluster_col]]))
  cluster_levels <- levels(overlay_df$.overlay_cluster)
  n_clusters <- length(cluster_levels)
  
  use_contour <- overlay_type %in% c("contour", "both")
  use_grid <- overlay_type %in% c("grid", "both")
  
  palette_vec <- .resolve_overlay_palette(palette, cluster_levels)
  shape_vec <- .resolve_overlay_symbols(overlay_symbols, cluster_levels)
  
  layers <- list()
  color_needed <- FALSE
  shape_needed <- FALSE
  
  if (use_contour) {
    # Loop over clusters to create separate layers
    # This avoids potential issues with ggplot dropping groups in a single layer
    for (cl in cluster_levels) {
      cl_df <- overlay_df[overlay_df$.overlay_cluster == cl, ]
      if (nrow(cl_df) < 5) next 
      
      layers <- c(layers, list(
        ggplot2::geom_density_2d(
          data = cl_df,
          ggplot2::aes(
            x = dim1,
            y = dim2,
            color = .overlay_cluster
          ),
          bins = overlay_bins,
          size = 0.45,
          inherit.aes = FALSE,
          show.legend = NA # Allow merging into color scale
        )
      ))
      color_needed <- TRUE
    }
  }
  
  if (overlay_label) {
    label_df <- overlay_df %>%
      dplyr::group_by(.overlay_cluster) %>%
      dplyr::summarise(
        dim1 = if (overlay_label_stat == "median") stats::median(dim1, na.rm = TRUE) else mean(dim1, na.rm = TRUE),
        dim2 = if (overlay_label_stat == "median") stats::median(dim2, na.rm = TRUE) else mean(dim2, na.rm = TRUE),
        .groups = "drop"
      )
    if (nrow(label_df) > 0) {
      layers <- c(layers, list(
        ggplot2::geom_text(
          data = label_df,
          ggplot2::aes(
            x = dim1,
            y = dim2,
            label = .overlay_cluster,
            color = .overlay_cluster
          ),
          size = overlay_label_size,
          fontface = "bold",
          inherit.aes = FALSE
        )
      ))
      color_needed <- TRUE
    }
  }
  
  if (use_grid) {
    grid_df <- .build_overlay_grid_points(
      overlay_df = overlay_df,
      cluster_levels = cluster_levels,
      grid_res = overlay_grid_res,
      min_density = overlay_grid_min_density,
      bandwidth = bandwidth
    )
    
    if (!is.null(grid_df) && nrow(grid_df) > 0) {
      layers <- c(layers, list(
        ggplot2::geom_point(
          data = grid_df,
          ggplot2::aes(
            x = dim1,
            y = dim2,
            color = .overlay_cluster,
            shape = .overlay_cluster
          ),
          size = overlay_symbol_size,
          alpha = overlay_symbol_alpha,
          inherit.aes = FALSE
        )
      ))
      color_needed <- TRUE
      shape_needed <- TRUE
    }
  }
  
  if (color_needed) {
    layers <- c(layers, list(
      ggplot2::scale_color_manual(
        values = palette_vec,
        name = cluster_col
      )
    ))
  }
  
  if (shape_needed) {
    layers <- c(layers, list(
      ggplot2::scale_shape_manual(
        values = shape_vec,
        name = cluster_col
      )
    ))
  }
  
  return(layers)
}

.build_overlay_grid_points <- function(overlay_df,
                                       cluster_levels,
                                       grid_res,
                                       min_density,
                                       bandwidth) {
  dim1_range <- range(overlay_df$dim1, na.rm = TRUE)
  dim2_range <- range(overlay_df$dim2, na.rm = TRUE)
  dens_list <- list()
  
  for (cl in cluster_levels) {
    cl_df <- overlay_df[overlay_df$.overlay_cluster == cl, , drop = FALSE]
    if (nrow(cl_df) < 5) {
      next
    }
    if (is.null(bandwidth)) {
      kd <- MASS::kde2d(
        cl_df$dim1,
        cl_df$dim2,
        n = grid_res,
        lims = c(dim1_range, dim2_range)
      )
    } else {
      kd <- MASS::kde2d(
        cl_df$dim1,
        cl_df$dim2,
        n = grid_res,
        lims = c(dim1_range, dim2_range),
        h = bandwidth
      )
    }
    
    grid_points <- expand.grid(
      dim1 = kd$x,
      dim2 = kd$y,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid_points$density <- as.vector(kd$z)
    grid_points$.overlay_cluster <- cl
    dens_list[[cl]] <- grid_points
  }
  
  if (length(dens_list) == 0) {
    return(NULL)
  }
  
  dens_df <- dplyr::bind_rows(dens_list)
  if (nrow(dens_df) == 0) {
    return(NULL)
  }
  
  dominant <- dens_df %>%
    dplyr::group_by(dim1, dim2) %>%
    dplyr::slice_max(order_by = density, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  max_density <- max(dominant$density, na.rm = TRUE)
  if (!is.finite(max_density) || max_density <= 0) {
    return(NULL)
  }
  
  threshold <- max_density * min_density
  dominant <- dominant[dominant$density >= threshold, , drop = FALSE]
  
  if (nrow(dominant) > 5000) {
    dominant <- dominant[seq(1, nrow(dominant), length.out = 5000), ]
  }
  
  return(dominant)
}

.resolve_overlay_palette <- function(palette, cluster_levels) {
  n_clusters <- length(cluster_levels)
  
  if (is.null(palette)) {
    base_pal <- RColorBrewer::brewer.pal(min(max(3, n_clusters), 8), "Set2")
    if (length(base_pal) < n_clusters) {
      palette <- grDevices::colorRampPalette(base_pal)(n_clusters)
    } else {
      palette <- base_pal
    }
  } else {
    if (length(palette) < n_clusters) {
      palette <- grDevices::colorRampPalette(palette)(n_clusters)
    }
  }
  
  if (is.null(names(palette))) {
    names(palette) <- cluster_levels
  }
  palette <- palette[cluster_levels]
  names(palette) <- cluster_levels
  return(palette)
}

.resolve_overlay_symbols <- function(symbols, cluster_levels) {
  n_clusters <- length(cluster_levels)
  if (length(symbols) < n_clusters) {
    symbols <- rep(symbols, length.out = n_clusters)
  }
  names(symbols) <- cluster_levels
  return(symbols[cluster_levels])
}

.resolve_gradient_palette <- function(palette, n = 256) {
  known_opts <- c("viridis", "magma", "inferno", "plasma", "cividis", "turbo")
  if (length(palette) == 1) {
    if (palette %in% known_opts) {
      return(viridisLite::viridis(n, option = palette))
    }
    if (palette == "hot") {
      return(grDevices::heat.colors(n))
    }
    # Fallback: build a gradient from white to the provided color
    return(grDevices::colorRampPalette(c("#fdfdfd", palette))(n))
  }
  if (length(palette) < n) {
    return(grDevices::colorRampPalette(palette)(n))
  }
  return(palette)
}

.resolve_zebra_palette <- function(palette, bins) {
  if (length(palette) == 1) {
    key <- tolower(palette)
    viridis_opts <- c("viridis", "magma", "inferno", "plasma", "cividis", "turbo")
    
    if (key %in% c("bw", "grayscale", "zebra")) {
      palette <- c("#fdfdfd", "#e5e5e5", "#cacaca", "#9d9d9d", "#6f6f6f", "#303030")
    } else if (key %in% viridis_opts) {
      palette <- viridisLite::viridis(bins, option = key)
    } else if (.is_valid_color(palette)) {
      palette <- grDevices::colorRampPalette(c("#ffffff", palette))(bins)
    } else {
      palette <- c("#fdfdfd", "#e5e5e5", "#cacaca", "#9d9d9d", "#6f6f6f", "#303030")
    }
  }
  if (length(palette) < bins) {
    palette <- rep(palette, length.out = bins)
  }
  return(palette)
}

.is_valid_color <- function(col) {
  tryCatch({
    grDevices::col2rgb(col)
    TRUE
  }, error = function(...) FALSE)
}

