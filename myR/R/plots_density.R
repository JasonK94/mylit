#' Kernel Density Plot for Milo Neighborhoods / Single Cells
#'
#' Computes and plots a 2D kernel density estimate on a dimension reduction layout (e.g. UMAP),
#' optionally weighted by a numeric metric (e.g. logFC, significance).
#'
#' @param object A `Milo` object, `Seurat` object, or `SingleCellExperiment`.
#' @param reduction Character. Name of the reduction slot to use (default: "UMAP").
#' @param weight_metric Numeric vector or character string.
#'   - If character, it should be a column name in `da_results` (for Milo) or `meta.data` (for Seurat).
#'   - If numeric vector, it must match the number of points (neighborhoods or cells).
#'   - If NULL, standard density (cell/nhood density) is computed.
#' @param da_results Data frame. Required if `object` is Milo and `weight_metric` is a column name.
#' @param bandwidth Numeric vector of length 2. Smoothing bandwidth (sigma) for x and y directions.
#'   If NULL, estimated automatically via `MASS::bandwidth.nrd`.
#' @param grid_size Integer. Number of grid points in each direction (default: 100).
#' @param palette Character vector. Color palette for density (default: viridis).
#' @param alpha_range Numeric vector of length 2. Alpha range for density layer (default: c(0, 1)).
#' @param show_points Logical. Overlay original points? (default: FALSE)
#' @param point_size Numeric. Size of overlaid points (default: 0.5).
#' @param point_alpha Numeric. Transparency of overlaid points (default: 0.1).
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#'
#' @importFrom MASS kde2d
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn geom_point theme_minimal labs scale_alpha_continuous
#' @importFrom viridisLite viridis
#' @export
plot_milo_density <- function(
    object,
    reduction = "UMAP",
    weight_metric = NULL,
    da_results = NULL,
    bandwidth = NULL,
    grid_size = 100,
    palette = NULL,
    alpha_range = c(0, 1),
    show_points = FALSE,
    point_size = 0.5,
    point_alpha = 0.1,
    title = NULL
) {
    # 1. Extract coordinates
    coords <- .extract_coords(object, reduction)
    if (is.null(coords)) stop(sprintf("Reduction '%s' not found.", reduction))
    
    # 2. Extract weights
    weights <- .extract_weights(object, weight_metric, da_results)
    
    # Strategy:
    # If weight is NULL -> Standard point density
    # If weight is supplied -> Weighted density.
    # Weights must be non-negative. If negative (e.g. logFC), we take absolute value
    # and warn, or user should pre-process. Here we take abs() if negative detected.
    
    is_weighted <- !is.null(weights)
    if (is_weighted) {
        if (any(weights < 0, na.rm = TRUE)) {
            warning("Negative weights detected. Taking absolute values for density calculation.")
            weights <- abs(weights)
        }
        # Remove NAs
        keep <- !is.na(weights) & !is.na(coords[, 1]) & !is.na(coords[, 2])
        coords <- coords[keep, , drop = FALSE]
        weights <- weights[keep]
    } else {
        keep <- !is.na(coords[, 1]) & !is.na(coords[, 2])
        coords <- coords[keep, , drop = FALSE]
        weights <- rep(1, nrow(coords))
    }
    
    # 3. Compute 2D Density
    dens_data <- .compute_weighted_density(
        x = coords[, 1],
        y = coords[, 2],
        w = weights,
        n = grid_size,
        h = bandwidth
    )
    
    # 4. Plot
    df <- data.frame(
        x = rep(dens_data$x, each = grid_size),
        y = rep(dens_data$y, times = grid_size),
        density = as.vector(dens_data$z)
    )
    
    if (is.null(palette)) {
        palette <- viridisLite::viridis(100)
    }
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_raster(ggplot2::aes(fill = .data$density)) +
        ggplot2::scale_fill_gradientn(colors = palette, name = "Density") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "UMAP 1", y = "UMAP 2", title = title)
    
    if (show_points) {
        pt_df <- data.frame(x = coords[, 1], y = coords[, 2])
        p <- p + ggplot2::geom_point(
            data = pt_df,
            ggplot2::aes(x = .data$x, y = .data$y),
            size = point_size,
            alpha = point_alpha,
            inherit.aes = FALSE
        )
    }
    
    p
}

.extract_coords <- function(object, reduction) {
    if (inherits(object, "Milo") || inherits(object, "SingleCellExperiment")) {
        if (reduction %in% SingleCellExperiment::reducedDimNames(object)) {
            return(SingleCellExperiment::reducedDim(object, reduction))
        }
    } else if (inherits(object, "Seurat")) {
        if (reduction %in% Seurat::Reductions(object)) {
            return(Seurat::Embeddings(object, reduction = reduction))
        }
    }
    NULL
}

.extract_weights <- function(object, metric, da_results) {
    if (is.null(metric)) return(NULL)
    
    if (is.numeric(metric)) return(metric)
    
    if (is.character(metric)) {
        # Check da_results first (Milo case)
        if (!is.null(da_results) && metric %in% names(da_results)) {
            return(da_results[[metric]])
        }
        # Check meta.data (Seurat case)
        if (inherits(object, "Seurat")) {
            if (metric %in% names(object@meta.data)) {
                return(object@meta.data[[metric]])
            }
        }
    }
    stop(sprintf("Metric '%s' not found in provided data.", metric))
}

.compute_weighted_density <- function(x, y, w, n = 100, h = NULL) {
    # Calculate limits with buffer
    lims <- c(range(x), range(y))
    lims <- lims + c(-1, 1, -1, 1) * 0.05 * diff(lims)[c(1, 1, 3, 3)]
    
    # Calculate bandwidth if not provided
    if (is.null(h)) {
        h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
        h[h == 0] <- 1
    }
    
    # Create grid
    gx <- seq(lims[1], lims[2], length.out = n)
    gy <- seq(lims[3], lims[4], length.out = n)
    
    # Bin points
    ax <- cut(x, breaks = gx, labels = FALSE, include.lowest = TRUE)
    ay <- cut(y, breaks = gy, labels = FALSE, include.lowest = TRUE)
    
    # Grid matrix
    z <- matrix(0, n, n)
    
    # Aggregate weights
    valid <- !is.na(ax) & !is.na(ay)
    if (any(valid)) {
        idx <- (ay[valid] - 1) * n + ax[valid]
        sums <- tapply(w[valid], idx, sum)
        z[as.integer(names(sums))] <- as.vector(sums)
    }
    
    # Gaussian smoothing
    hx_grid <- h[1] / (gx[2] - gx[1])
    hy_grid <- h[2] / (gy[2] - gy[1])
    
    rx <- ceiling(3 * hx_grid)
    ry <- ceiling(3 * hy_grid)
    
    kx <- dnorm(seq(-rx, rx), sd = hx_grid)
    ky <- dnorm(seq(-ry, ry), sd = hy_grid)
    kx <- kx / sum(kx)
    ky <- ky / sum(ky)
    
    z_smooth <- t(apply(z, 1, function(row) stats::filter(row, kx, sides = 2)))
    z_smooth[is.na(z_smooth)] <- 0
    z_smooth <- apply(z_smooth, 2, function(col) stats::filter(col, ky, sides = 2))
    z_smooth[is.na(z_smooth)] <- 0
    
    list(x = gx, y = gy, z = z_smooth)
}
