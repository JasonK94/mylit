#' Kernel Density / Spatial Average Plot for Milo Neighborhoods
#'
#' Computes and plots a 2D kernel density estimate OR spatial average on a dimension reduction layout.
#'
#' @param object A `Milo` object, `Seurat` object, or `SingleCellExperiment`.
#' @param reduction Character. Name of the reduction slot to use (default: "UMAP").
#' @param weight_metric Numeric vector or character string.
#' @param da_results Data frame.
#' @param bandwidth Numeric vector of length 2.
#' @param grid_size Integer. Number of grid points in each direction (default: 100).
#' @param palette Character vector.
#' @param alpha_range Numeric vector of length 2.
#' @param show_points Logical.
#' @param point_size Numeric.
#' @param point_alpha Numeric.
#' @param title Character.
#' @param mode Character. "density" (default) or "average". 
#'   - "density": Computes weighted density (sum of weights). Weights must be non-negative (abs used if negative).
#'   - "average": Computes spatially smoothed average of weights (sum(w*k) / sum(k)). Preserves sign.
#'
#' @return A ggplot object.
#'
#' @importFrom MASS kde2d bandwidth.nrd
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn geom_point theme_minimal labs scale_alpha_continuous scale_fill_gradient2
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
    title = NULL,
    mode = c("density", "average")
) {
    mode <- match.arg(mode)
    
    # 1. Extract coordinates
    coords <- .extract_coords(object, reduction)
    if (is.null(coords)) stop(sprintf("Reduction '%s' not found.", reduction))
    
    # 2. Extract weights
    weights <- .extract_weights(object, weight_metric, da_results)
    
    # If object is Milo and weights match nhoods, compute nhood centers
    if (inherits(object, "Milo") && !is.null(weights)) {
        n_nhoods <- tryCatch(ncol(miloR::nhoods(object)), error=function(e) 0)
        n_cells <- nrow(coords)
        
        if (length(weights) == n_nhoods && n_nhoods != n_cells) {
            nhood_centers <- .compute_nhood_centers(object, reduction)
            if (!is.null(nhood_centers)) {
                coords <- nhood_centers
            } else {
                warning("Could not compute neighborhood centers. Using cell coordinates.")
            }
        }
    }
    
    # Validation
    if (length(weights) != nrow(coords)) {
        if (length(weights) == nrow(.extract_coords(object, reduction))) {
            coords <- .extract_coords(object, reduction) # Fallback to cell coords
        } else {
            if (!is.null(weights)) {
                stop(sprintf("Length of weights (%d) does not match number of points (%d).", 
                             length(weights), nrow(coords)))
            }
        }
    }
    
    # Handle weights based on mode
    is_weighted <- !is.null(weights)
    
    if (mode == "density") {
        if (is_weighted) {
            if (any(weights < 0, na.rm = TRUE)) {
                warning("Negative weights detected. Taking absolute values for density calculation.")
                weights <- abs(weights)
            }
            keep <- !is.na(weights) & !is.na(coords[, 1]) & !is.na(coords[, 2])
            coords <- coords[keep, , drop = FALSE]
            weights <- weights[keep]
        } else {
            keep <- !is.na(coords[, 1]) & !is.na(coords[, 2])
            coords <- coords[keep, , drop = FALSE]
            weights <- rep(1, nrow(coords))
        }
    } else if (mode == "average") {
        if (!is_weighted) stop("Mode 'average' requires a weight metric.")
        keep <- !is.na(weights) & !is.na(coords[, 1]) & !is.na(coords[, 2])
        coords <- coords[keep, , drop = FALSE]
        weights <- weights[keep]
    }
    
    # 3. Compute 2D Surface
    if (mode == "density") {
        surface_data <- .compute_weighted_density(
            x = coords[, 1],
            y = coords[, 2],
            w = weights,
            n = grid_size,
            h = bandwidth
        )
        z_vals <- as.vector(surface_data$z)
        fill_label <- "Density"
    } else {
        # Average mode
        w_sum <- .compute_weighted_density(
            x = coords[, 1],
            y = coords[, 2],
            w = weights,
            n = grid_size,
            h = bandwidth
        )
        d_sum <- .compute_weighted_density(
            x = coords[, 1],
            y = coords[, 2],
            w = rep(1, length(weights)),
            n = grid_size,
            h = bandwidth
        )
        
        den <- as.vector(d_sum$z)
        num <- as.vector(w_sum$z)
        
        # Mask low density regions
        # Threshold: e.g., 1% of max density is too aggressive if peaks are sharp.
        # Use a lower threshold or quantile.
        # Let's use 0.1% (0.001) or even lower to avoid excessive masking.
        mask_threshold <- max(den) * 0.001
        mask <- den < mask_threshold
        z_vals <- num / den
        z_vals[mask] <- NA
        
        surface_data <- w_sum
        fill_label <- "Average"
    }
    
    # 4. Plot
    df <- data.frame(
        x = rep(surface_data$x, each = grid_size),
        y = rep(surface_data$y, times = grid_size),
        value = z_vals
    )
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_raster(ggplot2::aes(fill = .data$value)) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "UMAP 1", y = "UMAP 2", title = title)
    
    if (mode == "density") {
        if (is.null(palette)) palette <- viridisLite::viridis(100)
        p <- p + ggplot2::scale_fill_gradientn(colors = palette, name = fill_label)
    } else {
        p <- p + ggplot2::scale_fill_gradient2(
            low = "blue", mid = "white", high = "red", 
            midpoint = 0, 
            name = fill_label,
            na.value = "transparent"
        )
    }
    
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

.compute_nhood_centers <- function(milo, reduction) {
    if (!reduction %in% SingleCellExperiment::reducedDimNames(milo)) return(NULL)
    
    cell_coords <- SingleCellExperiment::reducedDim(milo, reduction)
    nhood_mat <- miloR::nhoods(milo)
    
    if (is.null(nhood_mat)) return(NULL)
    
    nhood_sums <- Matrix::t(nhood_mat) %*% as.matrix(cell_coords[, 1:2])
    nhood_counts <- Matrix::colSums(nhood_mat)
    nhood_counts[nhood_counts == 0] <- 1
    
    centers <- nhood_sums / nhood_counts
    return(as.matrix(centers))
}

.extract_weights <- function(object, metric, da_results) {
    if (is.null(metric)) return(NULL)
    if (is.numeric(metric)) return(metric)
    
    if (is.character(metric)) {
        if (!is.null(da_results) && metric %in% names(da_results)) {
            if (inherits(object, "Milo") || inherits(object, "SingleCellExperiment")) {
                n_nhoods <- ncol(miloR::nhoods(object))
                if ("Nhood" %in% names(da_results)) {
                    w <- rep(NA, n_nhoods)
                    idx <- as.integer(da_results$Nhood)
                    valid_idx <- idx >= 1 & idx <= n_nhoods
                    w[idx[valid_idx]] <- da_results[[metric]][valid_idx]
                    return(w)
                } else {
                    if (nrow(da_results) == n_nhoods) return(da_results[[metric]])
                    return(da_results[[metric]])
                }
            }
            return(da_results[[metric]])
        }
        if (inherits(object, "Seurat")) {
            if (metric %in% names(object@meta.data)) return(object@meta.data[[metric]])
        }
    }
    stop(sprintf("Metric '%s' not found.", metric))
}

.compute_weighted_density <- function(x, y, w, n = 100, h = NULL) {
    lims <- c(range(x), range(y))
    lims <- lims + c(-1, 1, -1, 1) * 0.05 * diff(lims)[c(1, 1, 3, 3)]
    
    if (is.null(h)) {
        h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
        h[h == 0] <- 1
    }
    
    gx <- seq(lims[1], lims[2], length.out = n)
    gy <- seq(lims[3], lims[4], length.out = n)
    
    ax <- cut(x, breaks = gx, labels = FALSE, include.lowest = TRUE)
    ay <- cut(y, breaks = gy, labels = FALSE, include.lowest = TRUE)
    
    z <- matrix(0, n, n)
    valid <- !is.na(ax) & !is.na(ay)
    if (any(valid)) {
        idx <- (ay[valid] - 1) * n + ax[valid]
        sums <- tapply(w[valid], idx, sum)
        z[as.integer(names(sums))] <- as.vector(sums)
    }
    
    hx_grid <- h[1] / (gx[2] - gx[1])
    hy_grid <- h[2] / (gy[2] - gy[1])
    
    rx <- min(ceiling(3 * hx_grid), floor((n - 1) / 2))
    ry <- min(ceiling(3 * hy_grid), floor((n - 1) / 2))
    
    kx <- dnorm(seq(-rx, rx), sd = hx_grid)
    ky <- dnorm(seq(-ry, ry), sd = hy_grid)
    kx <- kx / sum(kx)
    ky <- ky / sum(ky)
    
    z_smooth <- t(apply(z, 1, function(row) {
        if (length(kx) >= length(row)) return(row)
        res <- stats::filter(row, kx, sides = 2)
        res[is.na(res)] <- 0
        res
    }))
    
    z_smooth <- apply(z_smooth, 2, function(col) {
        if (length(ky) >= length(col)) return(col)
        res <- stats::filter(col, ky, sides = 2)
        res[is.na(res)] <- 0
        res
    })
    
    list(x = gx, y = gy, z = z_smooth)
}
