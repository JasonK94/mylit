# ==============================================================================
# Test Utilities for FGS/TML
# ==============================================================================

#' Smart Subset for Testing
#'
#' Subsets a Seurat object while preserving rare groups.
#' Instead of random sampling across all cells, it limits the maximum number
#' of cells per group (e.g., cluster). Rare groups are kept intact.
#'
#' @param sobj Seurat object
#' @param group_var Column name for grouping (e.g., cluster variable)
#' @param max_cells_per_group Maximum number of cells to keep per group
#' @param seed Random seed
#'
#' @return Subsets Seurat object
smart_subset <- function(sobj, group_var, max_cells_per_group = 200, seed = 42) {
    if (!group_var %in% colnames(sobj@meta.data)) {
        warning(sprintf("Group var '%s' not found. Using random subset.", group_var))
        set.seed(seed)
        if (ncol(sobj) > max_cells_per_group * 5) {
            cells <- sample(colnames(sobj), max_cells_per_group * 5)
            return(subset(sobj, cells = cells))
        }
        return(sobj)
    }

    # Ensure group var is factor
    groups <- as.factor(sobj@meta.data[[group_var]])
    unique_groups <- levels(groups)

    cells_to_keep <- character(0)

    set.seed(seed)

    cat(sprintf("Smart subsetting by '%s' (max %d/group):\n", group_var, max_cells_per_group))

    for (g in unique_groups) {
        g_cells <- rownames(sobj@meta.data)[groups == g]
        n_g <- length(g_cells)

        if (n_g <= max_cells_per_group) {
            # Keep all cells in rare groups
            keep <- g_cells
            cat(sprintf("  - %s: %d (All kept)\n", g, n_g))
        } else {
            # Downsample abundant groups
            keep <- sample(g_cells, max_cells_per_group)
            cat(sprintf("  - %s: %d -> %d\n", g, n_g, max_cells_per_group))
        }
        cells_to_keep <- c(cells_to_keep, keep)
    }

    cat(sprintf("Total cells: %d -> %d\n", ncol(sobj), length(cells_to_keep)))

    subset(sobj, cells = cells_to_keep)
}

#' Prepare Data for FGS Test
#'
#' Loads and preprocesses data based on config
prepare_test_data <- function(config, downsample = TRUE, max_cells_per_group = 200) {
    if (!file.exists(config$path)) {
        stop("Data file not found: ", config$path)
    }

    cat("Loading data from:", config$path, "\n")
    sobj <- qs::qread(config$path)

    # Filter Target NAs (specific to stroke dataset logic)
    if (!is.null(config$target_na_value)) {
        target_vals <- as.character(sobj@meta.data[[config$target_var]])
        valid_cells <- target_vals != config$target_na_value & !is.na(target_vals)

        n_removed <- sum(!valid_cells)
        if (n_removed > 0) {
            cat(sprintf("Removing %d cells with target '%s' or NA\n", n_removed, config$target_na_value))
            sobj <- subset(sobj, cells = rownames(sobj@meta.data)[valid_cells])
        }
    }

    # Downsample if requested
    if (downsample) {
        sobj <- smart_subset(sobj, config$cluster_var, max_cells_per_group)
    }

    return(sobj)
}
