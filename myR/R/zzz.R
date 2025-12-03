# ==============================================================================
# Package Initialization & Configuration
# ==============================================================================

.onLoad <- function(libname, pkgname) {
    # Set conflicted preferences automatically
    if (requireNamespace("conflicted", quietly = TRUE)) {
        # dplyr preferences
        conflicted::conflicts_prefer(dplyr::select)
        conflicted::conflicts_prefer(dplyr::summarise)
        conflicted::conflicts_prefer(dplyr::filter)
        conflicted::conflicts_prefer(dplyr::mutate)
        conflicted::conflicts_prefer(dplyr::arrange)
        conflicted::conflicts_prefer(dplyr::summarize)

        # base preferences
        conflicted::conflicts_prefer(base::intersect)
        conflicted::conflicts_prefer(base::setdiff)
        conflicted::conflicts_prefer(base::union)

        # Matrix vs base (Matrix package often conflicts)
        conflicted::conflicts_prefer(base::colMeans)
        conflicted::conflicts_prefer(base::rowMeans)
        conflicted::conflicts_prefer(base::colSums)
        conflicted::conflicts_prefer(base::rowSums)

        # Seurat vs others
        if (requireNamespace("Seurat", quietly = TRUE)) {
            # Add Seurat specific preferences if needed
        }
    }
}
