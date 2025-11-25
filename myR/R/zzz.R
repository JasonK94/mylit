# ==============================================================================
# Package Initialization & Configuration
# ==============================================================================

.onLoad <- function(libname, pkgname) {
    # Set conflicted preferences automatically
    if (requireNamespace("conflicted", quietly = TRUE)) {
        # dplyr preferences
        conflicted::conflicts_prefer(dplyr::select, quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::summarise, quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::filter, quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::mutate, quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::arrange, quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::summarize, quiet = TRUE)

        # base preferences
        conflicted::conflicts_prefer(base::intersect, quiet = TRUE)
        conflicted::conflicts_prefer(base::setdiff, quiet = TRUE)
        conflicted::conflicts_prefer(base::union, quiet = TRUE)

        # Matrix vs base (Matrix package often conflicts)
        conflicted::conflicts_prefer(base::colMeans, quiet = TRUE)
        conflicted::conflicts_prefer(base::rowMeans, quiet = TRUE)
        conflicted::conflicts_prefer(base::colSums, quiet = TRUE)
        conflicted::conflicts_prefer(base::rowSums, quiet = TRUE)

        # Seurat vs others
        if (requireNamespace("Seurat", quietly = TRUE)) {
            # Add Seurat specific preferences if needed
        }
    }
}
