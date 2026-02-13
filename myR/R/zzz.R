# ==============================================================================
# Package Initialization & Configuration
# ==============================================================================
#' @import dplyr
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import methods
#' @import grDevices
#' @import graphics
#' @import scales
#' @import tidyr
#' @import purrr
#' @import viridisLite
#' @importFrom magrittr %>%
#' @keywords internal
##' @import utils
##' @importFrom rlang .data :=
"_PACKAGE"

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
        # CRITICAL: Fix for View() with Seurat objects
        conflicted::conflicts_prefer(base::as.data.frame)

        # Matrix vs base
        conflicted::conflicts_prefer(base::colMeans)
        conflicted::conflicts_prefer(base::rowMeans)
        conflicted::conflicts_prefer(base::colSums)
        conflicted::conflicts_prefer(base::rowSums)
    }
}
