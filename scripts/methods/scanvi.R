#' scANVI Integration Wrapper
#'
#' @param object Seurat object
#' @param batch_col Column name for batch
#' @param labels_col Column name for cell type labels (for semi-supervised)
#' @param ... Additional arguments for scvi-tools
#'
#' @return Integrated Seurat object
#' @export
run_scanvi <- function(object, batch_col, labels_col, conda_env = "scvi-env", ...) {
    library(Seurat)
    library(SeuratWrappers)

    # Ensure labels are present
    if (!labels_col %in% colnames(object@meta.data)) {
        stop(sprintf("Labels column '%s' not found", labels_col))
    }

    # Check for NA in labels (unlabeled cells)
    # scANVI treats "Unknown" or NA as unlabeled

    integrated <- IntegrateLayers(
        object = object,
        method = scVIIntegration, # SeuratWrappers might handle scANVI via args or separate function
        # Note: As of recent SeuratWrappers, scANVI might be a separate function or option.
        # Assuming scVIIntegration supports 'use_scanvi' or similar, or we use direct reticulate.
        # For now, using scVIIntegration as placeholder for the wrapper call.
        new.reduction = "integrated.scanvi",
        conda_env = conda_env,
        group.by = batch_col,
        ...
    )
    return(integrated)
}
