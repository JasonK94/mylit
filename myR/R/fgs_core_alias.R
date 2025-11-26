#' Find Gene signature v5.4 (bug fixes for ranger/glmnet/NMF paths)
#'
#' @description Wrapper around [find_gene_signature_v5_impl] that enables the
#' `method_impl = "v5.4"` execution path, activating the latest stability
#' fixes for ranger/glmnet/NMF.
#'
#' @inheritParams find_gene_signature_v5.3
#' @return Same as [find_gene_signature_v5.3]
#' @seealso [find_gene_signature_v5.3]
#' @export
find_gene_signature_v5.4 <- function(...) {
  find_gene_signature_v5_impl(..., method_impl = "v5.4")
}


#' @export
find_gene_signature_v5.3 <- function(data,
                                     meta.data = NULL,
                                     target_var,
                                     target_group = NULL,
                                     control_vars = NULL,
                                     method = c(
                                       "random_forest", "random_forest_ranger",
                                       "lasso", "ridge", "elastic_net",
                                       "pca_loadings", "nmf_loadings",
                                       "gam", "limma", "wilcoxon",
                                       "xgboost"
                                     ),
                                     n_features = 50,
                                     test_n = NULL,
                                     preprocess = TRUE,
                                     min_cells = 10,
                                     min_pct = 0.01,
                                     return_model = FALSE,
                                     fgs_seed = 42,
                                     lambda_selection = "lambda.1se",
                                     enet.alpha = 0.5,
                                     pca.n_pcs = 1,
                                     gam.min_unique = 15,
                                     gam.k = NULL,
                                     gam.k_dynamic_factor = 5,
                                     ...) {
  find_gene_signature_v5_impl(
    data = data,
    meta.data = meta.data,
    target_var = target_var,
    target_group = target_group,
    control_vars = control_vars,
    method = method,
    n_features = n_features,
    test_n = test_n,
    preprocess = preprocess,
    min_cells = min_cells,
    min_pct = min_pct,
    return_model = return_model,
    fgs_seed = fgs_seed,
    lambda_selection = lambda_selection,
    enet.alpha = enet.alpha,
    pca.n_pcs = pca.n_pcs,
    gam.min_unique = gam.min_unique,
    gam.k = gam.k,
    gam.k_dynamic_factor = gam.k_dynamic_factor,
    method_impl = "v5.3",
    ...
  )
}
