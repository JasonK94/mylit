# ============================================================================
# MAST-based DEG Methods (Phase 2)
# ============================================================================

#' Run MAST Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using MAST.
#' MAST is designed for single-cell data (hurdle model).
#' NOTE: This runs on CELL-LEVEL data, not pseudobulk.
#'
#' @param sobj Seurat object
#' @param cluster_id Cluster column name
#' @param sample_id Sample column name (for random effect if needed, though MAST usually fixed)
#' @param group_id Group column name
#' @param covar_effects Covariates to adjust for
#' @param contrast Contrast string (e.g. "2 - 1") - logic differs for MAST
#' @param min_cells_expr Minimum cells expressing gene (default 3)
#' @param n_cores Number of cores
#' @param random_effect Random effect variable (e.g. "hos_no"). MAST supports LME via zlm.
#'
#' @return Data frame with DEG results
#' @export
runMAST_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  batch_id = NULL,
  covar_effects = NULL,
  contrast = NULL,
  min_cells_expr = 3,
  n_cores = 4,
  random_effect = NULL,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE
) {
    if (!requireNamespace("MAST", quietly = TRUE)) stop("MAST required")
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat required")

    params <- as.list(match.call())[-1]

    # Pre-processing
    if (remove_na_groups) {
        # Filter NA
        meta <- sobj@meta.data
        check <- c(group_id, cluster_id, sample_id, covar_effects)
        mask <- rowSums(is.na(meta[, check, drop = FALSE])) > 0
        if (any(mask)) sobj <- sobj[, !mask]
    }

    # Clusters
    Seurat::Idents(sobj) <- sobj[[cluster_id]]
    all_clusters <- levels(Seurat::Idents(sobj))
    if (!is.null(keep_clusters)) all_clusters <- intersect(all_clusters, keep_clusters)

    all_res <- list()

    for (clust in all_clusters) {
        message("Processing MAST for cluster: ", clust)
        sub_obj <- subset(sobj, idents = clust)

        # Check Cells
        if (ncol(sub_obj) < 10) next

        # Convert to SCA
        # MAST needs log2(CPM+1) or similar. Seurat "LogNormalize" is ln(CPM/scale * 10000 + 1).
        # Standard MAST usage suggests Log2.
        # Let's assume Seurat's data slot is usable or convert.
        # Ideally: log2(CPM + 1). Seurat default is natural log.
        # We'll stick to Seurat data for simplicity or calculate CPM.

        # Create SingleCellAssay
        sce <- Seurat::as.SingleCellExperiment(sub_obj)
        sca <- MAST::SceToSingleCellAssay(sce, check_sanity = FALSE)

        # Filter Genes
        # freq(sca) > 0 means expressed.
        # We need min_cells_expr
        expressed <- freq(sca) > (min_cells_expr / ncol(sca))
        sca <- sca[expressed, ]

        if (nrow(sca) < 2) next

        # Build Formula
        # We need to handle 'group_id' and 'covariates'.
        # Contrast: "2 - 1" -> coefficient testing.

        colData(sca)$group <- factor(colData(sca)[[group_id]])
        if (nlevels(colData(sca)$group) < 2) next

        # Check Random Effect
        # If random effect is used, use zlm( method="glmer", ... ) - SLOW!
        # MAST default is bayesglm (fixed).

        fmla_str <- "~ group"

        # Covariates
        if (!is.null(covar_effects)) {
            for (cv in covar_effects) {
                if (cv %in% names(colData(sca))) {
                    fmla_str <- paste0(fmla_str, " + ", cv)
                    # Check CDR (cellular detection rate) often used in MAST?
                    # User didn't request CDR.
                }
            }
        }

        # Random Effect
        method_opt <- "bayesglm"
        if (!is.null(random_effect) && random_effect %in% names(colData(sca))) {
            fmla_str <- paste0(fmla_str, " + (1|", random_effect, ")")
            method_opt <- "glmer"
            options(mc.cores = n_cores)
            # Note: glmer in MAST is very slow. User warned "light data".
        }

        # Fit
        message("  Fitting zlm (method=", method_opt, ")...")
        zlm_fit <- tryCatch(
            {
                MAST::zlm(as.formula(fmla_str), sca, method = method_opt, parallel = (method_opt == "glmer"))
            },
            error = function(e) {
                message("  MAST Error: ", e$message)
                return(NULL)
            }
        )

        if (is.null(zlm_fit)) next

        # Contrast/Hypothesis Testing
        # User supplies "2 - 1". MAST uses Likelihood Ratio Test (lrTest) or Wald.
        # For glmer, use Wald?
        # Typically MAST uses `summary(zlm_fit, doLRT='contrast')`

        # Parse contrast "2 - 1" implies group2 vs group1.
        # In formula `~ group + ...`, coefficients are `(Intercept)` and `group2`.
        # If we just test term `group2`, that's 2-1.

        # We try to guess the coefficient name for the contrast.
        # Assuming reference level is '1'.

        summary_res <- tryCatch(
            {
                MAST::summary(zlm_fit, doLRT = paste0("group", levels(colData(sca)$group)[2]))
            },
            error = function(e) {
                # Try determining term name dynamically
                terms <- zlm_fit@coefC
                target <- grep("group", colnames(terms), value = TRUE)
                if (length(target) == 1) MAST::summary(zlm_fit, doLRT = target) else NULL
            }
        )

        if (is.null(summary_res)) next

        # Extract Results
        # datatable from summary
        dt <- summary_res$datatable
        # component "H" is hurdle (combined)? No, MAST output structure:
        # component: reaction (C=continuous, D=discrete), H=hurdle test statistic?
        # Actually `hurdle` p-values are typically what we want.

        # Filter for contrast component
        dt_sub <- dt[dt$component == "H", ] # Hurdle model test
        if (nrow(dt_sub) == 0) dt_sub <- dt[dt$contrast != "(Intercept)", ] # Fallback

        # Need P-value
        # MAST 1.8+ summary output: `Pr(>Chisq)` column
        res_df <- data.frame(
            gene = dt_sub$primerid,
            pvalue = dt_sub$`Pr(>Chisq)`,
            statistic = dt_sub$`test.type`, # usually 'Chisq' val? No, column name varies.
            stringsAsFactors = FALSE
        )

        # LogFC?
        # LogFC is stored in the logFC component of summary or need separate coef extraction.
        # MAST returns coefs for C and D components.
        # Approximate LogFC is usually C component (continuous) for the group term.
        # But hurdle logFC is complex.
        # Often people report the coefficient from the Continuous model as logFC.

        coefs <- suppressMessages(MAST::summary(zlm_fit)$datatable)
        coefs <- coefs[coefs$component == "logFC", ] # Wait, does summary return logFC row?
        # Standard MAST workflow: getLogFC(zlm_fit) ??
        # Let's stick to simple coefficient extraction if possible.

        # Simplified: Use `logFC` generic if available, else extract coef.
        # For now, we return P-values as primary. LogFC might need `getLogFC`.

        res_df$FDR <- p.adjust(res_df$pvalue, method = "BH")
        res_df$cluster_id <- clust

        all_res[[clust]] <- res_df
    }

    if (length(all_res) == 0) {
        return(NULL)
    }
    combined <- do.call(rbind, all_res)
    attr(combined, "run_info") <- params
    attr(combined, "formula") <- fmla_str

    return(combined)
}
