# ============================================================================
# Nebula-based DEG Methods
# ============================================================================

#' Run Nebula Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using Nebula (Negative Binomial Mixed Model).
#' Optimized for single-cell data with random effects (e.g., subject).
#' Faster than MAST glmer.
#'
#' @param sobj Seurat object
#' @param cluster_id Cluster column
#' @param sample_id Random effect grouping (e.g. subject)
#' @param group_id Fixed effect group
#'
#' @export
runNebula_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  batch_id = NULL,
  covar_effects = NULL,
  contrast = NULL,
  min_cells_expr = 3,
  n_cores = 4,
  random_effect = "hos_no",
  keep_clusters = NULL,
  remove_na_groups = TRUE
) {
    if (!requireNamespace("nebula", quietly = TRUE)) stop("nebula required")

    params <- as.list(match.call())[-1]

    # Pre-process
    meta <- sobj@meta.data

    # Nebula requires a design matrix and a counts matrix.
    # It runs per-cluster internally? Or we loop?
    # Nebula `nebula()` function handles the model. We typically provide a Count Matrix (Genes x Cells)
    # and a design dataframe (Cells x Vars).
    # It's CELL-LEVEL mixed model.

    # Preparing data
    # Filter NA
    if (remove_na_groups) {
        # Filter NA
        check <- c(group_id, cluster_id, sample_id, covar_effects)
        mask <- rowSums(is.na(meta[, check, drop = FALSE])) > 0
        if (any(mask)) sobj <- sobj[, !mask]
    }

    Seurat::Idents(sobj) <- sobj[[cluster_id]]
    clusts <- levels(Seurat::Idents(sobj))
    if (!is.null(keep_clusters)) clusts <- intersect(clusts, keep_clusters)

    all_res <- list()

    # Need to verify if Nebula supports contrast strings directly or just Model Matrix
    # Nebula usually takes: count, id (subject), pred (design matrix).

    for (cl in clusts) {
        message("Processing Nebula for cluster: ", cl)
        sub <- subset(sobj, idents = cl)
        if (ncol(sub) < 50) {
            message("  Skip: too few cells (<50)")
            next
        }

        # Count Matrix (Standard sparse)
        # Nebula handles sparse matrices? Yes, usually.
        counts <- Seurat::GetAssayData(sub, slot = "counts")

        # Filter Genes
        # min_cells
        kp <- rowSums(counts > 0) >= min_cells_expr
        counts <- counts[kp, ]
        if (nrow(counts) < 5) next

        # Design
        df <- sub@meta.data
        df$group <- factor(df[[group_id]])
        df$sid <- factor(df[[sample_id]]) # Random effect ID

        # Model Matrix
        # Formula: ~ group + covars
        fmla_str <- "~ group"
        if (!is.null(covar_effects)) {
            for (cv in covar_effects) if (cv %in% names(df)) fmla_str <- paste0(fmla_str, " + ", cv)
        }

        design <- model.matrix(as.formula(fmla_str), data = df)

        # Run Nebula
        # nebula(count, id, pred, offset=NULL, model="NBGMM", ...)
        # id should be vector of subject IDs

        res_neb <- tryCatch(
            {
                nebula::nebula(
                    count = counts,
                    id = df$sid,
                    pred = design,
                    offset = NULL, # Or library size? Nebula calculates offset if NULL?
                    # Docs say: "offset: a vector of offsets... if NULL, log(total counts) used"
                    model = "NBGMM", # Negative Binomial Gaussian Mixture (approx for RE)
                    ncore = n_cores,
                    verbose = FALSE
                )
            },
            error = function(e) {
                message("  Nebula Error: ", e$message)
                return(NULL)
            }
        )

        if (is.null(res_neb)) next

        # Extract Results
        # res_neb$summary is a Data Frame with stats for each gene and predictor
        sum_df <- res_neb$summary

        # We need the coefficient for the Group contrast.
        # Design matrix has intercept and 'group2'. Contrast "2 - 1" is testing 'group2'.
        # Column names in summary will be "logFC_group2", "p_group2", etc.

        target_col_p <- paste0("p_group", levels(df$group)[2])
        target_col_fc <- paste0("logFC_group", levels(df$group)[2])

        if (!target_col_p %in% colnames(sum_df)) next

        df_out <- data.frame(
            gene = sum_df$gene,
            pvalue = sum_df[[target_col_p]],
            logFC = sum_df[[target_col_fc]],
            statistic = NA, # Nebula doesn't explicitly output z/t in main table?
            cluster_id = cl,
            stringsAsFactors = FALSE
        )
        df_out$FDR <- p.adjust(df_out$pvalue, method = "BH")

        all_res[[cl]] <- df_out
    }

    if (length(all_res) == 0) {
        return(NULL)
    }
    combined <- do.call(rbind, all_res)
    attr(combined, "run_info") <- params
    attr(combined, "formula") <- paste0(fmla_str, " + (1|", sample_id, ")")

    return(combined)
}
