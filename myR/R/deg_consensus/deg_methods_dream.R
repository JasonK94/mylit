# ============================================================================
# Dream-based DEG Methods (Phase 2)
# ============================================================================
# variancePartition 패키지의 voomWithDreamWeights 및 dream 활용
# ============================================================================

#' Run Dream Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using variancePartition::dream (LMM).
#' Uses native create_pseudobulk_v1.
#' Handles batch effects as Random Effects if specified.
#'
#' @inheritParams runLIMMA_voom_v1
#' @param use_batch_as_random Logical. If TRUE and batch_id is provided,
#'        batch is modeled as a random effect (1|batch). Default is TRUE.
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runDREAM_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  batch_id = NULL,
  covar_effects = NULL,
  contrast = NULL,
  pb_min_cells = 3,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE,
  use_batch_as_random = TRUE
) {
  params <- as.list(match.call())[-1]
  if (is.null(contrast)) stop("contrast required")

  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "limma", "dplyr", "variancePartition", "edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Packages missing: ", paste(miss, collapse = ", "))

  if (!exists("create_pseudobulk_v1")) stop("create_pseudobulk_v1 not found")

  # 0. Pre-check Metadata
  meta <- sobj@meta.data
  if (remove_na_groups) {
    cols_to_check <- c(group_id, cluster_id, sample_id, batch_id, covar_effects)
    na_mask <- rep(FALSE, nrow(meta))
    for (c in cols_to_check) if (!is.null(c) && c %in% colnames(meta)) na_mask <- na_mask | is.na(meta[[c]]) | meta[[c]] %in% c("NA", "na")
    if (sum(na_mask) > 0) sobj <- sobj[, !na_mask]
    meta <- sobj@meta.data
  }
  if (length(unique(meta[[group_id]])) < 2) stop("Group < 2")

  # 1. Pseudobulk
  message("Pseudobulking (Dream)...")
  pb <- create_pseudobulk_v1(sobj, cluster_id, sample_id, group_id)
  pb_meta <- SummarizedExperiment::colData(pb)

  extract_groups <- function(contrast_str, levels_available) {
    z <- gsub("\\s+", "", contrast_str)
    toks <- unique(gsub("^group(_id)?", "", unlist(strsplit(z, "[^A-Za-z0-9_]+"))))
    toks <- toks[nchar(toks) > 0]
    keep <- intersect(toks, levels_available)
    if (length(keep) < 1) {
      g2 <- levels_available[vapply(levels_available, function(g) grepl(g, z), logical(1))]
      keep <- unique(g2)
    }
    keep
  }

  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb))
  if (!is.null(keep_clusters)) cluster_names <- intersect(cluster_names, keep_clusters)

  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb)[[clust]] # Gene x Sample count matrix

    # Filter low expression
    # Dream usually expects DGEList inputs
    # Minimal Filtering
    y <- edgeR::DGEList(counts = pb_clust)
    keep <- edgeR::filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)

    if (ncol(y) < 3) {
      message(sprintf("  Cluster %s: too few samples (%d), skip.", clust, ncol(y)))
      next
    }

    # Metadata
    sample_ids <- colnames(y)
    pb_clust_meta <- pb_meta[sample_ids, , drop = FALSE]

    # Handle Design
    if (is.null(pb_clust_meta[[group_id]])) next
    grp_col <- droplevels(factor(pb_clust_meta[[group_id]]))
    if (nlevels(grp_col) < 2) next

    design_meta <- as.data.frame(pb_clust_meta)
    design_meta$group <- grp_col

    # 1. Fixed Effects
    fixed_terms <- c("group")
    random_terms <- c()

    if (!is.null(batch_id) && batch_id %in% names(design_meta)) {
      if (use_batch_as_random) {
        b_col <- droplevels(factor(design_meta[[batch_id]]))
        if (nlevels(b_col) > 1) {
          design_meta$batch <- b_col
          random_terms <- c(random_terms, "(1|batch)")
        }
      } else {
        # Fixed: assume batch_id column exists or map specific name?
        # handle_design expects exact column names for Fixed
        # Should we map batch_id -> 'batch'?
        # If we used 'batch' in fixed_terms, we must have 'batch' col in design_meta
        design_meta$batch <- droplevels(factor(design_meta[[batch_id]]))
        fixed_terms <- c(fixed_terms, "batch")
      }
    }

    if (!is.null(covar_effects)) {
      valid_cv <- covar_effects[covar_effects %in% names(design_meta)]
      fixed_terms <- c(fixed_terms, valid_cv)
    }

    # Run handle_design on Fixed
    hd <- handle_design_v1(design_meta, fixed_terms)

    # Construct Full Formula
    # Fixed part: ~ 0 + terms...
    fm_str <- paste("~ 0 +", paste(hd$final_terms, collapse = " + "))

    # Append Random
    if (length(random_terms) > 0) {
      fm_str <- paste(fm_str, "+", paste(random_terms, collapse = " + "))
    }

    design_formula <- stats::as.formula(fm_str)

    # Check rank? Dream handles collinearity reasonably but better check.
    # But with Random Effects, model.matrix is not straightforward to check rank.
    # We skip explicit rank check for LMM part, trusting variancePartition.

    # Voom with Dream Weights
    message(sprintf("  Cluster %s: voomWithDreamWeights...", clust))
    vobj <- tryCatch(
      {
        variancePartition::voomWithDreamWeights(y, design_formula, hd$meta, quiet = TRUE)
      },
      error = function(e) {
        message("    voomWithDreamWeights failed: ", e$message)
        return(NULL)
      }
    )
    if (is.null(vobj)) next

    # Fit Dream
    message(sprintf("  Cluster %s: dream fitting...", clust))
    fit <- tryCatch(
      {
        variancePartition::dream(vobj, design_formula, hd$meta, quiet = TRUE)
      },
      error = function(e) {
        message("    dream fit failed: ", e$message)
        return(NULL)
      }
    )
    if (is.null(fit)) next

    # Contrast
    # We used ~ 0 + group ..., so coefficients are groupLEVELS
    design_cols <- colnames(fit$design) # Fixed effects part

    # Fix contrast
    # e.g. "2 - 1". group levels are "1", "2". coefs are "group1", "group2".
    # We need to map "1" -> "group1".

    # Identify group coefficients
    grp_levels <- levels(grp)
    # The coefficients will be named 'group1', 'group2', etc. if factor.

    # Need to parse contrast string
    # Helper to map user contrast "2 - 1" to "group2 - group1"
    fix_contrast_dream <- function(c_str, d_cols) {
      # Reuse logic from deg_methods_edger/limma but adapted
      # If d_cols contains "groupLevel", assume standard mapping

      # Simple heuristic: replace words matching group levels with "group"+word
      # ONLY if "group"+word exists in d_cols

      final_c <- c_str
      for (l in grp_levels) {
        target <- paste0("group", l)
        if (target %in% d_cols) {
          # Replace whole word 'l' with 'target'
          # Use word boundaries
          final_c <- gsub(paste0("\\b", l, "\\b"), target, final_c)
        }
      }
      return(final_c)
    }

    final_contrast <- fix_contrast_dream(contrast, design_cols)

    # makeContrasts
    # Note: variancePartition::getContrast can be used or limma::makeContrasts
    # dream returns MArrayLM-like object, compatible with limma contrasts

    cm <- tryCatch(
      {
        limma::makeContrasts(contrasts = final_contrast, levels = design_cols)
        # Note: design_cols only includes Fixed Effects. That's correct for contrast.
      },
      error = function(e) {
        message("    makeContrasts failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(cm)) next

    fit2 <- limma::contrasts.fit(fit, cm)
    fit2 <- limma::eBayes(fit2)

    # TopTable
    res <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P")

    res_df <- as.data.frame(res)
    res_df$cluster_id <- clust
    res_df$gene <- rownames(res_df)
    rownames(res_df) <- NULL

    if ("adj.P.Val" %in% names(res_df)) res_df$FDR <- res_df$adj.P.Val
    if ("P.Value" %in% names(res_df)) res_df$pvalue <- res_df$P.Value
    # t-statistic is already 't'
    if ("t" %in% names(res_df)) res_df$statistic <- res_df$t

    all_results[[clust]] <- res_df
  }

  if (length(all_results) == 0) stop("All clusters failed")
  combined <- do.call(rbind, all_results)
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }
  attr(combined, "run_info") <- params
  attr(combined, "formula") <- "Varies per cluster (LMM)"
  return(combined)
}
