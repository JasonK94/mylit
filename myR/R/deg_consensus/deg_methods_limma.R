# ============================================================================
# Limma-based DEG Methods (Phase 2)
# ============================================================================
# Standalone functions using create_pseudobulk_v1
# ============================================================================

#' Run limma-voom Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using limma with voom transformation.
#' Uses native create_pseudobulk_v1.
#'
#' @inheritParams runLIMMA_voom_v1
#' @export
runLIMMA_voom_v1 <- function(
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
  remove_na_groups = TRUE
) {
  params <- as.list(match.call())[-1]
  if (is.null(contrast)) stop("contrast required")

  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "limma", "edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Packages missing: ", paste(miss, collapse = ", "))

  if (!exists("create_pseudobulk_v1")) stop("create_pseudobulk_v1 not found")

  # 0. NA
  meta <- sobj@meta.data
  if (remove_na_groups) {
    cols_to_check <- c(group_id, cluster_id, sample_id, batch_id, covar_effects)
    na_mask <- rep(FALSE, nrow(meta))
    for (c in cols_to_check) if (!is.null(c) && c %in% colnames(meta)) na_mask <- na_mask | is.na(meta[[c]]) | meta[[c]] %in% c("NA", "na")
    if (sum(na_mask) > 0) {
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    }
  }
  if (length(unique(meta[[group_id]])) < 2) stop("Group < 2")

  # 1. Pseudobulk
  message("Pseudobulking...")
  pb <- create_pseudobulk_v1(sobj, cluster_id, sample_id, group_id)
  pb_meta <- SummarizedExperiment::colData(pb)

  # 3. Contrast Func
  fix_contrast <- function(contrast_str, design_cols) {
    z <- gsub("\\s+", "", contrast_str)
    toks <- unlist(strsplit(z, "([+\\-])", perl = TRUE))
    ops <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl = TRUE)))
    rebuild <- function(tok) {
      tok <- gsub("^group(_id)?", "group", tok)
      if (!grepl("^group", tok)) tok <- paste0("group", tok)
      tok
    }
    toks2 <- vapply(toks, rebuild, character(1))
    out <- toks2[1]
    if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i + 1])
    out
  }

  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb))
  if (!is.null(keep_clusters)) cluster_names <- intersect(cluster_names, keep_clusters)

  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb)[[clust]]
    s_keep <- colSums(pb_clust) > 0
    if (sum(s_keep) < 2) next
    pb_clust <- pb_clust[, s_keep, drop = FALSE]
    pb_meta_sub <- pb_meta[colnames(pb_clust), , drop = FALSE]

    # Handle Design
    if (is.null(pb_meta_sub[[group_id]])) next
    grp_col <- droplevels(factor(pb_meta_sub[[group_id]]))
    if (nlevels(grp_col) < 2) next

    design_meta <- as.data.frame(pb_meta_sub)
    design_meta$group <- grp_col

    target_terms <- c("group")
    if (!is.null(batch_id) && batch_id %in% names(design_meta)) {
      target_terms <- c(target_terms, batch_id)
    }
    if (!is.null(covar_effects)) {
      valid_cv <- covar_effects[covar_effects %in% names(design_meta)]
      target_terms <- c(target_terms, valid_cv)
    }

    hd <- handle_design_v1(design_meta, target_terms)

    # Reconstruct formula for Limma (0 + group)
    fm_str <- paste("~ 0 +", paste(hd$final_terms, collapse = " + "))
    design <- stats::model.matrix(stats::as.formula(fm_str), data = hd$meta)

    if (qr(design)$rank < ncol(design)) next
    contrast_fixed <- fix_contrast(contrast, colnames(design))
    tryCatch(
      {
        cm <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)
      },
      error = function(e) NULL
    ) -> cm
    if (is.null(cm)) next

    dge <- edgeR::DGEList(counts = pb_clust)
    keep <- edgeR::filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    if (nrow(dge) == 0) next

    dge <- edgeR::calcNormFactors(dge)
    v <- limma::voom(dge, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    fit <- limma::contrasts.fit(fit, cm)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = Inf, sort.by = "none")
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL

    if ("P.Value" %in% names(res)) res$pvalue <- res$P.Value
    if ("adj.P.Val" %in% names(res)) res$FDR <- res$adj.P.Val
    if ("t" %in% names(res)) res$statistic <- res$t

    all_results[[clust]] <- res
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
  attr(combined, "formula") <- if (exists("fm_str")) fm_str else "Unknown"
  return(combined)
}


#' Run limma-trend Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using limma with trend.
#' Uses native create_pseudobulk_v1.
#'
#' @inheritParams runLIMMA_voom_v1
#' @export
runLIMMA_trend_v1 <- function(
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
  remove_na_groups = TRUE
) {
  params <- as.list(match.call())[-1]
  if (is.null(contrast)) stop("contrast required")

  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "limma", "edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Packages missing: ", paste(miss, collapse = ", "))

  if (!exists("create_pseudobulk_v1")) stop("create_pseudobulk_v1 not found")

  # 0. NA
  meta <- sobj@meta.data
  if (remove_na_groups) {
    cols_to_check <- c(group_id, cluster_id, sample_id, batch_id, covar_effects)
    na_mask <- rep(FALSE, nrow(meta))
    for (c in cols_to_check) if (!is.null(c) && c %in% colnames(meta)) na_mask <- na_mask | is.na(meta[[c]]) | meta[[c]] %in% c("NA", "na")
    if (sum(na_mask) > 0) {
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    }
  }
  if (length(unique(meta[[group_id]])) < 2) stop("Group < 2")

  # 1. Pseudobulk
  pb <- create_pseudobulk_v1(sobj, cluster_id, sample_id, group_id)
  pb_meta <- SummarizedExperiment::colData(pb)

  # 3. Contrast Func
  fix_contrast <- function(contrast_str, design_cols) {
    z <- gsub("\\s+", "", contrast_str)
    toks <- unlist(strsplit(z, "([+\\-])", perl = TRUE))
    ops <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl = TRUE)))
    rebuild <- function(tok) {
      tok <- gsub("^group(_id)?", "group", tok)
      if (!grepl("^group", tok)) tok <- paste0("group", tok)
      tok
    }
    toks2 <- vapply(toks, rebuild, character(1))
    out <- toks2[1]
    if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i + 1])
    out
  }

  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb))
  if (!is.null(keep_clusters)) cluster_names <- intersect(cluster_names, keep_clusters)

  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb)[[clust]]
    s_keep <- colSums(pb_clust) > 0
    if (sum(s_keep) < 2) next
    pb_clust <- pb_clust[, s_keep, drop = FALSE]
    pb_meta_sub <- pb_meta[colnames(pb_clust), , drop = FALSE]

    # Handle Design
    if (is.null(pb_meta_sub[[group_id]])) next
    grp_col <- droplevels(factor(pb_meta_sub[[group_id]]))
    if (nlevels(grp_col) < 2) next

    design_meta <- as.data.frame(pb_meta_sub)
    design_meta$group <- grp_col

    target_terms <- c("group")
    if (!is.null(batch_id) && batch_id %in% names(design_meta)) {
      target_terms <- c(target_terms, batch_id)
    }
    if (!is.null(covar_effects)) {
      valid_cv <- covar_effects[covar_effects %in% names(design_meta)]
      target_terms <- c(target_terms, valid_cv)
    }

    hd <- handle_design_v1(design_meta, target_terms)

    # Reconstruct formula for Limma (0 + group)
    fm_str <- paste("~ 0 +", paste(hd$final_terms, collapse = " + "))
    design <- stats::model.matrix(stats::as.formula(fm_str), data = hd$meta)

    if (qr(design)$rank < ncol(design)) next
    contrast_fixed <- fix_contrast(contrast, colnames(design))
    tryCatch(
      {
        cm <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)
      },
      error = function(e) NULL
    ) -> cm
    if (is.null(cm)) next

    dge <- edgeR::DGEList(counts = pb_clust)
    keep <- edgeR::filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    if (nrow(dge) == 0) next

    dge <- edgeR::calcNormFactors(dge)
    logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    fit <- limma::lmFit(logCPM, design)
    fit <- limma::contrasts.fit(fit, cm)
    fit <- limma::eBayes(fit, trend = TRUE)

    res <- limma::topTable(fit, number = Inf, sort.by = "none")
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL

    if ("P.Value" %in% names(res)) res$pvalue <- res$P.Value
    if ("adj.P.Val" %in% names(res)) res$FDR <- res$adj.P.Val
    if ("t" %in% names(res)) res$statistic <- res$t

    all_results[[clust]] <- res
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
  attr(combined, "formula") <- if (exists("fm_str")) fm_str else "Unknown"
  return(combined)
}
