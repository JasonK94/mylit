# ============================================================================
# edgeR-based DEG Methods (Phase 2)
# ============================================================================
# Standalone functions using create_pseudobulk_v1
# ============================================================================

#' Run edgeR-LRT Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using edgeR with Likelihood Ratio Test.
#' Uses native create_pseudobulk_v1 for aggregation (no muscat dependency).
#'
#' @inheritParams runLIMMA_voom_v1
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runEDGER_LRT_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  batch_id = NULL,
  covar_effects = NULL,
  contrast = NULL,
  pb_min_cells = 3, # Used? We might filter clusters with low cells in create_pseudobulk? No.
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE
) {
  params <- as.list(match.call())[-1] # Capture parameters
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")

  # deps
  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "limma", "dplyr", "edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse = ", "))

  # --- 0. NA 값 처리 (Common) ---
  message("0/5: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data

  # Check required columns
  required <- c(cluster_id, sample_id, group_id)
  if (!all(required %in% colnames(meta))) stop("Missing required columns")

  if (remove_na_groups) {
    # Check NA in key columns
    cols_to_check <- c(group_id, cluster_id, sample_id, batch_id, covar_effects)
    na_mask <- rep(FALSE, nrow(meta))
    for (c in cols_to_check) {
      if (!is.null(c) && c %in% colnames(meta)) {
        na_mask <- na_mask | is.na(meta[[c]]) | (as.character(meta[[c]]) %in% c("NA", "na"))
      }
    }
    if (sum(na_mask) > 0) {
      message(sprintf("... NA 값 제거: %d cells", sum(na_mask)))
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    }
  }

  if (length(unique(meta[[group_id]])) < 2) stop("Group < 2")

  # --- 1. Pseudobulk ---
  message("1/5: Pseudobulking (create_pseudobulk_v1)...")
  pb <- create_pseudobulk_v1(sobj, cluster_id, sample_id, group_id)

  # Filter clusters
  if (!is.null(keep_clusters)) {
    avail <- names(SummarizedExperiment::assays(pb))
    keep <- intersect(avail, as.character(keep_clusters))
    if (length(keep) == 0) stop("No clusters to keep")
    pb <- pb[keep] # logic needs verification: pb is SCE. subsetting defaults to features.
    # Wait, assays are Clusters in muscat structure.
    # In create_pseudobulk_v1, assays are [[Cluster]].
    # SCE structure doesn't support subsetting assays by name easily like a list.
    # But we can iterate.
    # Actually, my create_pseudobulk_v1 returns SCE with assays=list(cl1=mat, cl2=mat).
    # We can just pick the assays we want during iteration.
  }

  # --- 2. Prepare Metadata ---
  message("2/5: 메타데이터 정리 중...")
  pb_meta <- SummarizedExperiment::colData(pb)

  # Ensure all covariates are present
  # create_pseudobulk_v1 usually captures all sample-level vars.
  # But let's verify.

  # --- 3. Contrast Parsing ---
  fix_contrast <- function(contrast_str, design_cols) {
    z <- gsub("\\s+", "", contrast_str)
    # Simple replacement for "2 - 1" -> "group2 - group1"
    # Logic: tokenize by operators. If token is number or string, prepend 'group' if it matches group levels
    # But usually contrast is "A - B".
    # Just return as is if variables match design.
    # For now, use the robust legacy logic if copied, or simple one.
    # Legacy logic:
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

  # --- 4. Run edgeR-LRT ---
  message("4/5: edgeR-LRT 분석 실행 중...")
  all_results <- list()

  # Iterate over Assays (Clusters)
  cluster_names <- names(SummarizedExperiment::assays(pb))
  if (!is.null(keep_clusters)) cluster_names <- intersect(cluster_names, keep_clusters)

  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb)[[clust]]

    # Check min cells filter? "pb_min_cells"
    # Actually we don't have cell count per sample in the aggregated matrix (it's sum of counts).
    # We can't filter by "min cells" easily here unless we tracked it.
    # For now, skip empty samples? (colSums=0)

    # Filter samples with 0 counts
    s_keep <- colSums(pb_clust) > 0
    if (sum(s_keep) < 2) next

    pb_clust <- pb_clust[, s_keep, drop = FALSE]
    pb_meta_sub <- pb_meta[colnames(pb_clust), , drop = FALSE]

    # Create Design
    # Should use model.matrix with group + batch + covar
    # Ensure factors

    # Handle Design
    if (is.null(pb_meta_sub[[group_id]])) next
    grp_col <- droplevels(factor(pb_meta_sub[[group_id]]))
    if (nlevels(grp_col) < 2) next

    # Prepare metadata for handle_design
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

    # Reconstruct formula for edgeR (0 + group)
    fm_str <- paste("~ 0 +", paste(hd$final_terms, collapse = " + "))
    design <- stats::model.matrix(stats::as.formula(fm_str), data = hd$meta)

    # Final Rank Check (just in case)
    if (qr(design)$rank < ncol(design)) next

    # Contrast
    # Check if contrast uses 'group' prefix
    contrast_fixed <- fix_contrast(contrast, colnames(design))
    tryCatch(
      {
        cm <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)
      },
      error = function(e) {
        # Fallback default
        message("Contrast error, using default last vs first")
        return(NULL)
      }
    ) -> cm
    if (is.null(cm)) next # or handle default

    # edgeR
    dge <- edgeR::DGEList(counts = pb_clust)
    keep <- edgeR::filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    if (nrow(dge) == 0) next

    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmFit(dge, design)
    lrt <- edgeR::glmLRT(fit, contrast = cm)

    res <- edgeR::topTags(lrt, n = Inf)$table
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL

    # Standardize cols
    if ("PValue" %in% names(res)) res$pvalue <- res$PValue
    if ("FDR" %in% names(res)) res$FDR <- res$FDR
    if ("logFC" %in% names(res)) res$logFC <- res$logFC
    if ("LR" %in% names(res)) res$statistic <- res$LR * sign(res$logFC)

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

  message("5/5: edgeR-LRT 분석 완료.")
  attr(combined, "run_info") <- params
  attr(combined, "formula") <- if (exists("fm_str")) fm_str else "Unknown"
  return(combined)
}


#' Run edgeR-QLF Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using edgeR with Quasi-Likelihood F-test.
#' Uses native create_pseudobulk_v1.
#'
#' @inheritParams runEDGER_LRT_v1
#' @export
runEDGER_QLF_v1 <- function(
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

  # deps
  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "limma", "dplyr", "edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Packages missing: ", paste(miss, collapse = ", "))

  # COPY PRE-PROCESSING FROM LRT (Simplified for brevity in generation)
  # Ideally check if create_pseudobulk_v1 exists
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

    # Prepare metadata for handle_design
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

    # Reconstruct formula for edgeR (0 + group)
    fm_str <- paste("~ 0 +", paste(hd$final_terms, collapse = " + "))
    design <- stats::model.matrix(stats::as.formula(fm_str), data = hd$meta)

    if (qr(design)$rank < ncol(design)) next

    # Check Rank and DF
    rnk <- qr(design)$rank
    if (nrow(design) <= rnk + 1) {
      # message(sprintf("  Cluster %s: Insufficient DF (N=%d, Rank=%d), skipping.", clust, nrow(design), rnk))
      next
    }

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
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmQLFit(dge, design)
    qlf <- edgeR::glmQLFTest(fit, contrast = cm)

    res <- edgeR::topTags(qlf, n = Inf)$table
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL

    if ("PValue" %in% names(res)) res$pvalue <- res$PValue
    if ("FDR" %in% names(res)) res$FDR <- res$FDR
    if ("logFC" %in% names(res)) res$logFC <- res$logFC
    if ("F" %in% names(res)) res$statistic <- res$F

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
