# ============================================================================
# DESeq2-based DEG Methods (Phase 2)
# ============================================================================
# runMUSCAT 스타일로 표준화된 DESeq2 계열 방법론
# ============================================================================

#' Run DESeq2-Wald Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using DESeq2 with Wald test.
#' This function follows the same interface as runMUSCAT for consistency.
#'
#' @inheritParams runLIMMA_voom_v1
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
#' Run DESeq2-Wald Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using DESeq2 with Wald test.
#' Uses native create_pseudobulk_v1 for aggregation (no muscat dependency).
#'
#' @inheritParams runLIMMA_voom_v1
#' @export
runDESEQ2_Wald_v1 <- function(
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

  c_str <- contrast # preserve

  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "limma", "dplyr", "DESeq2")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Packages missing: ", paste(miss, collapse = ", "))

  if (!exists("create_pseudobulk_v1")) stop("create_pseudobulk_v1 not found")

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
  message("Pseudobulking (DESeq2)...")
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
    pb_clust <- SummarizedExperiment::assays(pb)[[clust]]
    s_keep <- colSums(pb_clust) > 0
    if (sum(s_keep) < 2) next
    pb_clust <- pb_clust[, s_keep, drop = FALSE]

    sample_ids <- colnames(pb_clust)
    pb_clust_meta <- pb_meta[sample_ids, , drop = FALSE]


    # Handle Design
    if (is.null(pb_clust_meta[[group_id]])) next
    grp_col <- droplevels(factor(pb_clust_meta[[group_id]]))
    if (nlevels(grp_col) < 2) next

    design_meta <- as.data.frame(pb_clust_meta)
    design_meta$group <- grp_col

    target_terms <- c("group")
    if (!is.null(batch_id) && batch_id %in% names(design_meta)) {
      target_terms <- c("batch", target_terms) # batch first if present? DESeq2 doesn't care strictly but good style
      if ("batch" %in% names(design_meta)) {
        # Rename provided batch column to 'batch' or map it?
        # existing code: df_design$batch <- b
        design_meta$batch <- droplevels(factor(design_meta[[batch_id]]))
      }
    }
    # Wait, existing code mapped batch_id -> 'batch' column in df_design.
    # And used 'batch' in formula.
    # To use handle_design_v1 properly with my previous logic, I should setup design_meta matching terms.

    # Let's align with handle_design_v1 expectation:
    # 1. Create design_meta with FINAL column names I want in formula.

    dm <- data.frame(row.names = sample_ids)
    dm$group <- grp_col

    t_terms <- c("group")

    if (!is.null(batch_id) && batch_id %in% names(pb_clust_meta)) {
      dm$batch <- droplevels(factor(pb_clust_meta[[batch_id]]))
      t_terms <- c("batch", t_terms)
    }

    if (!is.null(covar_effects)) {
      for (cv in covar_effects) {
        if (cv %in% names(pb_clust_meta)) {
          # Copy as is, handle_design will factor/scale
          dm[[cv]] <- pb_clust_meta[[cv]]
          t_terms <- c(t_terms, cv)
        }
      }
    }

    hd <- handle_design_v1(dm, t_terms)
    design_formula <- hd$formula

    # Round for DESeq2
    pb_clust_int <- round(as.matrix(pb_clust))

    dds <- tryCatch(
      {
        DESeq2::DESeqDataSetFromMatrix(countData = pb_clust_int, colData = hd$meta, design = design_formula)
      },
      error = function(e) NULL
    )
    if (is.null(dds)) next

    dds <- tryCatch(
      {
        DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
      },
      error = function(e) NULL
    )
    if (is.null(dds)) next


    # Save used formula string for this cluster (optional, but good for debug)
    # But attributes are for aggregated result.

    # Results
    grp_levels <- levels(hd$meta$group)
    contrast_groups <- extract_groups(c_str, grp_levels)
    if (length(contrast_groups) >= 2) {
      res <- tryCatch(
        {
          DESeq2::results(dds, contrast = c("group", contrast_groups[1], contrast_groups[2]))
        },
        error = function(e) NULL
      )
    } else {
      # Default last vs first
      res <- tryCatch(
        {
          DESeq2::results(dds, contrast = c("group", tail(grp_levels, 1), grp_levels[1]))
        },
        error = function(e) NULL
      )
    }
    if (is.null(res)) next

    res_df <- as.data.frame(res)
    res_df$cluster_id <- clust
    res_df$gene <- rownames(res_df)
    rownames(res_df) <- NULL

    if ("log2FoldChange" %in% names(res_df)) res_df$logFC <- res_df$log2FoldChange
    if ("padj" %in% names(res_df)) res_df$FDR <- res_df$padj
    if ("stat" %in% names(res_df)) res_df$statistic <- res_df$stat

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
  attr(combined, "formula") <- "Varies per cluster (Wald)"
  return(combined)
}

#' Run DESeq2-LRT Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using DESeq2 with Likelihood Ratio Test.
#' Uses native create_pseudobulk_v1.
#'
#' @inheritParams runDESEQ2_Wald_v1
#' @export
runDESEQ2_LRT_v1 <- function(
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

  req <- c("Seurat", "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "limma", "dplyr", "DESeq2")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Packages missing: ", paste(miss, collapse = ", "))

  if (!exists("create_pseudobulk_v1")) stop("create_pseudobulk_v1 not found")

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
  message("Pseudobulking (DESeq2-LRT)...")
  pb <- create_pseudobulk_v1(sobj, cluster_id, sample_id, group_id)
  pb_meta <- SummarizedExperiment::colData(pb)

  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb))
  if (!is.null(keep_clusters)) cluster_names <- intersect(cluster_names, keep_clusters)

  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb)[[clust]]
    s_keep <- colSums(pb_clust) > 0
    if (sum(s_keep) < 2) next
    pb_clust <- pb_clust[, s_keep, drop = FALSE]

    sample_ids <- colnames(pb_clust)
    pb_clust_meta <- pb_meta[sample_ids, , drop = FALSE]

    # Handle Design
    if (is.null(pb_clust_meta[[group_id]])) next
    grp_col <- droplevels(factor(pb_clust_meta[[group_id]]))
    if (nlevels(grp_col) < 2) next

    dm <- data.frame(row.names = sample_ids)
    dm$group <- grp_col

    t_terms <- c("group")

    if (!is.null(batch_id) && batch_id %in% names(pb_clust_meta)) {
      dm$batch <- droplevels(factor(pb_clust_meta[[batch_id]]))
      t_terms <- c("batch", t_terms)
    }

    if (!is.null(covar_effects)) {
      for (cv in covar_effects) {
        if (cv %in% names(pb_clust_meta)) {
          dm[[cv]] <- pb_clust_meta[[cv]]
          t_terms <- c(t_terms, cv)
        }
      }
    }

    hd <- handle_design_v1(dm, t_terms)
    design_formula <- hd$formula

    # Reduced formula: All terms except 'group'
    reduced_terms <- setdiff(hd$final_terms, "group")
    if (length(reduced_terms) > 0) {
      reduced_formula <- stats::as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
    } else {
      reduced_formula <- stats::as.formula("~ 1")
    }

    pb_clust_int <- round(as.matrix(pb_clust))

    dds <- tryCatch(
      {
        DESeq2::DESeqDataSetFromMatrix(countData = pb_clust_int, colData = hd$meta, design = design_formula)
      },
      error = function(e) NULL
    )
    if (is.null(dds)) next

    dds <- tryCatch(
      {
        DESeq2::DESeq(dds, test = "LRT", reduced = reduced_formula, quiet = TRUE)
      },
      error = function(e) NULL
    )
    if (is.null(dds)) next

    res <- tryCatch(
      {
        DESeq2::results(dds)
      },
      error = function(e) NULL
    )
    if (is.null(res)) next

    res_df <- as.data.frame(res)
    res_df$cluster_id <- clust
    res_df$gene <- rownames(res_df)
    rownames(res_df) <- NULL

    if ("log2FoldChange" %in% names(res_df)) res_df$logFC <- res_df$log2FoldChange
    if ("padj" %in% names(res_df)) res_df$FDR <- res_df$padj
    if ("stat" %in% names(res_df)) res_df$statistic <- res_df$stat * sign(res_df$logFC)

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
  attr(combined, "formula") <- "Varies per cluster (LRT)"
  return(combined)
}
