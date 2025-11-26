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
runDESEQ2_Wald_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,
  covar_effects = NULL,
  contrast   = NULL,
  pb_min_cells = 3,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")

  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr","DESeq2")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0-5: runEDGER_LRT_v1과 동일한 전처리 ---
  message("0/7: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  
  required_cols <- c(cluster_id, sample_id, group_id)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  # covar_effects 컬럼 확인
  if (!is.null(covar_effects)) {
    missing_covars <- covar_effects[!covar_effects %in% colnames(meta)]
    if (length(missing_covars) > 0) {
      stop(sprintf("covar_effects에 지정된 컬럼이 없습니다: %s", paste(missing_covars, collapse=", ")))
    }
  }
  
  if (remove_na_groups) {
    na_mask <- is.na(meta[[group_id]]) | 
               is.na(meta[[cluster_id]]) | 
               is.na(meta[[sample_id]])
    if (!is.null(batch_id) && batch_id %in% colnames(meta)) {
      na_mask <- na_mask | is.na(meta[[batch_id]])
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(meta)) {
          na_mask <- na_mask | is.na(meta[[covar]])
        }
      }
    }
    
    if (is.character(meta[[group_id]])) {
      na_mask <- na_mask | (meta[[group_id]] == "NA" | meta[[group_id]] == "na")
    }
    if (is.character(meta[[cluster_id]])) {
      na_mask <- na_mask | (meta[[cluster_id]] == "NA" | meta[[cluster_id]] == "na")
    }
    if (is.character(meta[[sample_id]])) {
      na_mask <- na_mask | (meta[[sample_id]] == "NA" | meta[[sample_id]] == "na")
    }
    if (!is.null(batch_id) && batch_id %in% colnames(meta) && is.character(meta[[batch_id]])) {
      na_mask <- na_mask | (meta[[batch_id]] == "NA" | meta[[batch_id]] == "na")
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(meta) && is.character(meta[[covar]])) {
          na_mask <- na_mask | (meta[[covar]] == "NA" | meta[[covar]] == "na")
        }
      }
    }
    
    n_na_cells <- sum(na_mask)
    if (n_na_cells > 0) {
      message(sprintf("... NA 값(또는 'NA' 문자열)이 있는 %d 개의 세포를 제거합니다.", n_na_cells))
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    } else {
      message("... NA 값이 없습니다.")
    }
  }
  
  if (length(unique(meta[[group_id]])) < 2) {
    stop(sprintf("group_id ('%s')에 최소 2개의 그룹이 필요합니다. 현재: %s",
                 group_id, paste(unique(meta[[group_id]]), collapse=", ")))
  }

  message("1/7: Seurat -> SCE 변환 중...")
  sce <- Seurat::as.SingleCellExperiment(sobj)
  sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

  message("2/7: 메타데이터 정리 중...")
  sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
  sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
  sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
    sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
  }
  if (!is.null(covar_effects)) {
    for (covar in covar_effects) {
      if (covar %in% colnames(SummarizedExperiment::colData(sce))) {
        sce[[covar]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[covar]]))
      }
    }
  }

  message("3/7: Pseudobulking 중...")
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  message("4/7: Pseudobulk 메타데이터 보강 중...")
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))

  if (!"sample_id" %in% names(pb_meta)) {
    first_assay <- names(SummarizedExperiment::assays(pb))[1]
    sid_guess <- colnames(SummarizedExperiment::assays(pb)[[first_assay]])
    if (is.null(sid_guess)) stop("pb에 sample_id가 없습니다.")
    pb_meta$sample_id <- sid_guess
    rownames(pb_meta) <- pb_meta$sample_id
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)
  }

  sce_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  map_cols <- c("sample_id","group_id")
  if (!is.null(batch_id) && batch_id %in% names(sce_meta)) map_cols <- c(map_cols, batch_id)
  if (!is.null(covar_effects)) {
    covar_cols <- covar_effects[covar_effects %in% names(sce_meta)]
    if (length(covar_cols) > 0) map_cols <- c(map_cols, covar_cols)
  }
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])
  sce_map <- sce_map[complete.cases(sce_map), ]

  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  need_fix <- (!"group_id" %in% names(pb_meta)) ||
              (length(unique(pb_meta$group_id)) < 2) ||
              (all(unique(pb_meta$group_id) %in% c("type","group","group_id", NA, "")))
  need_covar_fix <- !is.null(covar_effects) && any(!covar_effects %in% names(pb_meta))
  if (need_fix || (!is.null(batch_id) && !batch_id %in% names(pb_meta)) || need_covar_fix) {
    pb_meta2 <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
    if ("group_id.x" %in% names(pb_meta2) && "group_id.y" %in% names(pb_meta2)) {
      pb_meta2$group_id <- ifelse(is.na(pb_meta2$group_id.y), pb_meta2$group_id.x, pb_meta2$group_id.y)
      pb_meta2$group_id.x <- NULL; pb_meta2$group_id.y <- NULL
    }
    rownames(pb_meta2) <- rownames(pb_meta)
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta2)
  }

  pb$sample_id <- droplevels(factor(SummarizedExperiment::colData(pb)$sample_id))
  pb$group_id  <- droplevels(factor(SummarizedExperiment::colData(pb)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb))) {
    pb[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[batch_id]]))
  }
  if (!is.null(covar_effects)) {
    for (covar in covar_effects) {
      if (covar %in% colnames(SummarizedExperiment::colData(pb))) {
        pb[[covar]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[covar]]))
      }
    }
  }

  message("5/7: Contrast 그룹 필터링 중...")
  extract_groups <- function(contrast_str, levels_available){
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
  grp_lvls <- levels(SummarizedExperiment::colData(pb)$group_id)
  tg <- extract_groups(contrast, grp_lvls)
  if (length(tg) < 2) stop(sprintf("contrast에서 추출한 그룹이 부족합니다. contrast='%s', 사용가능레벨=%s",
                                   contrast, paste(grp_lvls, collapse=", ")))

  keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
  pb_sub <- pb[, keep_idx]
  pb_sub$group_id <- droplevels(factor(SummarizedExperiment::colData(pb_sub)$group_id))

  sce_sub <- sce[, sce$sample_id %in% SummarizedExperiment::colData(pb_sub)$sample_id &
                    sce$group_id  %in% tg]
  sce_sub$cluster_id <- droplevels(factor(sce_sub$cluster_id))
  sce_sub$sample_id  <- droplevels(factor(sce_sub$sample_id))
  sce_sub$group_id   <- droplevels(factor(sce_sub$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce_sub))) {
    sce_sub[[batch_id]] <- droplevels(factor(sce_sub[[batch_id]]))
  }
  if (!is.null(covar_effects)) {
    for (covar in covar_effects) {
      if (covar %in% colnames(SummarizedExperiment::colData(sce_sub))) {
        sce_sub[[covar]] <- droplevels(factor(sce_sub[[covar]]))
      }
    }
  }

  # 6) DESeq2-Wald 분석 (클러스터별)
  message("6/7: DESeq2-Wald 분석 실행 중...")
  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb_sub))
  
  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb_sub)[[clust]]
    
    if (ncol(pb_clust) < 2) {
      message(sprintf("  클러스터 %s: 샘플 수 부족 (%d), 건너뜁니다.", clust, ncol(pb_clust)))
      next
    }
    
    # 해당 클러스터의 샘플별 메타데이터 추출
    sample_ids <- colnames(pb_clust)
    pb_clust_meta <- SummarizedExperiment::colData(pb_sub)[match(sample_ids, rownames(SummarizedExperiment::colData(pb_sub))), , drop = FALSE]
    pb_clust_group <- droplevels(factor(pb_clust_meta$group_id))
    
    if (length(levels(pb_clust_group)) < 2) {
      message(sprintf("  클러스터 %s: 그룹 수 부족 (%d), 건너뜁니다.", clust, length(levels(pb_clust_group))))
      next
    }
    
    # DESeq2를 위한 메타데이터 준비
    pb_clust_meta_df <- as.data.frame(pb_clust_meta)
    pb_clust_meta_df$group <- pb_clust_group
    rownames(pb_clust_meta_df) <- sample_ids
    
    # DESeq2 design formula
    formula_terms <- c("group")
    if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta_df)) {
      pb_clust_meta_df$batch <- droplevels(factor(pb_clust_meta_df[[batch_id]]))
      formula_terms <- c("batch", formula_terms)  # batch를 먼저 (DESeq2 convention)
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(pb_clust_meta_df)) {
          pb_clust_meta_df[[covar]] <- droplevels(factor(pb_clust_meta_df[[covar]]))
          formula_terms <- c(formula_terms, covar)
        }
      }
    }
    
    formula_str <- paste("~", paste(formula_terms, collapse = " + "))
    design_formula <- stats::as.formula(formula_str)
    
    # DESeq2 데이터셋 생성
    # pseudobulk count를 정수로 변환 (DESeq2는 정수를 요구)
    pb_clust_int <- round(as.matrix(pb_clust))
    dds <- tryCatch({
      DESeq2::DESeqDataSetFromMatrix(
        countData = pb_clust_int,
        colData = pb_clust_meta_df,
        design = design_formula
      )
    }, error = function(e) {
      message(sprintf("  클러스터 %s: DESeqDataSetFromMatrix 실패: %s", clust, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(dds)) next
    
    # DESeq2 분석 (Wald test)
    dds <- tryCatch({
      DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
    }, error = function(e) {
      message(sprintf("  클러스터 %s: DESeq 실행 실패: %s", clust, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(dds)) next
    
    # 결과 추출 (contrast: 마지막 레벨 - 첫 레벨)
    grp_levels <- levels(pb_clust_group)
    if (length(grp_levels) >= 2) {
      # contrast에서 그룹 추출
      contrast_groups <- extract_groups(contrast, grp_levels)
      if (length(contrast_groups) >= 2) {
        res <- tryCatch({
          DESeq2::results(dds, contrast = c("group", contrast_groups[2], contrast_groups[1]))
        }, error = function(e) {
          message(sprintf("  클러스터 %s: results 추출 실패: %s", clust, conditionMessage(e)))
          return(NULL)
        })
      } else {
        # 기본: 마지막 레벨 vs 첫 레벨
        res <- tryCatch({
          DESeq2::results(dds, contrast = c("group", tail(grp_levels, 1), grp_levels[1]))
        }, error = function(e) {
          message(sprintf("  클러스터 %s: results 추출 실패: %s", clust, conditionMessage(e)))
          return(NULL)
        })
      }
    } else {
      next
    }
    
    if (is.null(res)) next
    
    # 결과를 데이터프레임으로 변환
    res_df <- as.data.frame(res)
    res_df$cluster_id <- clust
    res_df$gene <- rownames(res_df)
    rownames(res_df) <- NULL
    
    # 컬럼명 표준화
    if ("log2FoldChange" %in% colnames(res_df)) {
      res_df$logFC <- res_df$log2FoldChange
    }
    if ("pvalue" %in% colnames(res_df)) {
      res_df$pvalue <- res_df$pvalue
    }
    if ("padj" %in% colnames(res_df)) {
      res_df$FDR <- res_df$padj
    }
    if ("stat" %in% colnames(res_df)) {
      res_df$statistic <- res_df$stat
    }
    
    all_results[[clust]] <- res_df
  }
  
  if (length(all_results) == 0) {
    stop("모든 클러스터에서 분석 실패")
  }
  
  combined <- do.call(rbind, all_results)
  
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  message("7/7: DESeq2-Wald 분석 완료.")
  return(combined)
}

#' Run DESeq2-LRT Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using DESeq2 with Likelihood Ratio Test.
#'
#' @inheritParams runDESEQ2_Wald_v1
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runDESEQ2_LRT_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,
  covar_effects = NULL,
  contrast   = NULL,
  pb_min_cells = 3,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")

  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr","DESeq2")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0-5: runDESEQ2_Wald_v1과 동일한 전처리 ---
  message("0/7: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  
  required_cols <- c(cluster_id, sample_id, group_id)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  # covar_effects 컬럼 확인
  if (!is.null(covar_effects)) {
    missing_covars <- covar_effects[!covar_effects %in% colnames(meta)]
    if (length(missing_covars) > 0) {
      stop(sprintf("covar_effects에 지정된 컬럼이 없습니다: %s", paste(missing_covars, collapse=", ")))
    }
  }
  
  if (remove_na_groups) {
    na_mask <- is.na(meta[[group_id]]) | 
               is.na(meta[[cluster_id]]) | 
               is.na(meta[[sample_id]])
    if (!is.null(batch_id) && batch_id %in% colnames(meta)) {
      na_mask <- na_mask | is.na(meta[[batch_id]])
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(meta)) {
          na_mask <- na_mask | is.na(meta[[covar]])
        }
      }
    }
    
    if (is.character(meta[[group_id]])) {
      na_mask <- na_mask | (meta[[group_id]] == "NA" | meta[[group_id]] == "na")
    }
    if (is.character(meta[[cluster_id]])) {
      na_mask <- na_mask | (meta[[cluster_id]] == "NA" | meta[[cluster_id]] == "na")
    }
    if (is.character(meta[[sample_id]])) {
      na_mask <- na_mask | (meta[[sample_id]] == "NA" | meta[[sample_id]] == "na")
    }
    if (!is.null(batch_id) && batch_id %in% colnames(meta) && is.character(meta[[batch_id]])) {
      na_mask <- na_mask | (meta[[batch_id]] == "NA" | meta[[batch_id]] == "na")
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(meta) && is.character(meta[[covar]])) {
          na_mask <- na_mask | (meta[[covar]] == "NA" | meta[[covar]] == "na")
        }
      }
    }
    
    n_na_cells <- sum(na_mask)
    if (n_na_cells > 0) {
      message(sprintf("... NA 값(또는 'NA' 문자열)이 있는 %d 개의 세포를 제거합니다.", n_na_cells))
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    } else {
      message("... NA 값이 없습니다.")
    }
  }
  
  if (length(unique(meta[[group_id]])) < 2) {
    stop(sprintf("group_id ('%s')에 최소 2개의 그룹이 필요합니다. 현재: %s",
                 group_id, paste(unique(meta[[group_id]]), collapse=", ")))
  }

  message("1/7: Seurat -> SCE 변환 중...")
  sce <- Seurat::as.SingleCellExperiment(sobj)
  sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

  message("2/7: 메타데이터 정리 중...")
  sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
  sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
  sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
    sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
  }
  if (!is.null(covar_effects)) {
    for (covar in covar_effects) {
      if (covar %in% colnames(SummarizedExperiment::colData(sce))) {
        sce[[covar]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[covar]]))
      }
    }
  }

  message("3/7: Pseudobulking 중...")
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  message("4/7: Pseudobulk 메타데이터 보강 중...")
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))

  if (!"sample_id" %in% names(pb_meta)) {
    first_assay <- names(SummarizedExperiment::assays(pb))[1]
    sid_guess <- colnames(SummarizedExperiment::assays(pb)[[first_assay]])
    if (is.null(sid_guess)) stop("pb에 sample_id가 없습니다.")
    pb_meta$sample_id <- sid_guess
    rownames(pb_meta) <- pb_meta$sample_id
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)
  }

  sce_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  map_cols <- c("sample_id","group_id")
  if (!is.null(batch_id) && batch_id %in% names(sce_meta)) map_cols <- c(map_cols, batch_id)
  if (!is.null(covar_effects)) {
    covar_cols <- covar_effects[covar_effects %in% names(sce_meta)]
    if (length(covar_cols) > 0) map_cols <- c(map_cols, covar_cols)
  }
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])
  sce_map <- sce_map[complete.cases(sce_map), ]

  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  need_fix <- (!"group_id" %in% names(pb_meta)) ||
              (length(unique(pb_meta$group_id)) < 2) ||
              (all(unique(pb_meta$group_id) %in% c("type","group","group_id", NA, "")))
  need_covar_fix <- !is.null(covar_effects) && any(!covar_effects %in% names(pb_meta))
  if (need_fix || (!is.null(batch_id) && !batch_id %in% names(pb_meta)) || need_covar_fix) {
    pb_meta2 <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
    if ("group_id.x" %in% names(pb_meta2) && "group_id.y" %in% names(pb_meta2)) {
      pb_meta2$group_id <- ifelse(is.na(pb_meta2$group_id.y), pb_meta2$group_id.x, pb_meta2$group_id.y)
      pb_meta2$group_id.x <- NULL; pb_meta2$group_id.y <- NULL
    }
    rownames(pb_meta2) <- rownames(pb_meta)
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta2)
  }

  pb$sample_id <- droplevels(factor(SummarizedExperiment::colData(pb)$sample_id))
  pb$group_id  <- droplevels(factor(SummarizedExperiment::colData(pb)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb))) {
    pb[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[batch_id]]))
  }
  if (!is.null(covar_effects)) {
    for (covar in covar_effects) {
      if (covar %in% colnames(SummarizedExperiment::colData(pb))) {
        pb[[covar]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[covar]]))
      }
    }
  }

  message("5/7: Contrast 그룹 필터링 중...")
  extract_groups <- function(contrast_str, levels_available){
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
  grp_lvls <- levels(SummarizedExperiment::colData(pb)$group_id)
  tg <- extract_groups(contrast, grp_lvls)
  if (length(tg) < 2) stop(sprintf("contrast에서 추출한 그룹이 부족합니다. contrast='%s', 사용가능레벨=%s",
                                   contrast, paste(grp_lvls, collapse=", ")))

  keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
  pb_sub <- pb[, keep_idx]
  pb_sub$group_id <- droplevels(factor(SummarizedExperiment::colData(pb_sub)$group_id))

  sce_sub <- sce[, sce$sample_id %in% SummarizedExperiment::colData(pb_sub)$sample_id &
                    sce$group_id  %in% tg]
  sce_sub$cluster_id <- droplevels(factor(sce_sub$cluster_id))
  sce_sub$sample_id  <- droplevels(factor(sce_sub$sample_id))
  sce_sub$group_id   <- droplevels(factor(sce_sub$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce_sub))) {
    sce_sub[[batch_id]] <- droplevels(factor(sce_sub[[batch_id]]))
  }
  if (!is.null(covar_effects)) {
    for (covar in covar_effects) {
      if (covar %in% colnames(SummarizedExperiment::colData(sce_sub))) {
        sce_sub[[covar]] <- droplevels(factor(sce_sub[[covar]]))
      }
    }
  }

  # 6) DESeq2-LRT 분석 (클러스터별)
  message("6/7: DESeq2-LRT 분석 실행 중...")
  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb_sub))
  
  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb_sub)[[clust]]
    
    if (ncol(pb_clust) < 2) {
      message(sprintf("  클러스터 %s: 샘플 수 부족 (%d), 건너뜁니다.", clust, ncol(pb_clust)))
      next
    }
    
    sample_ids <- colnames(pb_clust)
    pb_clust_meta <- SummarizedExperiment::colData(pb_sub)[match(sample_ids, rownames(SummarizedExperiment::colData(pb_sub))), , drop = FALSE]
    pb_clust_group <- droplevels(factor(pb_clust_meta$group_id))
    
    if (length(levels(pb_clust_group)) < 2) {
      message(sprintf("  클러스터 %s: 그룹 수 부족 (%d), 건너뜁니다.", clust, length(levels(pb_clust_group))))
      next
    }
    
    pb_clust_meta_df <- as.data.frame(pb_clust_meta)
    pb_clust_meta_df$group <- pb_clust_group
    rownames(pb_clust_meta_df) <- sample_ids
    
    # Design formula 구성 (LRT는 reduced model도 필요)
    formula_terms <- c("group")
    reduced_terms <- character()
    
    if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta_df)) {
      pb_clust_meta_df$batch <- droplevels(factor(pb_clust_meta_df[[batch_id]]))
      reduced_terms <- c(reduced_terms, "batch")
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(pb_clust_meta_df)) {
          pb_clust_meta_df[[covar]] <- droplevels(factor(pb_clust_meta_df[[covar]]))
          formula_terms <- c(formula_terms, covar)
          reduced_terms <- c(reduced_terms, covar)
        }
      }
    }
    
    # Full model: batch (if present) + group + covar_effects
    if (length(reduced_terms) > 0) {
      formula_terms <- c(reduced_terms, formula_terms)
    }
    formula_str <- paste("~", paste(formula_terms, collapse = " + "))
    design_formula <- stats::as.formula(formula_str)
    
    # Reduced model: batch (if present) + covar_effects (group 제외)
    if (length(reduced_terms) > 0) {
      reduced_str <- paste("~", paste(reduced_terms, collapse = " + "))
    } else {
      reduced_str <- "~ 1"
    }
    reduced_formula <- stats::as.formula(reduced_str)
    
    # pseudobulk count를 정수로 변환 (DESeq2는 정수를 요구)
    pb_clust_int <- round(as.matrix(pb_clust))
    dds <- tryCatch({
      DESeq2::DESeqDataSetFromMatrix(
        countData = pb_clust_int,
        colData = pb_clust_meta_df,
        design = design_formula
      )
    }, error = function(e) {
      message(sprintf("  클러스터 %s: DESeqDataSetFromMatrix 실패: %s", clust, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(dds)) next
    
    # DESeq2 분석 (LRT)
    dds <- tryCatch({
      DESeq2::DESeq(dds, test = "LRT", reduced = reduced_formula, quiet = TRUE)
    }, error = function(e) {
      message(sprintf("  클러스터 %s: DESeq 실행 실패: %s", clust, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(dds)) next
    
    # LRT 결과 추출
    res <- tryCatch({
      DESeq2::results(dds)
    }, error = function(e) {
      message(sprintf("  클러스터 %s: results 추출 실패: %s", clust, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(res)) next
    
    res_df <- as.data.frame(res)
    res_df$cluster_id <- clust
    res_df$gene <- rownames(res_df)
    rownames(res_df) <- NULL
    
    if ("log2FoldChange" %in% colnames(res_df)) {
      res_df$logFC <- res_df$log2FoldChange
    }
    if ("pvalue" %in% colnames(res_df)) {
      res_df$pvalue <- res_df$pvalue
    }
    if ("padj" %in% colnames(res_df)) {
      res_df$FDR <- res_df$padj
    }
    if ("stat" %in% colnames(res_df)) {
      res_df$statistic <- res_df$stat
    }
    
    all_results[[clust]] <- res_df
  }
  
  if (length(all_results) == 0) {
    stop("모든 클러스터에서 분석 실패")
  }
  
  combined <- do.call(rbind, all_results)
  
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  message("7/7: DESeq2-LRT 분석 완료.")
  return(combined)
}

