# ============================================================================
# Limma-based DEG Methods (Phase 2)
# ============================================================================
# runMUSCAT 스타일로 표준화된 limma 계열 방법론
# ============================================================================

#' Run limma-voom Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using limma with voom transformation.
#' This function follows the same interface as runMUSCAT for consistency.
#'
#' @param sobj Seurat object
#' @param cluster_id Column name for cell type/cluster (default: "seurat_clusters")
#' @param sample_id Column name for sample ID (default: "hos_no")
#' @param group_id Column name for group/condition (default: "type")
#' @param batch_id Optional column name for batch variable
#' @param covar_effects Optional character vector of column names for additional covariates 
#'   to include as fixed effects (e.g., c("sex"))
#' @param contrast Contrast string (e.g., "IS - SAH")
#' @param pb_min_cells Minimum cells per pseudobulk sample (default: 3)
#' @param keep_clusters Optional vector of cluster IDs to keep
#' @param cluster_label_map Optional named vector mapping cluster IDs to labels
#' @param remove_na_groups Remove cells with NA in group_id (default: TRUE)
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runLIMMA_voom_v1 <- function(
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

  # deps
  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr","edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0. NA 값 처리 (runMUSCAT과 동일) ---
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

  # 1) Seurat -> SCE, prepSCE
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

  # 2) Pseudobulk
  message("3/7: Pseudobulking 중...")
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  # 2-1) pb 메타 보강
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
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])
  sce_map <- sce_map[complete.cases(sce_map), ]

  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  need_fix <- (!"group_id" %in% names(pb_meta)) ||
              (length(unique(pb_meta$group_id)) < 2) ||
              (all(unique(pb_meta$group_id) %in% c("type","group","group_id", NA, "")))
  if (need_fix || (!is.null(batch_id) && !batch_id %in% names(pb_meta))) {
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

  # 3) contrast 그룹만 자동 subset
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

  # 4) contrast 파싱 함수 (공통)
  fix_contrast <- function(contrast_str, design_cols){
    z <- gsub("\\s+", "", contrast_str)
    toks <- unlist(strsplit(z, "([+\\-])", perl=TRUE))
    ops  <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl=TRUE)))
    rebuild <- function(tok){
      tok <- gsub("^group(_id)?", "group", tok)
      if (!grepl("^group", tok)) tok <- paste0("group", tok)
      tok
    }
    toks2 <- vapply(toks, rebuild, character(1))
    out <- toks2[1]; if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i+1])
    out
  }

  # 5) limma-voom 분석 (클러스터별)
  message("6/7: limma-voom 분석 실행 중...")
  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb_sub))
  
  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb_sub)[[clust]]
    
    # 최소 샘플 수 확인
    if (ncol(pb_clust) < 2) {
      message(sprintf("  클러스터 %s: 샘플 수 부족 (%d), 건너뜁니다.", clust, ncol(pb_clust)))
      next
    }
    
    # 해당 클러스터의 샘플별 메타데이터 추출
    sample_ids <- colnames(pb_clust)
    pb_clust_meta <- SummarizedExperiment::colData(pb_sub)[match(sample_ids, rownames(SummarizedExperiment::colData(pb_sub))), , drop = FALSE]
    pb_clust_group <- droplevels(factor(pb_clust_meta$group_id))
    
    # 그룹 수 확인
    if (length(levels(pb_clust_group)) < 2) {
      message(sprintf("  클러스터 %s: 그룹 수 부족 (%d), 건너뜁니다.", clust, length(levels(pb_clust_group))))
      next
    }
    
    # 클러스터별 design matrix 생성
    pb_clust_meta$group <- pb_clust_group
    pb_clust_meta_df <- as.data.frame(pb_clust_meta)
    
    # Design formula 구성
    formula_terms <- c("group")
    if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta_df)) {
      pb_clust_meta_df$batch <- droplevels(factor(pb_clust_meta_df[[batch_id]]))
      formula_terms <- c(formula_terms, "batch")
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(pb_clust_meta_df)) {
          pb_clust_meta_df[[covar]] <- droplevels(factor(pb_clust_meta_df[[covar]]))
          formula_terms <- c(formula_terms, covar)
        }
      }
    }
    
    formula_str <- paste("~ 0 +", paste(formula_terms, collapse = " + "))
    design_clust <- stats::model.matrix(as.formula(formula_str), data = pb_clust_meta_df)
    
    # 클러스터별 contrast matrix 생성
    contrast_fixed <- fix_contrast(contrast, colnames(design_clust))
    contrast_matrix_clust <- limma::makeContrasts(contrasts = contrast_fixed, levels = design_clust)
    
    # Filter samples with zero library size BEFORE creating DGEList
    lib_sizes <- colSums(pb_clust)
    valid_samples <- lib_sizes > 0
    if (sum(valid_samples) < length(valid_samples)) {
      pb_clust <- pb_clust[, valid_samples, drop = FALSE]
      pb_clust_group <- pb_clust_group[valid_samples]
      pb_clust_meta <- pb_clust_meta[valid_samples, , drop = FALSE]
      # Recreate design matrix after filtering
      pb_clust_meta$group <- pb_clust_group
      pb_clust_meta_df <- as.data.frame(pb_clust_meta)
      
      # Rebuild formula terms
      formula_terms <- c("group")
      use_batch <- FALSE
      if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta_df)) {
        pb_clust_meta_df$batch <- droplevels(factor(pb_clust_meta_df[[batch_id]]))
        if (length(unique(pb_clust_meta_df$batch)) > 1 && length(unique(pb_clust_group)) > 1) {
          batch_group_table <- table(pb_clust_meta_df$batch, pb_clust_group)
          row_sums <- rowSums(batch_group_table > 0)
          col_sums <- colSums(batch_group_table > 0)
          if (!(all(row_sums == 1) || all(col_sums == 1))) {
            use_batch <- TRUE
            formula_terms <- c(formula_terms, "batch")
          }
        }
      }
      if (!is.null(covar_effects)) {
        for (covar in covar_effects) {
          if (covar %in% colnames(pb_clust_meta_df)) {
            pb_clust_meta_df[[covar]] <- droplevels(factor(pb_clust_meta_df[[covar]]))
            formula_terms <- c(formula_terms, covar)
          }
        }
      }
      
      formula_str <- paste("~ 0 +", paste(formula_terms, collapse = " + "))
      design_clust <- tryCatch({
        stats::model.matrix(as.formula(formula_str), data = pb_clust_meta_df)
      }, error = function(e) {
        stats::model.matrix(~ 0 + group, data = pb_clust_meta_df)
      })
      contrast_fixed <- fix_contrast(contrast, colnames(design_clust))
      contrast_matrix_clust <- limma::makeContrasts(contrasts = contrast_fixed, levels = design_clust)
    }
    
    if (ncol(pb_clust) < 2 || length(unique(pb_clust_group)) < 2) {
      message(sprintf("  클러스터 %s: 유효한 샘플이 부족합니다", clust))
      next
    }
    
    # DGEList 생성
    dge <- edgeR::DGEList(counts = pb_clust, group = pb_clust_group)
    
    # 필터링
    keep <- edgeR::filterByExpr(dge, group = pb_clust_group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    if (nrow(dge) == 0) {
      message(sprintf("  클러스터 %s: 필터링 후 유전자 없음, 건너뜁니다.", clust))
      next
    }
    
    # 정규화
    dge <- edgeR::calcNormFactors(dge)
    
    # voom 변환
    v <- limma::voom(dge, design_clust, plot = FALSE)
    
    # limma 분석
    fit <- limma::lmFit(v, design_clust)
    fit <- limma::contrasts.fit(fit, contrast_matrix_clust)
    fit <- limma::eBayes(fit)
    
    # 결과 추출
    res <- limma::topTable(fit, number = Inf, sort.by = "none")
    
    # 결과에 클러스터 정보 추가
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    all_results[[clust]] <- res
  }
  
  # 결과 결합
  if (length(all_results) == 0) {
    stop("모든 클러스터에서 분석 실패")
  }
  
  combined <- do.call(rbind, all_results)
  
  # 컬럼명 표준화 (runMUSCAT과 유사하게)
  if ("logFC" %in% colnames(combined)) {
    combined$logFC <- combined$logFC
  }
  if ("P.Value" %in% colnames(combined)) {
    combined$pvalue <- combined$P.Value
  }
  if ("adj.P.Val" %in% colnames(combined)) {
    combined$FDR <- combined$adj.P.Val
  }
  if ("t" %in% colnames(combined)) {
    combined$statistic <- combined$t
  }
  
  # cluster_label 추가
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  message("limma-voom 분석 완료.")
  return(combined)
}

#' Run limma-trend Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using limma with trend method.
#' This function follows the same interface as runMUSCAT for consistency.
#'
#' @inheritParams runLIMMA_voom_v1
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runLIMMA_trend_v1 <- function(
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

  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr","edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0-6: runLIMMA_voom_v1과 동일한 전처리 ---
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

  # 4) contrast 파싱 함수 (공통)
  fix_contrast <- function(contrast_str, design_cols){
    z <- gsub("\\s+", "", contrast_str)
    toks <- unlist(strsplit(z, "([+\\-])", perl=TRUE))
    ops  <- unlist(regmatches(z, gregexpr("([+\\-])", z, perl=TRUE)))
    rebuild <- function(tok){
      tok <- gsub("^group(_id)?", "group", tok)
      if (!grepl("^group", tok)) tok <- paste0("group", tok)
      tok
    }
    toks2 <- vapply(toks, rebuild, character(1))
    out <- toks2[1]; if (length(ops)) for (i in seq_along(ops)) out <- paste0(out, ops[i], toks2[i+1])
    out
  }

  message("6/7: limma-trend 분석 실행 중...")
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
    
    # 그룹 수 확인
    if (length(levels(pb_clust_group)) < 2) {
      message(sprintf("  클러스터 %s: 그룹 수 부족 (%d), 건너뜁니다.", clust, length(levels(pb_clust_group))))
      next
    }
    
    # 클러스터별 design matrix 생성
    pb_clust_meta$group <- pb_clust_group
    pb_clust_meta_df <- as.data.frame(pb_clust_meta)
    
    # Design formula 구성
    formula_terms <- c("group")
    if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta_df)) {
      pb_clust_meta_df$batch <- droplevels(factor(pb_clust_meta_df[[batch_id]]))
      formula_terms <- c(formula_terms, "batch")
    }
    if (!is.null(covar_effects)) {
      for (covar in covar_effects) {
        if (covar %in% colnames(pb_clust_meta_df)) {
          pb_clust_meta_df[[covar]] <- droplevels(factor(pb_clust_meta_df[[covar]]))
          formula_terms <- c(formula_terms, covar)
        }
      }
    }
    
    formula_str <- paste("~ 0 +", paste(formula_terms, collapse = " + "))
    design_clust <- stats::model.matrix(as.formula(formula_str), data = pb_clust_meta_df)
    
    # 클러스터별 contrast matrix 생성
    contrast_fixed <- fix_contrast(contrast, colnames(design_clust))
    contrast_matrix_clust <- limma::makeContrasts(contrasts = contrast_fixed, levels = design_clust)
    
    dge <- edgeR::DGEList(counts = pb_clust)
    keep <- edgeR::filterByExpr(dge, group = pb_clust_group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    if (nrow(dge) == 0) {
      message(sprintf("  클러스터 %s: 필터링 후 유전자 없음, 건너뜁니다.", clust))
      next
    }
    
    dge <- edgeR::calcNormFactors(dge)
    
    # logCPM 변환 (trend method)
    logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    
    # limma-trend 분석
    fit <- limma::lmFit(logCPM, design_clust)
    fit <- limma::contrasts.fit(fit, contrast_matrix_clust)
    fit <- limma::eBayes(fit, trend = TRUE)
    
    res <- limma::topTable(fit, number = Inf, sort.by = "none")
    
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    all_results[[clust]] <- res
  }
  
  if (length(all_results) == 0) {
    stop("모든 클러스터에서 분석 실패")
  }
  
  combined <- do.call(rbind, all_results)
  
  if ("logFC" %in% colnames(combined)) {
    combined$logFC <- combined$logFC
  }
  if ("P.Value" %in% colnames(combined)) {
    combined$pvalue <- combined$P.Value
  }
  if ("adj.P.Val" %in% colnames(combined)) {
    combined$FDR <- combined$adj.P.Val
  }
  if ("t" %in% colnames(combined)) {
    combined$statistic <- combined$t
  }
  
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  message("limma-trend 분석 완료.")
  return(combined)
}

