# ============================================================================
# edgeR-based DEG Methods (Phase 2)
# ============================================================================
# runMUSCAT2_v1 스타일로 표준화된 edgeR 계열 방법론
# ============================================================================

#' Run edgeR-LRT Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using edgeR with Likelihood Ratio Test.
#' This function follows the same interface as runMUSCAT2_v1 for consistency.
#'
#' @inheritParams runLIMMA_voom_v1
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runEDGER_LRT_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,
  contrast   = NULL,
  pb_min_cells = 3,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE,
  min_samples_per_group = 2
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")

  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr","edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0-3: runLIMMA_voom_v1과 동일한 전처리 ---
  message("0/7: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  
  required_cols <- c(cluster_id, sample_id, group_id)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  if (remove_na_groups) {
    na_mask <- is.na(meta[[group_id]]) | 
               is.na(meta[[cluster_id]]) | 
               is.na(meta[[sample_id]])
    if (!is.null(batch_id) && batch_id %in% colnames(meta)) {
      na_mask <- na_mask | is.na(meta[[batch_id]])
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

  # 4) contrast 파싱 함수
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

  # 5) edgeR-LRT 분석 (클러스터별)
  message("6/7: edgeR-LRT 분석 실행 중...")
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
    grp_counts <- table(pb_clust_group)
    if (any(grp_counts < min_samples_per_group)) {
      message(sprintf("  클러스터 %s: 그룹별 샘플 수 부족 (%s), 건너뜁니다.",
                      clust,
                      paste(sprintf("%s=%d", names(grp_counts), grp_counts), collapse = ", ")))
      next
    }
    
    # 클러스터별 design matrix 생성
    pb_clust_meta$group <- pb_clust_group
    if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta)) {
      pb_clust_meta$batch <- droplevels(factor(pb_clust_meta[[batch_id]]))
      design_clust <- stats::model.matrix(~ 0 + group + batch,
                                         data = as.data.frame(pb_clust_meta))
    } else {
      design_clust <- stats::model.matrix(~ 0 + group,
                                         data = as.data.frame(pb_clust_meta))
    }
    
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
      # Recreate design matrix after filtering (re-check batch confounding)
      pb_clust_meta$group <- pb_clust_group
      use_batch <- FALSE
      if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta)) {
        pb_clust_meta$batch <- droplevels(factor(pb_clust_meta[[batch_id]]))
        if (length(unique(pb_clust_meta$batch)) > 1 && length(unique(pb_clust_group)) > 1) {
          batch_group_table <- table(pb_clust_meta$batch, pb_clust_group)
          row_sums <- rowSums(batch_group_table > 0)
          col_sums <- colSums(batch_group_table > 0)
          if (!(all(row_sums == 1) || all(col_sums == 1))) {
            use_batch <- TRUE
          }
        }
      }
      design_clust <- tryCatch({
        if (use_batch) {
          stats::model.matrix(~ 0 + group + batch, data = as.data.frame(pb_clust_meta))
        } else {
          stats::model.matrix(~ 0 + group, data = as.data.frame(pb_clust_meta))
        }
      }, error = function(e) {
        stats::model.matrix(~ 0 + group, data = as.data.frame(pb_clust_meta))
      })
      contrast_fixed <- fix_contrast(contrast, colnames(design_clust))
      contrast_matrix_clust <- limma::makeContrasts(contrasts = contrast_fixed, levels = design_clust)
      grp_counts <- table(pb_clust_group)
      if (length(grp_counts) < 2 || any(grp_counts < min_samples_per_group)) {
        message(sprintf("  클러스터 %s: 유효 샘플 재확인 후에도 부족, 건너뜁니다.", clust))
        next
      }
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
    
    # Dispersion 추정
    dge <- edgeR::estimateDisp(dge, design_clust)
    
    # LRT (Likelihood Ratio Test)
    fit <- edgeR::glmFit(dge, design_clust)
    lrt <- edgeR::glmLRT(fit, contrast = contrast_matrix_clust)
    
    # 결과 추출
    res <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
    
    # 결과에 클러스터 정보 추가
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    # 컬럼명 표준화
    if ("logFC" %in% colnames(res)) {
      res$logFC <- res$logFC
    }
    if ("PValue" %in% colnames(res)) {
      res$pvalue <- res$PValue
    }
    if ("FDR" %in% colnames(res)) {
      res$FDR <- res$FDR
    }
    if ("LR" %in% colnames(res)) {
      res$statistic <- res$LR
    }
    
    all_results[[clust]] <- res
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

  message("7/7: edgeR-LRT 분석 완료.")
  return(combined)
}

#' Run edgeR-QLF Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using edgeR with Quasi-Likelihood F-test.
#'
#' @inheritParams runEDGER_LRT_v1
#'
#' @return Data frame with differential expression results per cluster
#'
#' @export
runEDGER_QLF_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,
  contrast   = NULL,
  pb_min_cells = 3,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE,
  min_samples_per_group = 2
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")

  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr","edgeR")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0-5: runEDGER_LRT_v1과 동일한 전처리 (간단히 함수로 추출 가능하지만 일단 복사) ---
  message("0/7: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  
  required_cols <- c(cluster_id, sample_id, group_id)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  if (remove_na_groups) {
    na_mask <- is.na(meta[[group_id]]) | 
               is.na(meta[[cluster_id]]) | 
               is.na(meta[[sample_id]])
    if (!is.null(batch_id) && batch_id %in% colnames(meta)) {
      na_mask <- na_mask | is.na(meta[[batch_id]])
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

  # 5) edgeR-QLF 분석 (클러스터별)
  message("6/7: edgeR-QLF 분석 실행 중...")
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
    grp_counts <- table(pb_clust_group)
    if (any(grp_counts < min_samples_per_group)) {
      message(sprintf("  클러스터 %s: 그룹별 샘플 수 부족 (%s), 건너뜁니다.",
                      clust,
                      paste(sprintf("%s=%d", names(grp_counts), grp_counts), collapse = ", ")))
      next
    }
    
    pb_clust_meta$group <- pb_clust_group
    if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta)) {
      pb_clust_meta$batch <- droplevels(factor(pb_clust_meta[[batch_id]]))
      design_clust <- stats::model.matrix(~ 0 + group + batch,
                                         data = as.data.frame(pb_clust_meta))
    } else {
      design_clust <- stats::model.matrix(~ 0 + group,
                                         data = as.data.frame(pb_clust_meta))
    }
    
    contrast_fixed <- fix_contrast(contrast, colnames(design_clust))
    contrast_matrix_clust <- limma::makeContrasts(contrasts = contrast_fixed, levels = design_clust)
    
    lib_sizes <- colSums(pb_clust)
    valid_samples <- lib_sizes > 0
    if (sum(valid_samples) < length(valid_samples)) {
      pb_clust <- pb_clust[, valid_samples, drop = FALSE]
      pb_clust_group <- pb_clust_group[valid_samples]
      pb_clust_meta <- pb_clust_meta[valid_samples, , drop = FALSE]
      pb_clust_meta$group <- pb_clust_group
      if (!is.null(batch_id) && batch_id %in% colnames(pb_clust_meta)) {
        pb_clust_meta$batch <- droplevels(factor(pb_clust_meta[[batch_id]]))
        design_clust <- stats::model.matrix(~ 0 + group + batch,
                                            data = as.data.frame(pb_clust_meta))
      } else {
        design_clust <- stats::model.matrix(~ 0 + group,
                                            data = as.data.frame(pb_clust_meta))
      }
      contrast_fixed <- fix_contrast(contrast, colnames(design_clust))
      contrast_matrix_clust <- limma::makeContrasts(contrasts = contrast_fixed, levels = design_clust)
      grp_counts <- table(pb_clust_group)
      if (length(grp_counts) < 2 || any(grp_counts < min_samples_per_group)) {
        message(sprintf("  클러스터 %s: 유효 샘플 재확인 후에도 부족, 건너뜁니다.", clust))
        next
      }
    }
    
    dge <- edgeR::DGEList(counts = pb_clust)
    keep <- edgeR::filterByExpr(dge, group = pb_clust_group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    if (nrow(dge) == 0) {
      message(sprintf("  클러스터 %s: 필터링 후 유전자 없음, 건너뜁니다.", clust))
      next
    }
    
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design_clust)
    
    # QLF (Quasi-Likelihood F-test)
    fit <- edgeR::glmQLFit(dge, design_clust)
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_matrix_clust)
    
    res <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    
    res$cluster_id <- clust
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    if ("logFC" %in% colnames(res)) {
      res$logFC <- res$logFC
    }
    if ("PValue" %in% colnames(res)) {
      res$pvalue <- res$PValue
    }
    if ("FDR" %in% colnames(res)) {
      res$FDR <- res$FDR
    }
    if ("F" %in% colnames(res)) {
      res$statistic <- res$F
    }
    
    all_results[[clust]] <- res
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

  message("7/7: edgeR-QLF 분석 완료.")
  return(combined)
}

