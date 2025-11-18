# ============================================================================
# Dream-based DEG Method (variancePartition)
# ============================================================================

#' Run Dream (variancePartition) Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using the dream workflow from the
#' variancePartition package. This method models sample-level random effects
#' (e.g., patient IDs) via linear mixed models on pseudobulk counts.
#'
#' @inheritParams runLIMMA_voom_v1
#' @param patient_col Column name to use for the random effect (default: `sample_id`)
#' @param min_samples_per_group Minimum number of pseudobulk samples required
#'   per contrast group (default: 2). Clusters that do not meet this requirement
#'   are skipped.
#'
#' @return Data frame of differential expression results per cluster.
#' @export
runDREAM_v1 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,
  patient_col = sample_id,
  contrast   = NULL,
  pb_min_cells = 3,
  keep_clusters = NULL,
  cluster_label_map = NULL,
  remove_na_groups = TRUE,
  min_samples_per_group = 2
) {
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")
  
  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment",
           "S4Vectors","limma","edgeR","variancePartition")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse = ", "))
  
  message("0/7: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  required_cols <- c(cluster_id, sample_id, group_id, patient_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  if (remove_na_groups) {
    na_mask <- is.na(meta[[group_id]]) |
      is.na(meta[[cluster_id]]) |
      is.na(meta[[sample_id]]) |
      is.na(meta[[patient_col]])
    if (!is.null(batch_id) && batch_id %in% colnames(meta)) {
      na_mask <- na_mask | is.na(meta[[batch_id]])
    }
    n_na <- sum(na_mask)
    if (n_na > 0) {
      message(sprintf("... NA 값이 있는 %d 개의 세포를 제거합니다.", n_na))
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    } else {
      message("... NA 값이 없습니다.")
    }
  }
  
  if (length(unique(meta[[group_id]])) < 2) {
    stop(sprintf("group_id ('%s')에 최소 2개의 그룹이 필요합니다.", group_id))
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
  
normalize_sample_column <- function(df) {
    if (!"sample_id" %in% names(df)) {
      sample_cols <- grep("^sample_id", names(df), value = TRUE)
      if (length(sample_cols) > 0) {
        df$sample_id <- df[[sample_cols[1]]]
      }
    }
    df
  }
  
  sample_lookup <- data.frame(
    sample_id = as.character(meta[[sample_id]]),
    group_id = meta[[group_id]],
    stringsAsFactors = FALSE
  )
  sample_lookup[[patient_col]] <- meta[[patient_col]]
  if (!is.null(batch_id)) {
    sample_lookup[[batch_id]] <- meta[[batch_id]]
  }
  sample_lookup <- sample_lookup[!is.na(sample_lookup$sample_id), , drop = FALSE]
  sample_lookup$sample_id <- as.character(sample_lookup$sample_id)
  sample_lookup <- dplyr::distinct(sample_lookup, sample_id, .keep_all = TRUE)
  
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
  if (!patient_col %in% map_cols) map_cols <- c(map_cols, patient_col)
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])
  sce_map <- sce_map[complete.cases(sce_map), ]
  
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  pb_meta <- normalize_sample_column(pb_meta)
  pb_meta2 <- dplyr::left_join(pb_meta, sce_map, by = "sample_id")
  pb_meta2 <- normalize_sample_column(pb_meta2)
  if (!patient_col %in% names(pb_meta2)) {
    if ("sample_id" %in% names(pb_meta2)) {
      pb_meta2[[patient_col]] <- pb_meta2$sample_id
    } else {
      pb_meta2[[patient_col]] <- rownames(pb_meta2)
    }
  }
  rownames(pb_meta2) <- rownames(pb_meta)
  SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta2)
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
  if (length(tg) < 2) stop(sprintf("contrast에서 추출한 그룹이 부족합니다. contrast='%s'", contrast))
  
  keep_idx <- SummarizedExperiment::colData(pb)$group_id %in% tg
  pb_sub <- pb[, keep_idx]
  pb_sub$group_id <- droplevels(factor(SummarizedExperiment::colData(pb_sub)$group_id))
  
  all_results <- list()
  cluster_names <- names(SummarizedExperiment::assays(pb_sub))
  
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
  
  message("6/7: dream 분석 실행 중...")
  pb_sub_meta <- as.data.frame(SummarizedExperiment::colData(pb_sub))
  pb_sub_meta <- normalize_sample_column(pb_sub_meta)
  pb_sub_meta$sample_id <- as.character(pb_sub_meta$sample_id)
  
  for (clust in cluster_names) {
    pb_clust <- SummarizedExperiment::assays(pb_sub)[[clust]]
    if (ncol(pb_clust) < 2) {
      message(sprintf("  클러스터 %s: 샘플 수 부족 (%d), 건너뜁니다.", clust, ncol(pb_clust)))
      next
    }
    
    sample_ids <- colnames(pb_clust)
    pb_clust_meta <- dplyr::left_join(
      data.frame(sample_id = sample_ids, stringsAsFactors = FALSE),
      sample_lookup,
      by = "sample_id"
    )
    pb_clust_meta <- pb_clust_meta[!is.na(pb_clust_meta$group_id), , drop = FALSE]
    if (nrow(pb_clust_meta) < 2) {
      message(sprintf("  클러스터 %s: 그룹 정보가 부족하여 건너뜁니다.", clust))
      next
    }
    keep_samples <- sample_ids[sample_ids %in% pb_clust_meta$sample_id]
    if (length(keep_samples) < 2) {
      message(sprintf("  클러스터 %s: sample_id 매칭 실패, 건너뜁니다.", clust))
      next
    }
    pb_clust <- pb_clust[, keep_samples, drop = FALSE]
    pb_clust_meta <- pb_clust_meta[match(keep_samples, pb_clust_meta$sample_id), , drop = FALSE]
    pb_clust_meta$group_id <- droplevels(factor(pb_clust_meta$group_id))
    pb_clust_meta$group <- pb_clust_meta$group_id
    
    group_counts <- table(pb_clust_meta$group)
    if (length(group_counts) < 2 || any(group_counts < min_samples_per_group)) {
      message(sprintf("  클러스터 %s: 그룹별 샘플 수 부족 (%s), 건너뜁니다.",
                      clust,
                      paste(names(group_counts), group_counts, sep=":", collapse=", ")))
      next
    }
    
    lib_sizes <- colSums(pb_clust)
    valid_samples <- lib_sizes > 0
    if (sum(valid_samples) < 2) {
      message(sprintf("  클러스터 %s: 유효 샘플 수 부족, 건너뜁니다.", clust))
      next
    }
    pb_clust <- round(pb_clust[, valid_samples, drop = FALSE])
    storage.mode(pb_clust) <- "integer"
    pb_clust_meta <- pb_clust_meta[valid_samples, , drop = FALSE]
    pb_clust_meta$group <- droplevels(factor(pb_clust_meta$group))
    
    patient_vec <- NULL
    if (!is.null(patient_col) && patient_col %in% colnames(pb_clust_meta)) {
      patient_vec <- pb_clust_meta[[patient_col]]
    } else if (sample_id %in% colnames(pb_clust_meta)) {
      patient_vec <- pb_clust_meta[[sample_id]]
    } else if ("sample_id" %in% colnames(pb_clust_meta)) {
      patient_vec <- pb_clust_meta[["sample_id"]]
    }
    if (is.null(patient_vec) || length(patient_vec) == 0) {
      patient_vec <- rownames(pb_clust_meta)
    }
    if (is.null(patient_vec) || length(patient_vec) == 0) {
      message(sprintf("  클러스터 %s: random effect 컬럼을 찾을 수 없습니다, 건너뜁니다.", clust))
      next
    }
    pb_clust_meta$patient_effect <- droplevels(factor(patient_vec))
    if (length(levels(pb_clust_meta$patient_effect)) < 2) {
      message(sprintf("  클러스터 %s: random effect 레벨 수 부족, 건너뜁니다.", clust))
      next
    }
    
    form <- stats::as.formula("~ 0 + group + (1|patient_effect)")
    vobj <- tryCatch(
      variancePartition::voomWithDreamWeights(
        pb_clust,
        form,
        data = as.data.frame(pb_clust_meta),
        quiet = TRUE
      ),
      error = function(e) {
        message(sprintf("  클러스터 %s: voomWithDreamWeights 실패 (%s)", clust, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(vobj)) next
    
    fit <- tryCatch(
      variancePartition::dream(
        vobj,
        form,
        data = as.data.frame(pb_clust_meta),
        quiet = TRUE
      ),
      error = function(e) {
        message(sprintf("  클러스터 %s: dream 실패 (%s)", clust, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) next
    
    design <- stats::model.matrix(~ 0 + group, data = as.data.frame(pb_clust_meta))
    contrast_fixed <- fix_contrast(contrast, colnames(design))
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_fixed, levels = colnames(design))
    
    fit_contr <- variancePartition::contrastCoefficients(fit, contrast_matrix)
    tt <- limma::topTable(fit_contr, coef = 1, number = Inf, sort.by = "none")
    tt$cluster_id <- clust
    if (!is.null(cluster_label_map)) {
      tt$cluster_label <- cluster_label_map[as.character(clust)]
      tt$cluster_label[is.na(tt$cluster_label)] <- as.character(clust)
    }
    all_results[[clust]] <- tt
  }
  
  if (length(all_results) == 0) {
    warning("dream 분석 결과가 없습니다.")
    return(data.frame())
  }
  
  message("7/7: 결과 정리 중...")
  res_df <- do.call(rbind, all_results)
  return(res_df)
}


