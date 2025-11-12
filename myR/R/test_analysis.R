# ============================================================================
# Differential Expression Analysis Functions
# ============================================================================
# Consolidated from test_working.R (priority 1), test.R (priority 2), 
# test_cursor.R (priority 3)
# Version names are explicitly included to avoid conflicts
# ============================================================================

#' Run MAST Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using MAST (Model-based Analysis
#' of Single-cell Transcriptomics). MAST uses a hurdle model that separately
#' models the probability of expression (logistic component) and the level of
#' expression (continuous component).
#'
#' @param sobj Seurat object
#' @param formula Character string or formula object specifying the model
#'   (e.g., "~ g3" or "~ g3 + batch")
#' @param min_cells_expr Minimum number of cells expressing a gene for it to
#'   be included in the analysis (default: 10)
#' @param n_cores Number of cores for parallel execution (default: 4)
#' @param lrt_variable Variable name for likelihood ratio test (e.g., "g3")
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item primerid: Gene identifier
#'     \item p_value_hurdle: P-value from hurdle model
#'     \item coef: Coefficient estimate
#'     \item ci.hi: Upper confidence interval
#'     \item ci.lo: Lower confidence interval
#'   }
#'
#' @note
#' This function operates on single-cell level data. For pseudobulked analysis,
#' consider using \code{runMUSCAT_v5} or applying pseudobulking before calling
#' this function.
#'
#' @export
runMAST_v1 <- function(sobj,
                    formula,
                    min_cells_expr = 10,
                    n_cores = 4,
                    lrt_variable = NULL) {
  
  if (!requireNamespace("MAST", quietly = TRUE)) stop("MAST 패키지가 필요합니다.")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("SCE 패키지가 필요합니다.")
  if (is.null(lrt_variable)) stop("'lrt_variable' 인자 (예: 'g3')를 지정해야 합니다.")
  
  if (is.character(formula)) {
    formula_obj <- as.formula(formula)
  } else if (inherits(formula, "formula")) {
    formula_obj <- formula
  } else {
    stop("'formula'는 문자열 또는 formula 객체여야 합니다.")
  }
  
  # --- 1. SCA 객체 생성 및 정규화 ---
  message("1/5: Seurat -> SingleCellAssay(SCA) 객체 변환 중...")
  
  # MAST requires SingleCellAssay, not SingleCellExperiment
  # Convert via SingleCellExperiment then to SCA using SceToSingleCellAssay
  sce <- Seurat::as.SingleCellExperiment(sobj)
  
  # Use MAST::SceToSingleCellAssay for conversion
  sca <- tryCatch({
    MAST::SceToSingleCellAssay(sce)
  }, error = function(e) {
    # Fallback: use FromMatrix if SceToSingleCellAssay fails
    warning("SceToSingleCellAssay failed, trying FromMatrix...")
    counts_mat <- SummarizedExperiment::assay(sce, "counts")
    if (is.null(counts_mat) || nrow(counts_mat) == 0) {
      assay_names <- SummarizedExperiment::assayNames(sce)
      if (length(assay_names) > 0) {
        counts_mat <- SummarizedExperiment::assay(sce, assay_names[1])
      } else {
        stop("No counts matrix found in Seurat object")
      }
    }
    if (inherits(counts_mat, "sparseMatrix")) {
      counts_mat <- as.matrix(counts_mat)
    }
    cdata <- SummarizedExperiment::colData(sce)
    fdata <- SummarizedExperiment::rowData(sce)
    if (nrow(fdata) == 0 || !"primerid" %in% colnames(fdata)) {
      fdata <- S4Vectors::DataFrame(primerid = rownames(counts_mat))
      rownames(fdata) <- rownames(counts_mat)
    }
    MAST::FromMatrix(exprsArray = counts_mat, cData = cdata, fData = fdata, check_sanity = FALSE)
  })
  
  # --- 2. 유전자 필터링 ---
  message(sprintf("2/5: 유전자 필터링 (min %d cells)...", min_cells_expr))
  # freq()는 발현 비율 (0~1)을 반환
  keep_genes <- (MAST::freq(sca) * ncol(sca)) >= min_cells_expr
  sca_filtered <- sca[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(sca)))
  
  # --- 3. 정규화 (MAST는 log2cpm 사용) ---
  message("3/5: Log2(CPM+1) 정규화 중...")
  SummarizedExperiment::assay(sca_filtered, "logcpm") <- MAST::cpm(sca_filtered, log = TRUE)
  
  # --- 4. zlm (Hurdle LMM) 실행 ---
  message(sprintf("4/5: MAST::zlm 실행 (Cores: %d). 시간이 오래 걸릴 수 있습니다...", n_cores))
  
  zfit <- MAST::zlm(formula_obj, 
                    sca = sca_filtered, 
                    method = "glmer", 
                    parallel = TRUE,
                    nCores = n_cores)
  
  # --- 5. LRT 결과 요약 ---
  message(sprintf("5/5: LRT 검정 수행 (변수: %s)...", lrt_variable))
  summary_res <- summary(zfit, doLRT = lrt_variable)
  summary_dt <- summary_res$datatable
  
  results_df <- merge(
      summary_dt[component == 'H', .(primerid, `Pr(>Chisq)`)], # Hurdle (logistic)
      summary_dt[component == 'logcpm', .(primerid, coef, ci.hi, ci.lo)], # Continuous
      by = 'primerid'
  )
  colnames(results_df)[2] <- "p_value_hurdle"
  results_df <- results_df[order(p_value_hurdle), ]
  
  message("MAST 분석 완료.")
  return(results_df)
}

#' @export
runMAST <- function(...) {
  .Deprecated("runMAST_v1", package = "myR")
  runMAST_v1(...)
}


#' Run NEBULA Analysis (v1)
#'
#' @description
#' Performs differential expression analysis using NEBULA (Negative Binomial
#' mixed-effects model). NEBULA accounts for patient/sample-level random effects
#' and is suitable for multi-level experimental designs.
#'
#' @param sobj Seurat object
#' @param layer Assay layer to use (default: "counts")
#' @param fixed_effects Character vector of fixed effect variables
#'   (e.g., c("target_col", "celltype_col"))
#' @param covar_effects Character vector of covariate variables
#'   (e.g., c("batch_col"))
#' @param patient_col Column name for patient/sample ID (default: "patient_col")
#' @param offset Column name for offset variable (default: "nCount_RNA")
#' @param min_count Minimum number of cells expressing a gene (default: 10)
#'
#' @return NEBULA result object from \code{nebula::nebula()}
#'
#' @note
#' This function operates on single-cell level data. Pseudobulking can be
#' applied before calling this function using \code{create_pseudobulk()} or
#' similar utilities.
#'
#' @export
runNEBULA_v1 <- function(sobj,
                                layer = "counts",
                                fixed_effects = c("target_col", "celltype_col"),
                                covar_effects = c("batch_col"),
                                patient_col = "patient_col",
                                offset = "nCount_RNA",
                                min_count = 10) {
  
  # --- 0. 데이터 추출 ---
  meta <- sobj@meta.data
  counts <- GetAssayData(sobj, layer = layer) # dgCMatrix (희소 행렬)
  
  # --- 1. 유전자 필터링 ---
  message(sprintf("1/6: 유전자 필터링 (min %d cells)...", min_count))
  keep_genes <- rowSums(counts > 0) >= min_count
  counts_filtered <- counts[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(counts)))
  
  # --- 2. NA 값 확인 및 제거 ---
  # [수정] NEBULA는 'covar_effects'도 고정 효과로 처리해야 함
  all_fixed_vars <- c(fixed_effects, covar_effects)
  
  # [수정] NA 체크할 모든 변수 목록
  vars_to_check <- c(all_fixed_vars, patient_col, offset)
  
  message("2/6: 모델 변수에서 NA 값 확인 중...")
  message(paste("... 확인 대상:", paste(vars_to_check, collapse = ", ")))
  
  keep_cells_idx <- complete.cases(meta[, vars_to_check])
  n_removed <- sum(!keep_cells_idx)
  
  if (n_removed > 0) {
    message(sprintf("... NA 값으로 인해 %d 개의 세포를 제거합니다.", n_removed))
  }
  
  # NA가 없는 '깨끗한' 데이터 생성
  meta_clean <- meta[keep_cells_idx, ]
  counts_clean <- counts_filtered[, keep_cells_idx]
  
  message(sprintf("... 최종 분석 대상 세포: %d 개", nrow(meta_clean)))
  
  # --- 3. 디자인 행렬 및 벡터 생성 ---
  message("3/6: 디자인 행렬 및 벡터 생성 중...")
  
  # [수정] 범주형 변수와 숫자형 변수 분리
  # 'offset'은 숫자로 유지해야 함!
  factor_vars <- c(all_fixed_vars, patient_col)
  
  # [수정] lapply + as.factor 적용 (오타 수정 as.factors -> as.factor)
  meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], as.factor)
  
  # [수정] 'offset'은 numeric으로 강제 변환 (안전장치)
  meta_clean[[offset]] <- as.numeric(meta_clean[[offset]])
  if (any(is.na(meta_clean[[offset]]))) {
      stop(sprintf("'%s' 컬럼에 NA가 아닌 숫자형 값만 있어야 합니다.", offset))
  }
  
  # [수정] 'formula' 인수를 없애고, 변수명으로부터 동적 생성
  formula_str <- paste("~", paste(all_fixed_vars, collapse = " + "))
  message(sprintf("... 사용할 고정 효과 포뮬러: %s", formula_str))
  
  design_matrix <- model.matrix(as.formula(formula_str), data = meta_clean)
  
  # id 및 offset 벡터 추출
  id_vector <- meta_clean[[patient_col]]
  offset_vector <- meta_clean[[offset]]
  
  # --- 4. group_cell()로 데이터 정렬 ---
  message(sprintf("4/6: NEBULA를 위해 id (%s) 기준으로 데이터 정렬 중...", patient_col))
  data_grouped <- nebula::group_cell(
    count = counts_clean,
    id = id_vector,
    pred = design_matrix,
    offset = offset_vector
  )
  
  # --- 5. NEBULA 실행 ---
  message("5/6: NEBULA 실행 중 (NBLMM)...")
  # Add stability options to prevent convergence issues
  re_nebula <- tryCatch({
    nebula::nebula(
      count = data_grouped$count,
      id = data_grouped$id,
      pred = data_grouped$pred,
      offset = data_grouped$offset,
      model = "NBLMM", # lme4::glmer.nb와 유사한 Lognormal random effect 권장
      method = "HL" # Use HL (Hessian-Laplace) method for stability
    )
  }, error = function(e) {
    # If HL fails, try with fewer genes or different method
    warning("NEBULA with HL method failed, trying with default settings...")
    tryCatch({
      nebula::nebula(
        count = data_grouped$count,
        id = data_grouped$id,
        pred = data_grouped$pred,
        offset = data_grouped$offset,
        model = "NBLMM"
      )
    }, error = function(e2) {
      stop(sprintf("NEBULA failed: %s. Try reducing number of genes or checking data quality.", 
                   conditionMessage(e2)))
    })
  })
  
  # --- 6. 결과 반환 ---
  message("6/6: 분석 완료.")
  
  # [수정] 함수는 결과를 'return' 해야 함
  return(re_nebula)
}

#' @export
runNEBULA <- function(...) {
  .Deprecated("runNEBULA_v1", package = "myR")
  runNEBULA_v1(...)
}


#' Run MUSCAT Analysis (v5 - Pseudobulking)
#'
#' @description
#' Performs differential expression analysis using MUSCAT (Multi-sample
#' Multi-condition DE analysis). This function automatically performs
#' pseudobulking by aggregating cells within each cluster and sample, then
#' applies edgeR, DESeq2, or limma-voom for differential expression testing.
#'
#' @param sobj Seurat object
#' @param cluster_id Column name for cell type/cluster (default: "seurat_clusters")
#' @param sample_id Column name for sample ID (default: "hos_no")
#' @param group_id Column name for group/condition (default: "type")
#' @param batch_id Optional column name for batch variable
#' @param contrast Contrast string (e.g., "IS - SAH")
#' @param method DE method: "edgeR", "DESeq2", "limma-trend", or "limma-voom"
#'   (default: "edgeR")
#' @param pb_min_cells Minimum cells per pseudobulk sample (default: 3)
#' @param filter_genes Filtering method: "none", "genes", "both", or "edgeR"
#'   (default: "edgeR")
#' @param keep_clusters Optional vector of cluster IDs to keep
#' @param cluster_label_map Optional named vector mapping cluster IDs to labels
#'
#' @return Data frame with differential expression results per cluster
#'
#' @note
#' This function automatically performs pseudobulking via
#' \code{muscat::aggregateData()}. For single-cell level analysis, use
#' \code{runMAST_v1} or \code{runNEBULA_v1}.
#'
#' @export
runMUSCAT_v5 <- function(
  sobj,
  cluster_id = "seurat_clusters",
  sample_id  = "hos_no",
  group_id   = "type",
  batch_id   = NULL,                 # ex) "exp_batch"
  contrast   = NULL,                 # ex) "IS - SAH"
  method     = "edgeR",
  pb_min_cells = 3,
  filter_genes = c("none","genes","both","edgeR"),
  keep_clusters = NULL,
  cluster_label_map = NULL
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")
  filter_genes <- match.arg(filter_genes)

  # deps
  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # 1) Seurat -> SCE, prepSCE
  sce <- Seurat::as.SingleCellExperiment(sobj)
  sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

  # factor 보장
  sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
  sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
  sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
    sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
  }

  # 2) Pseudobulk
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  # (선택) 특정 클러스터만
  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  # 2-1) pb 메타 보강 (sample_id / group_id / batch)
  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))

  # sample_id 없으면 assay의 colnames로 복구
  if (!"sample_id" %in% names(pb_meta)) {
    first_assay <- names(SummarizedExperiment::assays(pb))[1]
    sid_guess <- colnames(SummarizedExperiment::assays(pb)[[first_assay]])
    if (is.null(sid_guess)) stop("pb에 sample_id가 없습니다.")
    pb_meta$sample_id <- sid_guess
    rownames(pb_meta) <- pb_meta$sample_id
    SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(pb_meta)
  }

  # sce에서 (sample_id -> group_id / batch) map
  sce_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  map_cols <- c("sample_id","group_id")
  if (!is.null(batch_id) && batch_id %in% names(sce_meta)) map_cols <- c(map_cols, batch_id)
  sce_map <- unique(sce_meta[, map_cols, drop=FALSE])

  # pb에 group_id / batch 보강
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

  # factor화
  pb$sample_id <- droplevels(factor(SummarizedExperiment::colData(pb)$sample_id))
  pb$group_id  <- droplevels(factor(SummarizedExperiment::colData(pb)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb))) {
    pb[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(pb)[[batch_id]]))
  }

  # 3) contrast 그룹만 자동 subset
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

  # **sce도 동일 기준으로 subset (resDS용 필수)**
  sce_sub <- sce[, sce$sample_id %in% SummarizedExperiment::colData(pb_sub)$sample_id &
                    sce$group_id  %in% tg]
  sce_sub$cluster_id <- droplevels(factor(sce_sub$cluster_id))
  sce_sub$sample_id  <- droplevels(factor(sce_sub$sample_id))
  sce_sub$group_id   <- droplevels(factor(sce_sub$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce_sub))) {
    sce_sub[[batch_id]] <- droplevels(factor(sce_sub[[batch_id]]))
  }

  # 4) design/contrast (batch는 'batch'로 복사해서 사용)
  pb_sub$group <- pb_sub$group_id
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(pb_sub))) {
    pb_sub$batch <- droplevels(factor(SummarizedExperiment::colData(pb_sub)[[batch_id]]))
    design <- stats::model.matrix(~ 0 + group + batch,
                                  data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
  } else {
    design <- stats::model.matrix(~ 0 + group,
                                  data = as.data.frame(SummarizedExperiment::colData(pb_sub)))
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
  contrast_fixed <- fix_contrast(contrast, colnames(design))
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_fixed, levels = design)

  # 5) pbDS
  res <- muscat::pbDS(
    pb_sub,
    design    = design,
    method    = method,
    contrast  = contrast_matrix,
    min_cells = pb_min_cells,
    filter    = filter_genes,
    verbose   = TRUE
  )

  # 6) 결과 평탄화: **sce_sub가 먼저, res가 다음**
  combined <- muscat::resDS(sce_sub, res)

  # cluster_id 정리 + 라벨
  if ("cluster" %in% names(combined) && !"cluster_id" %in% names(combined)) {
    combined$cluster_id <- combined$cluster
  }
  if (!"cluster_id" %in% names(combined)) stop("resDS 결과에 'cluster_id'가 없습니다.")
  if (!is.null(cluster_label_map)) {
    combined$cluster_label <- cluster_label_map[as.character(combined$cluster_id)]
    combined$cluster_label[is.na(combined$cluster_label)] <- as.character(combined$cluster_id)
  } else {
    combined$cluster_label <- as.character(combined$cluster_id)
  }

  return(combined)
}

#' @export
runMUSCAT <- function(...) {
  .Deprecated("runMUSCAT_v5", package = "myR")
  runMUSCAT_v5(...)
}

