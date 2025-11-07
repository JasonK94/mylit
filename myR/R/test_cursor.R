

#' @title NEBULA 파이프라인 함수
#' @description Seurat 객체를 받아 NA 처리, 데이터 정렬, NEBULA GLMM을 실행.
#'              Modeling에서, y~NB(mu,phi)이며, log(mu)=log(o)+BX+random effect 에서 "o"는 offset, random effect는 patient_col, X는 fixed effects와 covar_effects.
#'
#' @param sobj (Seurat) Seurat 객체
#' @param layer (character) Raw count가 저장된 layer (예: "counts")
#' @param fixed_effects (character vector) 고정 효과로 사용할 meta.data 컬럼명
#'                    (예: c("target_col", "celltype_col"))
#' @param covar_effects (character vector) 임의 효과처럼 보정하고 싶은 공변량.
#'                    NEBULA의 한계로 인해 '고정 효과'로 처리됩니다.
#'                    (예: c("batch_col", "percent.mt"))
#' @param patient_col (character) 환자 ID 컬럼. NEBULA의 'id' (유일한 임의 효과)로 사용됨.
#' @param offset (character) Offset으로 사용할 'numeric' 컬럼 (예: "nCount_RNA")
#' @param min_count (numeric) 유전자 필터링 기준 (최소 발현 세포 수)
#'
#' @return (nebula) NEBULA 실행 결과 객체
#'
#' @export 
runNEBULA <- function(sobj,
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
  re_nebula <- nebula::nebula(
    count = data_grouped$count,
    id = data_grouped$id,
    pred = data_grouped$pred,
    offset = data_grouped$offset,
    model = "NBLMM" # lme4::glmer.nb와 유사한 Lognormal random effect 권장
  )
  
  # --- 6. 결과 반환 ---
  message("6/6: 분석 완료.")
  
  # [수정] 함수는 결과를 'return' 해야 함
  return(re_nebula)
}

#' @title MAST 파이프라인 함수
#' @description lme4의 다중 임의 효과를 지원하는 zlm (Hurdle-LMM) 실행
#'
#' @param sobj (Seurat) Seurat 객체
#' @param formula (formula) lme4 문법의 포뮬러.
#'              (예: ~ g3 + (1|hos_no) + (1|GEM))
#' @param min_cells_expr (numeric) 유전자 필터링 기준 (최소 발현 세포 수)
#' @param n_cores (numeric) 병렬 처리 코어 수
#' @param lrt_variable (character) LRT 검정을 수행할 변수명 (예: "g3")
#'
#' @return (data.table) MAST 요약 결과
#'
#' @export
runMAST_v1 <- function(sobj,
                              formula,
                              min_cells_expr = 10,
                              n_cores = 4,
                              lrt_variable = NULL) {
  
  if (!requireNamespace("MAST", quietly = TRUE)) stop("MAST 패키지가 필요합니다.")
  if (is.null(lrt_variable)) stop("'lrt_variable' 인자 (예: 'g3')를 지정해야 합니다.")
  
  # --- 1. SCA 객체 생성 및 정규화 ---
  message("1/5: Seurat -> SingleCellAssay(SCA) 객체 변환 중...")
  sca <- MAST::FromSeurat(sobj, cData = TRUE, assay = "RNA", slot = "counts")
  
  # --- 2. 유전자 필터링 ---
  message(sprintf("2/5: 유전자 필터링 (min %d cells)...", min_cells_expr))
  # freq()는 발현 비율 (0~1)을 반환
  keep_genes <- (MAST::freq(sca) * ncol(sca)) >= min_cells_expr
  sca_filtered <- sca[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(sca)))

  # --- 3. 정규화 (MAST는 log2cpm 사용) ---
  message("3/5: Log2(CPM+1) 정규화 중...")
  # 이 단계에서 밀집 행렬이 생성될 수 있어 메모리 주의
  SummarizedExperiment::assay(sca_filtered, "logcpm") <- MAST::cpm(sca_filtered, log = TRUE)
  
  # --- 4. zlm (Hurdle LMM) 실행 ---
  message(sprintf("4/5: MAST::zlm 실행 (Cores: %d). 시간이 오래 걸릴 수 있습니다...", n_cores))
  
  # method="glmer"가 임의 효과를 처리
  zfit <- MAST::zlm(formula, 
                    sca = sca_filtered, 
                    method = "glmer", 
                    parallel = TRUE,
                    nCores = n_cores)
  
  # --- 5. LRT 결과 요약 ---
  message(sprintf("5/5: LRT 검정 수행 (변수: %s)...", lrt_variable))
  summary_res <- summary(zfit, doLRT = lrt_variable)
  summary_dt <- summary_res$datatable
  
  # p-value (Hurdle p-value) 기준 정렬
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

#' @title Muscat (Pseudo-bulking) 파이프라인 함수
#' @description 세포 유형별, 환자별로 count를 합산(Pseudo-bulk)하여 DE 분석
#'
#' @param sobj (Seurat) Seurat 객체. meta.data에 아래 3개 컬럼 포함 필수
#' @param cluster_id (character) 세포 유형 컬럼명 (예: "cell_type")
#' @param sample_id (character) 환자/샘플 ID 컬럼명 (예: "hos_no")
#' @param group_id (character) 비교할 그룹 컬럼명 (예: "g3")
#' @param formula_str (character) DE 분석에 사용할 포뮬러 (limma/edgeR/DESeq2)
#'                 (예: "~ 0 + group_id + batch_col")
#' @param contrast (character) 비교할 contrast (limma::makeContrasts 문법)
#'                 (예: "group_idLevelA-group_idLevelB")
#' @param method (character) DE 엔진: "dream", "edgeR", "DESeq2"
#'
#' @return (list) muscat::pbDS 결과 (세포 유형별 DE 테이블 리스트)
#'
#' @export
runMUSCAT_v1 <- function(sobj,
                                   cluster_id = "cell_type",
                                   sample_id = "hos_no",
                                   group_id = "g3",
                                   formula_str = "~ 0 + group + batch", # 'group'은 group_id를 의미
                                   contrast = NULL,
                                   method = "dream") {
  
  if (!requireNamespace("muscat", quietly = TRUE)) stop("muscat 패키지가 필요합니다.")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("SCE 패키지가 필요합니다.")
  if (is.null(contrast)) stop("'contrast' (예: 'groupA-groupB')를 지정해야 합니다.")

  # --- 1. SCE 객체 변환 ---
  message("1/4: Seurat -> SingleCellExperiment(SCE) 변환 중...")
  sce <- as.SingleCellExperiment(sobj)

  # --- 2. Pseudo-bulking 실행 ---
  message(sprintf("2/4: Pseudo-bulking (by %s, %s)...", cluster_id, sample_id))
  
  # muscat 포맷에 맞게 객체 준비
  # (sample_id, cluster_id, group_id를 특정 이름으로 지정)
  sce$sample_id <- sce[[sample_id]]
  sce$cluster_id <- sce[[cluster_id]]
  sce$group_id <- sce[[group_id]]
  
  pb <- muscat::aggregateData(sce,
                              assay = "counts",
                              by = c("cluster_id", "sample_id"),
                              sample_id = "sample_id",
                              cluster_id = "cluster_id")
  
  # --- 3. DE 분석 (pbDS) ---
  message(sprintf("3/4: muscat::pbDS 실행 (method: %s)...", method))
  
  # pbDS는 포뮬러의 'group_id'를 'group'으로 내부적으로 참조함
  # 따라서 formula_str에서 group_id 대신 'group' 사용
  final_formula_str <- gsub(group_id, "group", formula_str)
  final_contrast_str <- gsub(group_id, "group", contrast)

  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  design <- model.matrix(as.formula(final_formula_str), data = pb_meta)
  
  # limma::makeContrasts를 위한 환경 생성
  contrast_matrix <- limma::makeContrasts(contrasts = final_contrast_str, levels = design)

  # pbDS 실행
  res <- muscat::pbDS(pb,
                      formula = as.formula(final_formula_str),
                      design = design,
                      method = method,
                      contrast = contrast_matrix,
                      run_voom = (method == "dream"), # dream일 때만 voom
                      verbose = TRUE)
  
  # --- 4. 결과 반환 ---
  message("4/4: 분석 완료. 결과 테이블(res$table) 확인.")
  return(res)
}

runMAST <- function(sobj,
                              formula,
                              min_cells_expr = 10,
                              n_cores = 4,
                              lrt_variable = NULL) {
  
  if (!requireNamespace("MAST", quietly = TRUE)) stop("MAST 패키지가 필요합니다.")
  if (is.null(lrt_variable)) stop("'lrt_variable' 인자 (예: 'g3')를 지정해야 합니다.")
  # [신규] Formula 유연성 확보
  if (is.character(formula)) {
    formula_obj <- as.formula(formula)
  } else if (inherits(formula, "formula")) {
    formula_obj <- formula
  } else {
    stop("'formula'는 문자열 또는 formula 객체여야 합니다.")
  }
  # --- 1. SCA 객체 생성 및 정규화 ---
  message("1/5: Seurat -> SingleCellAssay(SCA) 객체 변환 중...")
  sca <- MAST::FromSeurat(sobj, cData = TRUE, assay = "RNA", slot = "counts")
  
  # --- 2. 유전자 필터링 ---
  message(sprintf("2/5: 유전자 필터링 (min %d cells)...", min_cells_expr))
  # freq()는 발현 비율 (0~1)을 반환
  keep_genes <- (MAST::freq(sca) * ncol(sca)) >= min_cells_expr
  sca_filtered <- sca[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(sca)))

  # --- 3. 정규화 (MAST는 log2cpm 사용) ---
  message("3/5: Log2(CPM+1) 정규화 중...")
  # 이 단계에서 밀집 행렬이 생성될 수 있어 메모리 주의
  SummarizedExperiment::assay(sca_filtered, "logcpm") <- MAST::cpm(sca_filtered, log = TRUE)
  
  # --- 4. zlm (Hurdle LMM) 실행 ---
  message(sprintf("4/5: MAST::zlm 실행 (Cores: %d). 시간이 오래 걸릴 수 있습니다...", n_cores))
  
  # method="glmer"가 임의 효과를 처리
  zfit <- MAST::zlm(formula_obj, 
                    sca = sca_filtered, 
                    method = "glmer", 
                    parallel = TRUE,
                    nCores = n_cores)
  
  # --- 5. LRT 결과 요약 ---
  message(sprintf("5/5: LRT 검정 수행 (변수: %s)...", lrt_variable))
  summary_res <- summary(zfit, doLRT = lrt_variable)
  summary_dt <- summary_res$datatable
  
  # p-value (Hurdle p-value) 기준 정렬
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

runMUSCAT <- function(sobj,
                                   cluster_id = "cell_type",
                                   sample_id = "hos_no",
                                   group_id = "g3",
                                   formula_str = "~ 0 + group + batch", # 'group'은 group_id를 의미
                                   contrast = NULL,
                                   method = "dream",
                                   cell_level_covars = NULL) {
  
  if (!requireNamespace("muscat", quietly = TRUE)) stop("muscat 패키지가 필요합니다.")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("SCE 패키지가 필요합니다.")
  if (is.null(contrast)) stop("'contrast' (예: 'groupA-groupB')를 지정해야 합니다.")
# --- 1. SCE 객체 변환 및 메타데이터 준비 ---
  message("1/5: SCE 변환 및 메타데이터 준비 중...")
  sce <- as.SingleCellExperiment(sobj)
  
  # muscat 포맷에 맞게 객체 준비
  sce$sample_id <- sce[[sample_id]]
  sce$cluster_id <- sce[[cluster_id]]
  sce$group_id <- sce[[group_id]]
  
  # [신규] 세포 수준 공변량 집계 (Patient-Cluster-Level)
  meta_df <- as.data.frame(SummarizedExperiment::colData(sce))
  
  if (!is.null(cell_level_covars)) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr 패키지가 필요합니다.")
    
    message(sprintf("2/5: 세포 수준 공변량 (%s) 집계 중...", 
                    paste(cell_level_covars, collapse = ", ")))
                    
    # sample_id와 cluster_id로 그룹화하여 평균 계산
    agg_covars <- meta_df %>%
      dplyr::group_by(sample_id, cluster_id) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(cell_level_covars), mean, .names = "{.col}"), .groups = "drop")
  } else {
    message("2/5: 세포 수준 공변량 없음, 건너뜀...")
  }

  # --- 3. Pseudo-bulking 실행 ---
  message("3/5: Pseudo-bulking (by cluster, sample)...")
  pb <- muscat::aggregateData(sce,
                              assay = "counts",
                              by = c("cluster_id", "sample_id"),
                              sample_id = "sample_id",
                              cluster_id = "cluster_id")
  
  # [신규] 집계된 공변량을 'pb' 객체의 colData에 병합
  if (!is.null(cell_level_covars)) {
    pb_coldata <- as.data.frame(SummarizedExperiment::colData(pb))
    
    # pb$sample_id와 pb$cluster_id를 기준으로 병합
    pb_coldata_merged <- pb_coldata %>%
      dplyr::left_join(agg_covars, by = c("sample_id", "cluster_id"))
    
    # Rowname 순서가 중요하므로, 원본 순서 보장
    rownames(pb_coldata_merged) <- rownames(pb_coldata)
    SummarizedExperiment::colData(pb) <- DataFrame(pb_coldata_merged)
  }

  # --- 4. DE 분석 (pbDS) ---
  message(sprintf("4/5: muscat::pbDS 실행 (method: %s)...", method))
  
  # [수정] formula 유연성 (string -> formula)
  if (is.character(formula_str)) {
    formula_str <- gsub(group_id, "group", formula_str)
    formula_obj <- as.formula(formula_str)
  } else {
    stop("formula_str은 문자열(character)로 제공해야 합니다.")
  }
  
  if (is.character(contrast)) {
    contrast <- gsub(group_id, "group", contrast)
  } else {
    stop("contrast는 문자열(character)로 제공해야 합니다.")
  }

  pb_meta <- as.data.frame(SummarizedExperiment::colData(pb))
  design <- model.matrix(formula_obj, data = pb_meta)
  contrast_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)

  res <- muscat::pbDS(pb,
                      formula = formula_obj,
                      design = design,
                      method = method,
                      contrast = contrast_matrix,
                      run_voom = (method == "dream"),
                      verbose = TRUE)
  
  # --- 5. 결과 반환 ---
  message("5/5: 분석 완료.")
  return(res)
}