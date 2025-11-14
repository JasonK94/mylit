# ============================================================================
# LDS (Limma-Dream-SVA) 파이프라인 모듈화
# ============================================================================

#' @title LDS 단계 1: 데이터 추출 및 검증
#' @description Seurat/DGEList/Matrix에서 count 데이터와 메타데이터 추출
#' @param sobj Seurat 객체, DGEList 객체, 또는 Raw Count Matrix
#' @param meta.data 메타데이터 (Matrix 입력 시 필수)
#' @param layer Seurat 객체의 count layer
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(counts_matrix, meta.data, n_cells, n_genes)
lds_01_extract_data <- function(sobj,
                                 meta.data = NULL,
                                 layer = "counts",
                                 save_intermediate = FALSE,
                                 output_dir = tempdir(),
                                 prefix = "lds") {
  message("1/7: 입력 데이터 처리 중...")
  
  if (inherits(sobj, "Seurat")) {
    counts_matrix <- Seurat::GetAssayData(sobj, layer = layer)
    if (is.null(meta.data)) {
      meta.data <- sobj@meta.data
    }
  } else if (inherits(sobj, "DGEList")) {
    counts_matrix <- sobj$counts
    if (is.null(meta.data)) {
      meta.data <- sobj$samples
    }
  } else if (is.matrix(sobj) || inherits(sobj, "dgCMatrix")) {
    counts_matrix <- sobj
    if (is.null(meta.data)) {
      stop("`sobj`가 Matrix일 경우, `meta.data`를 반드시 제공해야 합니다.")
    }
  } else {
    stop("`sobj`는 Seurat, DGEList, 또는 Matrix여야 합니다.")
  }
  
  if (ncol(counts_matrix) != nrow(meta.data)) {
    stop("Count matrix의 열 수(샘플)와 meta.data의 행 수(샘플)가 일치하지 않습니다.")
  }
  
  result <- list(
    counts_matrix = counts_matrix,
    meta.data = meta.data,
    n_cells = ncol(counts_matrix),
    n_genes = nrow(counts_matrix)
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_01_extract_data.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 1b: NA 필터링
#' @description 포뮬러 변수의 NA 값을 가진 셀 제거
#' @param counts_matrix Count matrix
#' @param meta.data 메타데이터
#' @param formula 포뮬러
#' @param remove_na NA 제거 여부
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(counts_matrix, meta.data, n_removed)
lds_01b_filter_na <- function(counts_matrix,
                               meta.data,
                               formula,
                               remove_na = TRUE,
                               save_intermediate = FALSE,
                               output_dir = tempdir(),
                               prefix = "lds") {
  if (!remove_na) {
    return(list(
      counts_matrix = counts_matrix,
      meta.data = meta.data,
      n_removed = 0
    ))
  }
  
  message("1b/7: NA 값 확인 및 필터링 중...")
  
  # 포뮬러에서 사용되는 모든 변수 추출
  all_vars <- all.vars(formula)
  # 메타데이터에 실제로 존재하는 변수만
  vars_to_check <- all_vars[all_vars %in% colnames(meta.data)]
  
  n_removed <- 0
  
  if (length(vars_to_check) > 0) {
    # NA가 있는 행 찾기
    na_rows <- !complete.cases(meta.data[, vars_to_check, drop = FALSE])
    n_na <- sum(na_rows)
    
    if (n_na > 0) {
      message(sprintf("... NA 값으로 인해 %d 개의 셀을 제거합니다.", n_na))
      message(sprintf("... 확인한 변수: %s", paste(vars_to_check, collapse = ", ")))
      
      # NA가 없는 인덱스
      keep_cells <- !na_rows
      
      # counts_matrix와 meta.data 모두 필터링
      counts_matrix <- counts_matrix[, keep_cells, drop = FALSE]
      meta.data <- meta.data[keep_cells, , drop = FALSE]
      
      n_removed <- n_na
      message(sprintf("... 필터링 후: %d 개의 셀", ncol(counts_matrix)))
    } else {
      message("... NA 값이 없습니다.")
    }
  } else {
    message("... 포뮬러 변수가 메타데이터에 없어 NA 체크를 스킵합니다.")
  }
  
  result <- list(
    counts_matrix = counts_matrix,
    meta.data = meta.data,
    n_removed = n_removed
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_01b_filter_na.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 2: 포뮬러 파싱
#' @description LMM 포뮬러에서 고정 효과와 임의 효과 분리
#' @param formula LMM 포뮬러
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(fixed_effects_formula, original_formula)
lds_02_parse_formula <- function(formula,
                                  save_intermediate = FALSE,
                                  output_dir = tempdir(),
                                  prefix = "lds") {
  message("2/7: 포뮬러 파싱 중...")
  
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("LMM 포뮬러 파싱을 위해 'lme4' 패키지가 필요합니다.")
  }
  
  # `formula`에서 고정 효과(Fixed Effects)만 추출
  fixed_effects_formula <- lme4::nobars(formula)
  
  if (is.null(fixed_effects_formula)) {
    # 예: ~ (1|emrid) + (1|set) 처럼 고정 효과가 아예 없는 경우
    fixed_effects_formula <- ~ 1 
  }
  
  result <- list(
    fixed_effects_formula = fixed_effects_formula,
    original_formula = formula
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_02_parse_formula.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 3: DGEList 생성, 필터링, 정규화
#' @description edgeR를 사용한 전처리
#' @param counts_matrix Count matrix
#' @param meta.data 메타데이터
#' @param fixed_effects_formula 고정 효과 포뮬러
#' @param min.count filterByExpr의 min.count
#' @param min.total.count filterByExpr의 min.total.count
#' @param min.prop filterByExpr의 min.prop
#' @param large.n filterByExpr의 large.n
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(dge, keep_genes, n_genes_before, n_genes_after)
lds_03_preprocess_dge <- function(counts_matrix,
                                   meta.data,
                                   fixed_effects_formula,
                                   min.count = 10,
                                   min.total.count = 15,
                                   min.prop = 0.1,
                                   large.n = 10,
                                   save_intermediate = FALSE,
                                   output_dir = tempdir(),
                                   prefix = "lds") {
  message("3/7: DGEList 생성 및 필터링 중...")
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("'edgeR' 패키지가 필요합니다.")
  }
  
  # `filterByExpr`를 위한 디자인 행렬 (고정 효과 기반)
  design_for_filter <- model.matrix(fixed_effects_formula, data = meta.data)
  
  dge <- edgeR::DGEList(counts_matrix, samples = meta.data)
  n_genes_before <- nrow(dge)
  
  # filterByExpr 옵션 적용
  keep_genes <- edgeR::filterByExpr(
    dge, 
    design = design_for_filter,
    min.count = min.count,
    min.total.count = min.total.count,
    min.prop = min.prop,
    large.n = large.n
  )
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  n_genes_after <- sum(keep_genes)
  
  if (n_genes_after == 0) {
    stop("모든 유전자가 필터링되었습니다. `filterByExpr` 조건을 확인하십시오.")
  }
  
  dge <- edgeR::calcNormFactors(dge)
  
  message(sprintf("... 유전자 필터링 완료: %d / %d 개 통과", 
                  n_genes_after, n_genes_before))
  
  result <- list(
    dge = dge,
    keep_genes = keep_genes,
    n_genes_before = n_genes_before,
    n_genes_after = n_genes_after
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_03_preprocess_dge.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 4: SVA 실행
#' @description Surrogate Variable Analysis 실행
#' @param dge DGEList 객체
#' @param meta.data 메타데이터
#' @param fixed_effects_formula 고정 효과 포뮬러
#' @param n_sv 사용할 SV 개수 (NULL이면 자동 결정)
#' @param sv_var_cutoff 잔차 분산 설명 비율
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(sva_obj, svs_final, n_sv_final, meta.data_with_sv)
lds_04_run_sva <- function(dge,
                           meta.data,
                           fixed_effects_formula,
                           n_sv = NULL,
                           sv_var_cutoff = 0.5,
                           save_intermediate = FALSE,
                           output_dir = tempdir(),
                           prefix = "lds") {
  message("4/7: SVA 실행 (숨겨진 변동성 탐색)...")
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("'limma' 패키지가 필요합니다.")
  }
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("'sva' 패키지가 필요합니다.")
  }
  
  mod_sva <- model.matrix(fixed_effects_formula, data = meta.data)
  mod0_sva <- model.matrix(~ 1, data = meta.data)
  
  # SVA는 voom 변환된 데이터에 실행
  v_sva <- limma::voom(dge, mod_sva, plot = FALSE)
  
  # 1) n.sv=NULL로 SVA를 실행하여 *최대* SV 개수 및 SV *모두* 찾기
  sva_obj <- sva::sva(v_sva$E, mod = mod_sva, mod0 = mod0_sva, n.sv = NULL)
  
  n_sv_max <- sva_obj$n.sv
  
  if (n_sv_max == 0) {
    message("... SVA가 유의미한 대리 변수(SV)를 찾지 못했습니다.")
    svs_final <- NULL
    n_sv_final <- 0
    meta.data_with_sv <- meta.data
  } else {
    message(sprintf("... SVA가 최대 %d개의 SV를 탐지했습니다.", n_sv_max))
    
    # 2) 사용할 SV 개수 결정
    if (!is.null(n_sv)) {
      # 사용자가 SV 개수 명시
      n_sv_final <- min(n_sv, n_sv_max)
    } else if (is.null(sv_var_cutoff)) {
      # sv_var_cutoff가 NULL이면 전체 SV 사용
      n_sv_final <- n_sv_max
      message(sprintf("... sv_var_cutoff=NULL이므로 전체 %d개 SV를 사용합니다.", n_sv_final))
    } else {
      # 잔차 분산 기반으로 SV 개수 자동 결정
      message(sprintf("... 잔차 분산의 %.0f%%를 설명하는 SV 개수 자동 탐색 중...", 
                      sv_var_cutoff * 100))
      
      # SVA가 사용한 것과 동일한 잔차(Residuals)를 수동 계산
      res_matrix <- t(resid(lm.fit(mod_sva, t(v_sva$E))))
      
      # 잔차 행렬의 SVD
      svd_res <- svd(res_matrix)
      
      # 각 SV가 설명하는 분산 비율
      percent_var_explained <- (svd_res$d^2) / sum(svd_res$d^2)
      cumulative_var <- cumsum(percent_var_explained)
      
      # Cutoff를 만족하는 최소 SV 개수 찾기
      n_sv_auto <- which(cumulative_var >= sv_var_cutoff)[1]
      
      if (is.na(n_sv_auto)) {
        n_sv_final <- n_sv_max
      } else {
        n_sv_final <- n_sv_auto
      }
      
      message(sprintf("... SV %d개가 잔차 분산의 %.1f%%를 설명합니다.", 
                      n_sv_final, cumulative_var[n_sv_final] * 100))
    }
    
    # 3) 최종 SV 추출 및 메타데이터에 추가
    if (n_sv_final > 0) {
      svs_final <- sva_obj$sv[, 1:n_sv_final, drop = FALSE]
      colnames(svs_final) <- paste0("SV", 1:n_sv_final)
      meta.data_with_sv <- cbind(meta.data, svs_final)
    } else {
      svs_final <- NULL
      meta.data_with_sv <- meta.data
    }
  }
  
  result <- list(
    sva_obj = sva_obj,
    svs_final = svs_final,
    n_sv_final = n_sv_final,
    meta.data_with_sv = meta.data_with_sv,
    v_sva = v_sva
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_04_run_sva.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 5: 최종 포뮬러 생성
#' @description 원본 포뮬러에 SV 추가
#' @param original_formula 원본 포뮬러
#' @param svs_final 사용된 SV 매트릭스
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(final_formula, final_formula_str)
lds_05_build_final_formula <- function(original_formula,
                                        svs_final,
                                        save_intermediate = FALSE,
                                        output_dir = tempdir(),
                                        prefix = "lds") {
  message("5/7: 최종 모델 포뮬러 생성 중...")
  
  original_formula_str <- paste(deparse(original_formula), collapse = "")
  
  if (!is.null(svs_final) && ncol(svs_final) > 0) {
    sv_str <- paste(colnames(svs_final), collapse = " + ")
    # 포뮬러에 SV 추가
    final_formula_str <- paste(original_formula_str, sv_str, sep = " + ")
  } else {
    final_formula_str <- original_formula_str
  }
  
  final_formula <- as.formula(final_formula_str)
  message(sprintf("... 최종 포뮬러: %s", final_formula_str))
  
  result <- list(
    final_formula = final_formula,
    final_formula_str = final_formula_str
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_05_build_final_formula.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 6: Limma-Dream 파이프라인 실행
#' @description voomWithDreamWeights, dream, eBayes 실행
#' @param dge DGEList 객체
#' @param final_formula 최종 포뮬러
#' @param meta.data 메타데이터 (SV 포함)
#' @param n_cores 병렬 처리 코어 수
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(fit_ebayes, v_dream, fit_dream)
lds_06_run_dream <- function(dge,
                              final_formula,
                              meta.data,
                              n_cores = 4,
                              save_intermediate = FALSE,
                              output_dir = tempdir(),
                              prefix = "lds") {
  message(sprintf("6/7: limma-dream 실행 (Core: %d개)...", n_cores))
  
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("'BiocParallel' 패키지가 필요합니다.")
  }
  if (!requireNamespace("variancePartition", quietly = TRUE)) {
    stop("'variancePartition' 패키지가 필요합니다.")
  }
  
  BPPARAM_SETUP <- BiocParallel::MulticoreParam(n_cores)
  
  # 1) 유효 가중치 계산 (voomWithDreamWeights)
  v_dream <- variancePartition::voomWithDreamWeights(
    dge, 
    final_formula, 
    meta.data, 
    BPPARAM = BPPARAM_SETUP
  )
  
  # 2) LMM 피팅 (dream)
  # dream은 contrast 없이 실행하고, topTable에서 coef를 지정하여 p-value 계산
  fit_dream <- variancePartition::dream(
    v_dream, 
    final_formula, 
    meta.data, 
    BPPARAM = BPPARAM_SETUP
  )
  
  # 3) Empirical Bayes 조정
  fit_ebayes <- limma::eBayes(fit_dream)
  
  # 4) 고정 효과 변수 추출
  fixed_vars_in_formula <- all.vars(final_formula)
  fixed_vars_in_formula <- fixed_vars_in_formula[
    !grepl("^SV\\d+$", fixed_vars_in_formula) & 
    !grepl("^\\(1\\|", fixed_vars_in_formula) &
    fixed_vars_in_formula != "(Intercept)"
  ]
  
  # 5) p-value 계산 (dream 결과에서 coefficients와 standard errors를 사용)
  # dream의 경우 eBayes 후에도 t-statistics가 NA일 수 있으므로 직접 계산
  if (length(fixed_vars_in_formula) > 0) {
    # p.value 매트릭스 초기화 (없는 경우)
    if (is.null(fit_ebayes$p.value)) {
      fit_ebayes$p.value <- matrix(NA, nrow = nrow(fit_ebayes$coefficients), 
                                    ncol = ncol(fit_ebayes$coefficients))
      colnames(fit_ebayes$p.value) <- colnames(fit_ebayes$coefficients)
      rownames(fit_ebayes$p.value) <- rownames(fit_ebayes$coefficients)
    }
    
    # df.total 확인 (벡터 또는 스칼라)
    df_total <- fit_ebayes$df.total
    n_genes <- nrow(fit_ebayes$coefficients)
    
    # t-statistics 매트릭스 초기화 (없는 경우)
    if (is.null(fit_ebayes$t)) {
      fit_ebayes$t <- matrix(NA, nrow = nrow(fit_ebayes$coefficients), 
                             ncol = ncol(fit_ebayes$coefficients))
      colnames(fit_ebayes$t) <- colnames(fit_ebayes$coefficients)
      rownames(fit_ebayes$t) <- rownames(fit_ebayes$coefficients)
    }
    
    # sigma (표준 오차) 확인
    sigma <- fit_ebayes$sigma
    if (length(sigma) == 1) {
      sigma <- rep(sigma, n_genes)
    } else if (length(sigma) != n_genes) {
      sigma <- rep(mean(sigma, na.rm = TRUE), n_genes)
    }
    
    # 각 고정 효과 변수에 대해 t-statistics와 p-value 계산
    for (coef_name in fixed_vars_in_formula) {
      if (coef_name %in% colnames(fit_ebayes$coefficients)) {
        coefs <- fit_ebayes$coefficients[, coef_name]
        
        # stdev.unscaled에서 표준 오차 계산
        if (!is.null(fit_ebayes$stdev.unscaled) && coef_name %in% colnames(fit_ebayes$stdev.unscaled)) {
          stdev_unscaled <- fit_ebayes$stdev.unscaled[, coef_name]
          # 표준 오차 = sigma * stdev.unscaled
          se <- sigma * stdev_unscaled
        } else {
          # stdev.unscaled가 없으면 NA로 설정
          se <- rep(NA, n_genes)
        }
        
        # t-statistics = coefficients / standard error
        t_stats <- coefs / se
        fit_ebayes$t[, coef_name] <- t_stats
        
        # df_total 처리
        if (length(df_total) == 1) {
          df_vec <- rep(df_total, n_genes)
        } else if (length(df_total) == n_genes) {
          df_vec <- df_total
        } else {
          df_vec <- rep(mean(df_total, na.rm = TRUE), n_genes)
        }
        
        # p-value 계산 (양측 검정)
        p_values <- rep(NA, n_genes)
        valid_idx <- !is.na(t_stats) & !is.na(df_vec) & is.finite(t_stats) & is.finite(df_vec) & !is.na(se) & se > 0
        if (sum(valid_idx) > 0) {
          p_values[valid_idx] <- 2 * pt(abs(t_stats[valid_idx]), df = df_vec[valid_idx], lower.tail = FALSE)
        }
        fit_ebayes$p.value[, coef_name] <- p_values
      }
    }
    message("... p-value 계산 완료 (coefficients 및 standard errors 기반)")
  } else {
    message("... dream 피팅 완료 (p-value는 topTable에서 계산됨)")
  }
  
  result <- list(
    fit_ebayes = fit_ebayes,
    v_dream = v_dream,
    fit_dream = fit_dream,
    contrast_applied = FALSE,  # contrast 없이 실행
    fixed_vars = fixed_vars_in_formula
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_06_run_dream.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 단계 7: SVA와 메타데이터 상관관계 분석
#' @description SVA와 메타데이터 변수 간 상관관계 계산
#' @param svs_final 사용된 SV 매트릭스
#' @param meta.data 메타데이터
#' @param plot_sva_correlation 상관관계 분석 수행 여부
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list(correlation_matrix, metadata_vars, sv_names) 또는 NULL
lds_07_analyze_sva_correlation <- function(svs_final,
                                            meta.data,
                                            plot_sva_correlation = TRUE,
                                            save_intermediate = FALSE,
                                            output_dir = tempdir(),
                                            prefix = "lds") {
  if (!plot_sva_correlation || is.null(svs_final) || ncol(svs_final) == 0) {
    return(NULL)
  }
  
  message("7b/7: SVA와 메타데이터 상관관계 분석 중...")
  
  # 메타데이터에서 숫자형/팩터형 변수만 선택
  meta_vars <- colnames(meta.data)
  # SV 변수 제외
  meta_vars <- meta_vars[!grepl("^SV\\d+$", meta_vars)]
  # 숫자형 또는 팩터형 변수만
  numeric_vars <- sapply(meta_vars, function(v) {
    is.numeric(meta.data[[v]]) || is.factor(meta.data[[v]])
  })
  meta_vars_numeric <- meta_vars[numeric_vars]
  
  if (length(meta_vars_numeric) == 0) {
    message("... 숫자형/팩터형 메타데이터 변수가 없어 상관관계 분석을 스킵합니다.")
    return(NULL)
  }
  
  # 상관관계 계산
  sv_cor_matrix <- matrix(NA, nrow = ncol(svs_final), ncol = length(meta_vars_numeric))
  rownames(sv_cor_matrix) <- colnames(svs_final)
  colnames(sv_cor_matrix) <- meta_vars_numeric
  
  for (i in 1:ncol(svs_final)) {
    for (j in seq_along(meta_vars_numeric)) {
      var_name <- meta_vars_numeric[j]
      var_data <- meta.data[[var_name]]
      
      # 팩터는 숫자로 변환
      if (is.factor(var_data)) {
        var_data <- as.numeric(var_data)
      }
      
      # NA 제거 후 상관관계 계산
      valid_idx <- !is.na(svs_final[, i]) & !is.na(var_data)
      if (sum(valid_idx) > 1) {
        sv_cor_matrix[i, j] <- cor(svs_final[valid_idx, i], var_data[valid_idx], method = "pearson")
      }
    }
  }
  
  result <- list(
    correlation_matrix = sv_cor_matrix,
    metadata_vars = meta_vars_numeric,
    sv_names = colnames(svs_final)
  )
  
  message(sprintf("... %d개의 메타데이터 변수와 상관관계 계산 완료", length(meta_vars_numeric)))
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_07_sva_correlation.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  return(result)
}

#' @title LDS 통합 파이프라인
#' @description 모든 단계를 순차적으로 실행하는 메인 함수
#' @param sobj Seurat 객체, DGEList 객체, 또는 Raw Count Matrix
#' @param formula 모델 포뮬러
#' @param meta.data 메타데이터
#' @param layer Seurat 객체의 count layer
#' @param n_sv 사용할 SV 개수
#' @param sv_var_cutoff 잔차 분산 설명 비율 (NULL이면 전체 SV 사용)
#' @param n_cores 병렬 처리 코어 수
#' @param remove_na NA 제거 여부
#' @param min.count filterByExpr의 min.count
#' @param min.total.count filterByExpr의 min.total.count
#' @param min.prop filterByExpr의 min.prop
#' @param large.n filterByExpr의 large.n
#' @param plot_sva_correlation SVA 상관관계 분석 여부
#' @param save_intermediate 중간 결과 저장 여부
#' @param output_dir 저장 디렉터리
#' @param prefix 파일명 접두사
#' @return list (fit, voom, sva_obj, svs_used, final_formula, dge, sv_correlation)
LDS <- function(sobj,
                formula,
                meta.data = NULL,
                layer = "counts",
                n_sv = NULL,
                sv_var_cutoff = 0.5,
                n_cores = max(1, parallel::detectCores() - 2),
                remove_na = TRUE,
                min.count = 10,
                min.total.count = 15,
                min.prop = 0.1,
                large.n = 10,
                plot_sva_correlation = TRUE,
                save_intermediate = FALSE,
                output_dir = tempdir(),
                prefix = "lds") {
  
  # 필수 패키지 확인
  required_packages <- c("BiocParallel", "edgeR", "limma", "sva", "variancePartition", "lme4", "qs")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("필수 패키지가 설치되지 않았습니다: ", paste(missing_packages, collapse = ", "))
  }
  
  # --- 단계 1: 데이터 추출 ---
  step1 <- lds_01_extract_data(
    sobj = sobj,
    meta.data = meta.data,
    layer = layer,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 1b: NA 필터링 ---
  step1b <- lds_01b_filter_na(
    counts_matrix = step1$counts_matrix,
    meta.data = step1$meta.data,
    formula = formula,
    remove_na = remove_na,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 2: 포뮬러 파싱 ---
  step2 <- lds_02_parse_formula(
    formula = formula,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 3: DGEList 전처리 ---
  step3 <- lds_03_preprocess_dge(
    counts_matrix = step1b$counts_matrix,
    meta.data = step1b$meta.data,
    fixed_effects_formula = step2$fixed_effects_formula,
    min.count = min.count,
    min.total.count = min.total.count,
    min.prop = min.prop,
    large.n = large.n,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 4: SVA 실행 ---
  step4 <- lds_04_run_sva(
    dge = step3$dge,
    meta.data = step1b$meta.data,
    fixed_effects_formula = step2$fixed_effects_formula,
    n_sv = n_sv,
    sv_var_cutoff = sv_var_cutoff,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 5: 최종 포뮬러 생성 ---
  step5 <- lds_05_build_final_formula(
    original_formula = step2$original_formula,
    svs_final = step4$svs_final,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 6: Dream 실행 ---
  step6 <- lds_06_run_dream(
    dge = step3$dge,
    final_formula = step5$final_formula,
    meta.data = step4$meta.data_with_sv,
    n_cores = n_cores,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 7: SVA 상관관계 분석 ---
  step7 <- lds_07_analyze_sva_correlation(
    svs_final = step4$svs_final,
    meta.data = step4$meta.data_with_sv,
    plot_sva_correlation = plot_sva_correlation,
    save_intermediate = save_intermediate,
    output_dir = output_dir,
    prefix = prefix
  )
  
  # --- 단계 8: SVA 상관관계 Heatmap 생성 ---
  step8 <- NULL
  if (plot_sva_correlation && !is.null(step4$svs_final) && ncol(step4$svs_final) > 0) {
    step8 <- lds_08_create_heatmaps(
      svs_final = step4$svs_final,
      meta.data = step4$meta.data_with_sv,
      output_dir = output_dir,
      prefix = prefix,
      method = "spearman",
      top_n = 10,
      p_threshold = 0.05,
      save_intermediate = save_intermediate
    )
  }
  
  message("7/7: 분석 완료.")
  
  # --- 최종 결과 반환 ---
  result_list <- list(
    fit = step6$fit_ebayes,
    voom = step6$v_dream,
    sva_obj = step4$sva_obj,
    svs_used = step4$svs_final,
    final_formula = step5$final_formula,
    dge = step3$dge,
    fixed_vars = step6$fixed_vars,
    contrast_applied = step6$contrast_applied
  )
  
  if (!is.null(step7)) {
    result_list$sv_correlation <- step7
  }
  
  if (!is.null(step8)) {
    result_list$heatmaps <- step8
  }
  
  return(result_list)
}

