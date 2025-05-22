#' Statistical Analysis for scRNA-seq Metadata
#'
#' This function performs comprehensive statistical analysis on scRNA-seq data
#' comparing groups across different categorical variables.
#' @import dplyr ggplot2
#' @param transpose whether to... transpose
#'
#' @return a simple form of scatter plot with regression fit. 
#' @export
scatter_smooth = function(sobj, feature, clinical_variable = "nih_change", transpose = FALSE) {
  
  # 1) Make a per-cell data.frame with the expression for 'feature'
  df <- data.frame(
    sample_no = sobj@meta.data$sample_no,
    FEATURE = as.numeric(FetchData(sobj, vars = feature)[, feature])
  )
  
  # 2) Compute average expression per patient
  df_avg <- df %>%
    group_by(sample_no) %>%
    summarise(avg_FEATURE = mean(FEATURE, na.rm = TRUE))
  
  # 3) Merge with your clinical variable
  meta_patient <- sobj@meta.data %>%
    select(sample_no, all_of(clinical_variable)) %>%
    distinct(sample_no, .keep_all = TRUE)
  
  df_merged <- left_join(df_avg, meta_patient, by = "sample_no") %>%
    mutate(
      nih_change_numeric = as.numeric(as.character(.data[[clinical_variable]]))
    )
  
  # 4) Branch logic for transpose:
  #    transpose = FALSE -> x = clinical_variable, y = avg_FEATURE
  #    transpose = TRUE  -> x = avg_FEATURE,         y = clinical_variable
  
  if (!transpose) {
    # Model: avg_FEATURE ~ nih_change_numeric
    model <- lm(avg_FEATURE ~ nih_change_numeric, data = df_merged)
    summary(model)
    
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    pval <- signif(summary(model)$coefficients[2, 4], 3)
    label_text <- paste0("y = ", intercept, " + ", slope, " * x\np = ", pval)
    
    # Plot: x = nih_change_numeric, y = avg_FEATURE
    p = ggplot(df_merged, aes(x = nih_change_numeric, y = avg_FEATURE)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      annotate(
        "text",
        x = min(df_merged$nih_change_numeric, na.rm = TRUE),
        y = max(df_merged$avg_FEATURE, na.rm = TRUE),
        label = label_text, hjust = 0, vjust = 1, size = 5
      ) +
      theme_bw() +
      xlab(clinical_variable) +
      ylab(paste0("Average ", feature, " Expression"))
    
  } else {
    # transpose = TRUE -> we flip the roles:
    # Model: nih_change_numeric ~ avg_FEATURE
    model <- lm(nih_change_numeric ~ avg_FEATURE, data = df_merged)
    summary(model)
    
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    pval <- signif(summary(model)$coefficients[2, 4], 3)
    label_text <- paste0("y = ", intercept, " + ", slope, " * x\np = ", pval)
    
    # Plot: x = avg_FEATURE, y = nih_change_numeric
    p = ggplot(df_merged, aes(x = avg_FEATURE, y = nih_change_numeric)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      annotate(
        "text",
        x = min(df_merged$avg_FEATURE, na.rm = TRUE),
        y = max(df_merged$nih_change_numeric, na.rm = TRUE),
        label = label_text, hjust = 0, vjust = 1, size = 5
      ) +
      theme_bw() +
      xlab(paste0("Average ", feature, " Expression")) +
      ylab(clinical_variable)
  }
  
  return(p)
}


#' 유전자 발현 데이터에 대한 선형 회귀 분석 수행
#'
#' 이 함수는 Seurat 객체, 유전자 목록, 샘플 그룹을 지정하는 메타데이터 열,
#' 그리고 수치형 결과 변수에 대한 메TA데이터 열을 입력으로 받습니다.
#' 각 샘플 그룹별로 각 유전자의 평균 발현량을 계산하고, 이를 결과 변수와 병합한 후
#' 각 유전자에 대해 결과 변수를 예측 변수로 하는 단순 선형 회귀를 수행합니다.
#'
#' @param sobj Seurat 객체.
#' @param genes 분석할 유전자 이름의 문자형 벡터.
#' @param sample_col 샘플 또는 그룹을 식별하는 메타데이터 열 이름을 지정하는 문자열입니다.
#'   기본값은 "sample"입니다.
#' @param numeric_predictor 선형 모델에서 예측 변수(독립 변수)로 사용될 수치형
#'   결과 변수에 대한 메타데이터 열 이름을 지정하는 문자열입니다.
#'   기본값은 "severity_score"입니다. 이 변수는 모델에서 Y ~ X 일 때 X 역할을 합니다.
#' @return 각 유전자에 대한 선형 모델 결과(절편, 기울기, p-값, R-제곱)를 포함하는 데이터 프레임.
#'   각 모델은 `gene_expression ~ numeric_predictor` 형태로 적합됩니다.
#'
#' @importFrom dplyr group_by summarise across all_of select distinct left_join
#' @importFrom rlang sym
#' @importFrom stats lm coef
#' @import Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # 예제 Seurat 객체 (pbmc_small) 및 데이터 생성
#' if (requireNamespace("SeuratData", quietly = TRUE) &&
#'     requireNamespace("Seurat", quietly = TRUE)) {
#'   suppressWarnings(SeuratData::InstallData("pbmc3k"))
#'   data("pbmc3k.final")
#'   sobj <- pbmc3k.final
#'
#'   # 임의의 'sample' 및 'severity_score' 메타데이터 추가
#'   set.seed(123)
#'   n_cells <- ncol(sobj)
#'   # 각 세포에 고유한 샘플 ID 부여 (실제로는 세포보다 샘플 수가 적을 것)
#'   # 이 예제에서는 각 세포가 고유한 '샘플'이라고 가정하여 그룹화의 의미를 보여줌
#'   # 실제 사용 시에는 환자 ID 또는 실험 배치 ID 등이 될 것임
#'   # 여기서는 설명을 위해 10개의 가상 샘플을 만듦
#'   sample_ids <- factor(paste0("sample_", rep(1:10, length.out = n_cells)))
#'   sobj$sample <- sample_ids
#'
#'   # 각 '샘플'에 대한 severity_score 생성 (샘플별로 하나의 값을 가져야 함)
#'   unique_samples <- unique(sobj$sample)
#'   severity_scores <- runif(length(unique_samples), 1, 10)
#'   names(severity_scores) <- unique_samples
#'   sobj$severity_score <- severity_scores[sobj$sample]
#'
#'   # 분석할 유전자 선택 (예: 상위 2개 변동 유전자)
#'   sobj <- FindVariableFeatures(sobj, verbose = FALSE)
#'   genes_to_analyze <- head(VariableFeatures(sobj), 2)
#'
#'   # 함수 실행
#'   results <- pseudobulk_linear_fit(sobj,
#'                                  genes = genes_to_analyze,
#'                                  sample_col = "sample",
#'                                  numeric_predictor = "severity_score")
#'   print(results)
#' }
#' }
pseudobulk_linear_fit_legacy <- function(sobj, genes, sample_col = "sample", numeric_predictor = "severity_score") {
  
  # 입력값 유효성 검사 (선택 사항이지만 권장)
  if (!inherits(sobj, "Seurat")) {
    stop("`sobj`는 Seurat 객체여야 합니다.")
  }
  if (!is.character(genes) || length(genes) == 0) {
    stop("`genes`는 하나 이상의 유전자 이름을 포함하는 문자형 벡터여야 합니다.")
  }
  if (!sample_col %in% colnames(sobj@meta.data)) {
    stop(paste0("`sample_col` '", sample_col, "'이 Seurat 객체의 메타데이터에 없습니다."))
  }
  if (!numeric_predictor %in% colnames(sobj@meta.data)) {
    stop(paste0("`numeric_predictor` '", numeric_predictor, "'이 Seurat 객체의 메타데이터에 없습니다."))
  }
  if (!is.numeric(sobj@meta.data[[numeric_predictor]])) {
    # 시도: factor를 numeric으로 변환 (주의해서 사용)
    if (is.factor(sobj@meta.data[[numeric_predictor]])) {
      warning(paste0("`numeric_predictor` '", numeric_predictor, "'이 factor형입니다. 수치형으로 변환을 시도합니다."))
      # Factor의 level을 숫자로 변환하면 문제가 생길 수 있으므로, character로 변환 후 numeric으로 변환
      original_levels <- levels(sobj@meta.data[[numeric_predictor]])
      converted_numeric <- suppressWarnings(as.numeric(as.character(sobj@meta.data[[numeric_predictor]])))
      if(any(is.na(converted_numeric)) && !any(is.na(sobj@meta.data[[numeric_predictor]]))) {
        stop(paste0("`numeric_predictor` '", numeric_predictor, "'를 수치형으로 변환하는 데 실패했습니다. Factor 레벨이 숫자가 아닐 수 있습니다: ", paste(original_levels, collapse=", ")))
      }
      sobj@meta.data[[numeric_predictor]] <- converted_numeric
    } else {
      stop(paste0("`numeric_predictor` '", numeric_predictor, "'이 메타데이터에서 수치형이 아닙니다."))
    }
  }
  
  
  # 지정된 유전자에 대한 발현 데이터 가져오기
  # gene_list는 더 이상 필요 없음. `genes`를 직접 사용.
  expr_matrix <- FetchData(sobj, vars = genes)
  
  # 발현 데이터에 샘플 식별자 열 추가
  expr_matrix[[sample_col]] <- sobj@meta.data[[sample_col]]
  
  # 샘플(환자)별 평균 유전자 발현량 계산
  # !!sym(sample_col)을 사용하여 프로그래밍 방식으로 열 이름을 지정
  avg_expr_by_sample <- expr_matrix %>%
    group_by(!!sym(sample_col)) %>%
    summarise(across(all_of(genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  # .groups = "drop" 추가: group_by 후 summarise의 기본 동작 변경에 따른 경고 방지 및 명시적 ungrouping
  
  
  # 메타데이터에서 샘플 식별자와 수치형 예측 변수(결과) 열만 선택하고 중복 제거
  meta_sample_info <- sobj@meta.data %>%
    select(all_of(c(sample_col, numeric_predictor))) %>%
    distinct()
  
  # 주석 처리된 'type unifying' 코드는 현재 로직에서 불필요해 보이므로 제거합니다.
  # 만약 특정 상황에서 타입 변환이 필요하다면, 해당 부분에 명시적으로 추가하는 것이 좋습니다.
  
  # 평균 발현 데이터와 수치형 예측 변수(결과) 데이터 병합
  df_merged <- left_join(avg_expr_by_sample, meta_sample_info, by = sample_col)
  
  # 수치형 예측 변수(결과)를 명시적으로 수치형으로 변환 (이미 factor 등이 숫자형으로 변환되었을 수 있지만, 확인 차원)
  # 이 단계는 위에서 입력값 유효성 검사를 통해 이미 처리되었을 수 있습니다.
  # 하지만 left_join 과정에서 타입이 변경될 가능성에 대비해 한 번 더 확인하는 것도 나쁘지 않습니다.
  if (!is.numeric(df_merged[[numeric_predictor]])) {
    df_merged[[numeric_predictor]] <- as.numeric(as.character(df_merged[[numeric_predictor]]))
    if (any(is.na(df_merged[[numeric_predictor]]))) {
      warning(paste0("병합 후 '", numeric_predictor, "' 열에 NA가 생성되었습니다. 타입 변환을 확인하세요."))
    }
  }
  
  
  # 각 유전자에 대해 선형 모델(lm) 수행
  # 모델: gene_expression ~ numeric_predictor
  results_lm <- lapply(genes, function(current_gene) {
    # 포뮬러 생성: Y ~ X  (여기서 Y는 current_gene, X는 numeric_predictor)
    formula_str <- paste0("`", current_gene, "` ~ `", numeric_predictor, "`") # 백틱으로 특수문자 포함 열 이름 처리
    model <- lm(as.formula(formula_str), data = df_merged)
    summary_model <- summary(model)
    
    # 결과 추출 (계수 이름이 numeric_predictor와 일치하는지 확인)
    # 독립 변수가 하나이므로, 두 번째 계수가 numeric_predictor의 기울기임
    slope_coef_name <- numeric_predictor # 포뮬러에서 사용된 이름
    p_value_index <- 2 # 일반적으로 (Intercept) 다음이 예측 변수
    slope_value <- NA
    p_value_value <- NA
    
    if (slope_coef_name %in% rownames(summary_model$coefficients)) {
      slope_value <- coef(model)[slope_coef_name]
      p_value_value <- summary_model$coefficients[slope_coef_name, 4]
    } else if (nrow(summary_model$coefficients) >= 2) { # 만약 이름 매칭이 안되면 두번째 행을 가정
      slope_value <- coef(model)[2]
      p_value_value <- summary_model$coefficients[2, 4]
      warning(paste0("유전자 '", current_gene, "' 모델에서 예측변수 '", numeric_predictor, "'의 계수 이름을 직접 찾지 못했습니다. 두 번째 계수를 사용합니다."))
    }
    
    
    data.frame(
      gene = current_gene,
      intercept = coef(model)[1], # 절편
      slope = slope_value,         # 기울기
      p_value = p_value_value,     # p-값
      r_squared = summary_model$r.squared # R-제곱
    )
  })
  
  # 결과를 단일 데이터 프레임으로 결합
  results_lm_df <- do.call(rbind, results_lm)
  
  return(results_lm_df)
}

#' @title 유전자 발현 데이터에 대한 그룹별 또는 전체 선형 회귀 분석 수행 (최적화됨)
#' @description
#' 이 함수는 Seurat 객체, 유전자 목록, (선택적) 그룹 지정 메타데이터 열,
#' 샘플 식별 메타데이터 열, 그리고 수치형 예측 변수 메타데이터 열을 입력으로 받습니다.
#' 먼저 전체 샘플에 대한 유사벌크 평균 발현량을 계산하고 메타데이터와 병합합니다.
#' 그 후, 각 그룹별(또는 전체)로 각 유전자에 대해 수치형 예측 변수를 사용한 선형 회귀를 수행합니다.
#'
#' @param sobj Seurat 객체.
#' @param genes 분석할 유전자 이름의 문자형 벡터.
#' @param sample_col 샘플 또는 유사벌크 단위를 식별하는 메타데이터 열 이름. 기본값은 "sample".
#' @param numeric_predictor 선형 모델에서 예측 변수(독립 변수, X)로 사용될 수치형
#'   메타데이터 열 이름. 기본값은 "severity_score".
#' @param group_col (선택 사항) 그룹별 분석을 위한 메타데이터 열 이름.
#'   `NULL` (기본값)이면 전체 샘플에 대해 분석을 수행합니다.
#' @param p_adjust_method `stats::p.adjust`에 사용될 p-값 보정 방법. 기본값은 "BH".
#' @param min_samples_per_group 그룹별 또는 전체 분석 시 필요한 최소 고유 샘플 수. 기본값은 3.
#' @param min_distinct_predictor_values 회귀 분석을 위해 필요한 `numeric_predictor`의 최소 고유 값 수. 기본값은 2.
#'
#' @return 각 유전자(및 각 그룹)에 대한 선형 모델 결과 (절편, 기울기, 기울기 표준오차,
#'   모델에 사용된 샘플 수, p-값, 보정된 p-값, R-제곱)를 포함하는 데이터 프레임.
#'   `group_col`이 지정된 경우, 결과에 'group' 열이 추가됩니다.
#'
#' @importFrom dplyr group_by summarise across all_of select distinct left_join arrange bind_rows filter n_distinct
#' @importFrom rlang sym !! :=
#' @importFrom stats lm coef summary.lm p.adjust as.formula na.omit
#' @import Seurat
#' @importFrom tibble rownames_to_column
#'
#' @export
pseudobulk_linear_fit_legacy2 <- function(sobj,
                                          genes,
                                          sample_col = "sample",
                                          numeric_predictor = "severity_score",
                                          group_col = NULL,
                                          p_adjust_method = "BH",
                                          min_samples_per_group = 3,
                                          min_distinct_predictor_values = 2) {
  
  # --- 1. 입력값 유효성 검사 ---
  if (!inherits(sobj, "Seurat")) stop("`sobj`는 Seurat 객체여야 합니다.")
  if (!is.character(genes) || length(genes) == 0) stop("`genes`는 하나 이상의 유전자 이름을 포함하는 문자형 벡터여야 합니다.")
  if (!all(genes %in% rownames(sobj))) {
    missing_genes <- genes[!genes %in% rownames(sobj)]
    stop(paste0("일부 유전자가 Seurat 객체에 없습니다: ", paste(missing_genes, collapse=", ")))
  }
  meta_cols <- colnames(sobj@meta.data)
  if (!sample_col %in% meta_cols) stop(paste0("`sample_col` '", sample_col, "'이 메타데이터에 없습니다."))
  if (!numeric_predictor %in% meta_cols) stop(paste0("`numeric_predictor` '", numeric_predictor, "'이 메타데이터에 없습니다."))
  if (!is.null(group_col) && !group_col %in% meta_cols) stop(paste0("`group_col` '", group_col, "'이 메타데이터에 없습니다."))
  
  # numeric_predictor를 numeric으로 변환 시도 (원본 sobj@meta.data 수정)
  predictor_data_orig <- sobj@meta.data[[numeric_predictor]]
  if (!is.numeric(predictor_data_orig)) {
    warning_msg <- paste0("`numeric_predictor` '", numeric_predictor, "'이(가) ", class(predictor_data_orig)[1], "형입니다. 수치형으로 변환을 시도합니다.")
    original_NAs <- is.na(predictor_data_orig)
    if (is.factor(predictor_data_orig)) {
      converted_numeric <- suppressWarnings(as.numeric(as.character(predictor_data_orig)))
    } else if (is.character(predictor_data_orig)) {
      converted_numeric <- suppressWarnings(as.numeric(predictor_data_orig))
    } else {
      stop(paste0("`numeric_predictor` '", numeric_predictor, "'의 타입(", class(predictor_data_orig)[1], ")을 수치형으로 자동 변환할 수 없습니다."))
    }
    new_NAs_introduced <- any(is.na(converted_numeric) & !original_NAs)
    if (all(is.na(converted_numeric)) && !all(original_NAs)) stop(paste0("`numeric_predictor` 변환 실패: 모든 값이 NA가 되었습니다."))
    if (new_NAs_introduced) warning(paste0(warning_msg, " 일부 값이 NA로 변환되었습니다.")) else warning(warning_msg)
    sobj@meta.data[[numeric_predictor]] <- converted_numeric
    if (!is.numeric(sobj@meta.data[[numeric_predictor]])) stop(paste0("`numeric_predictor`를 수치형으로 최종 변환 실패."))
  }
  
  # --- 2. 사전 그룹별/전체 샘플 수 및 예측 변수 유효성 검사 (메타데이터 기반) ---
  meta_df_check <- sobj@meta.data %>%
    select(all_of(c(sample_col, numeric_predictor, group_col))) %>%
    filter(!is.na(.data[[numeric_predictor]])) # 예측 변수 NA인 샘플은 제외하고 카운트
  
  valid_groups <- c()
  if (is.null(group_col)) {
    n_unique_samples <- n_distinct(meta_df_check[[sample_col]])
    n_unique_preds <- n_distinct(meta_df_check[[numeric_predictor]])
    if (n_unique_samples >= min_samples_per_group && n_unique_preds >= min_distinct_predictor_values) {
      valid_groups <- "all_samples"
    } else {
      stop(paste0("전체 샘플에 대해 유효한 샘플 수(", n_unique_samples, "<", min_samples_per_group, ") 또는 예측 변수 고유값 수(", n_unique_preds, "<", min_distinct_predictor_values,")가 부족합니다."))
    }
  } else {
    meta_df_check_grouped <- meta_df_check %>%
      group_by(!!sym(group_col)) %>%
      summarise(
        n_unique_samples = n_distinct(.data[[sample_col]]),
        n_unique_preds = n_distinct(.data[[numeric_predictor]]),
        .groups = "drop"
      ) %>%
      filter(n_unique_samples >= min_samples_per_group, n_unique_preds >= min_distinct_predictor_values)
    
    valid_groups <- unique(as.character(meta_df_check_grouped[[group_col]]))
    all_defined_groups <- unique(as.character(sobj@meta.data[[group_col]]))
    invalid_groups <- setdiff(all_defined_groups[!is.na(all_defined_groups)], valid_groups)
    if (length(invalid_groups) > 0) {
      warning(paste0("다음 그룹들은 샘플 수 또는 예측 변수 다양성 부족으로 분석에서 제외됩니다: ", paste(invalid_groups, collapse=", ")))
    }
    if (length(valid_groups) == 0) {
      stop("모든 그룹이 분석을 위한 최소 샘플 수 또는 예측 변수 다양성 기준을 충족하지 못합니다.")
    }
  }
  
  # --- 3. 전체 샘플에 대한 유사벌크 및 메타데이터 준비 ---
  expr_data <- GetAssayData(sobj, assay = DefaultAssay(sobj), slot = "data")[genes, , drop = FALSE]
  expr_df_long <- as.data.frame(t(as.matrix(expr_data))) %>%
    rownames_to_column("cell_id_temp_") %>% # 실제 cell id를 사용하기보다 sobj@meta.data의 sample_col을 사용
    mutate(!!sym(sample_col) := sobj@meta.data[cell_id_temp_temp_, sample_col]) %>% #sobj@meta.data에서 sample_col 가져오기
    select(-cell_id_temp_)
  
  # 수정된 부분: sobj@meta.data에서 직접 sample_col을 가져와서 expr_df_long에 추가
  # 위 코드는 GetAssayData 순서와 sobj@meta.data 순서가 동일하다고 가정하는데, 명시적으로 cell id로 매칭하는 것이 더 안전.
  # 다만, GetAssayData의 colnames는 sobj의 colnames와 일치하므로, 다음과 같이 할 수 있음:
  cell_barcodes <- colnames(expr_data) # GetAssayData에서 얻은 cell barcode
  expr_df_transposed <- as.data.frame(t(as.matrix(expr_data))) # 행: cell, 열: gene
  expr_df_transposed[[sample_col]] <- sobj@meta.data[cell_barcodes, sample_col]
  
  
  avg_expr_all_samples <- expr_df_transposed %>%
    group_by(!!sym(sample_col)) %>%
    summarise(across(all_of(genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  cols_to_select <- c(sample_col, numeric_predictor)
  if (!is.null(group_col)) cols_to_select <- c(cols_to_select, group_col)
  
  meta_data_full <- sobj@meta.data %>%
    select(all_of(cols_to_select)) %>%
    distinct()
  
  df_merged_full <- left_join(avg_expr_all_samples, meta_data_full, by = sample_col)
  
  # numeric_predictor가 최종적으로 numeric인지 확인
  if (!is.numeric(df_merged_full[[numeric_predictor]])) {
    df_merged_full[[numeric_predictor]] <- as.numeric(as.character(df_merged_full[[numeric_predictor]]))
    # 여기서 NA 생성 경고는 이미 위에서 처리했으므로 생략, 단 최종 확인은 필요
    if (all(is.na(df_merged_full[[numeric_predictor]]))) {
      stop("최종 병합 데이터에서 numeric_predictor가 모두 NA가 되었습니다.")
    }
  }
  
  
  # --- 4. 내부 함수: 단일 그룹 (또는 전체)에 대한 분석 수행 (이미 병합된 데이터 사용) ---
  .run_lm_on_merged_data <- function(current_df_merged, current_genes, current_numeric_predictor) {
    if (nrow(current_df_merged) == 0) return(NULL)
    
    results_list <- lapply(current_genes, function(gene) {
      formula_str <- paste0("`", gene, "` ~ `", current_numeric_predictor, "`")
      
      # NA 값 처리 및 유효 데이터 확인
      model_data <- current_df_merged[, c(gene, current_numeric_predictor), drop = FALSE]
      model_data_complete <- stats::na.omit(model_data)
      
      n_samples_for_model <- nrow(model_data_complete)
      n_distinct_preds_in_model_data <- length(unique(model_data_complete[[current_numeric_predictor]]))
      
      if (n_samples_for_model < min_samples_per_group || n_distinct_preds_in_model_data < min_distinct_predictor_values) {
        warning(paste0("유전자 '", gene, "'에 대해 유효 데이터 부족 (샘플수:", n_samples_for_model, ", 예측변수 고유값:", n_distinct_preds_in_model_data, ")으로 모델을 적합할 수 없습니다."))
        return(data.frame(
          gene = gene, intercept = NA_real_, slope = NA_real_, slope_se = NA_real_,
          n_samples_in_model = n_samples_for_model, p_value = NA_real_, r_squared = NA_real_
        ))
      }
      
      model <- lm(as.formula(formula_str), data = model_data_complete)
      summary_model <- summary(model)
      coefs <- coef(summary_model)
      
      intercept_val <- NA_real_; slope_val <- NA_real_; slope_se_val <- NA_real_; p_val <- NA_real_
      if ("(Intercept)" %in% rownames(coefs)) intercept_val <- coefs["(Intercept)", "Estimate"]
      
      predictor_in_model <- names(coef(model))[2]
      if (!is.na(predictor_in_model) && predictor_in_model %in% rownames(coefs)) {
        slope_val <- coefs[predictor_in_model, "Estimate"]
        slope_se_val <- coefs[predictor_in_model, "Std. Error"]
        p_val <- coefs[predictor_in_model, "Pr(>|t|)"]
      } else if (nrow(coefs) >= 2) {
        warning(paste0("유전자 '", gene, "' 모델에서 예측변수 '", current_numeric_predictor, 
                       "'의 계수 이름을 명시적으로 찾지 못했습니다. 모델의 두 번째 계수를 사용합니다. (실제 모델 계수명: '", rownames(coefs)[2], "')"))
        slope_val <- coefs[2, "Estimate"]; slope_se_val <- coefs[2, "Std. Error"]; p_val <- coefs[2, "Pr(>|t|)"]
      }
      
      data.frame(
        gene = gene, intercept = intercept_val, slope = slope_val, slope_se = slope_se_val,
        n_samples_in_model = n_samples_for_model, p_value = p_val, r_squared = summary_model$r.squared
      )
    })
    bind_rows(results_list)
  }
  
  # --- 5. 분석 실행 (그룹별 또는 전체) ---
  final_results_list <- list()
  for (grp_val in valid_groups) {
    current_data_subset <- df_merged_full
    if (!is.null(group_col) && grp_val != "all_samples") { # "all_samples"는 group_col이 NULL일 때의 플레이스홀더
      current_data_subset <- df_merged_full %>% filter(!!sym(group_col) == grp_val)
    }
    
    # 각 그룹/전체에 대해 numeric_predictor 값이 충분히 다양한지 최종 확인 (df_merged_full에서 이미 필터링 되었을 수 있지만 안전장치)
    if(length(unique(na.omit(current_data_subset[[numeric_predictor]]))) < min_distinct_predictor_values ||
       nrow(distinct(na.omit(current_data_subset[, c(sample_col, numeric_predictor)]), .data[[sample_col]])) < min_samples_per_group ) {
      warning(paste0(ifelse(is.null(group_col), "전체 데이터셋", paste0("그룹 '", grp_val, "'")),
                     "에서 최종 분석 데이터의 예측변수 다양성 또는 샘플 수가 부족하여 건너뜁니다."))
      next # 다음 그룹으로
    }
    
    group_results_df <- .run_lm_on_merged_data(current_data_subset, genes, numeric_predictor)
    if (!is.null(group_results_df) && nrow(group_results_df) > 0) {
      group_results_df$group <- grp_val
      final_results_list[[grp_val]] <- group_results_df
    }
  }
  final_results_df <- bind_rows(final_results_list)
  
  # --- 6. 결과 정리 및 보정된 p-값 계산 ---
  if (is.null(final_results_df) || nrow(final_results_df) == 0) {
    warning("분석 결과 데이터 프레임이 비어있습니다.")
    return(data.frame(group=character(), gene=character(), intercept=numeric(), slope=numeric(), 
                      slope_se=numeric(), n_samples_in_model=integer(), p_value=numeric(), 
                      adj_p_value=numeric(), r_squared=numeric()))
  }
  
  valid_p_indices <- !is.na(final_results_df$p_value)
  if (any(valid_p_indices)) {
    final_results_df$adj_p_value <- NA_real_
    final_results_df$adj_p_value[valid_p_indices] <- stats::p.adjust(
      final_results_df$p_value[valid_p_indices], method = p_adjust_method
    )
  } else {
    final_results_df$adj_p_value <- NA_real_
  }
  
  cols_ordered <- c("group", "gene", "intercept", "slope", "slope_se", "n_samples_in_model", "p_value", "adj_p_value", "r_squared")
  final_results_df <- final_results_df[, intersect(cols_ordered, names(final_results_df)), drop = FALSE]
  
  return(final_results_df)
}


#' @title 그룹 간 선형 회귀 기울기 사후 비교
#' @description
#' `pseudobulk_linear_fit` 함수의 결과를 사용하여 각 유전자에 대해 그룹 간 기울기를 비교합니다.
#' 두 그룹의 경우 t-검정을 수행하고, 세 그룹 이상인 경우 모든 쌍에 대해 t-검정을 수행하고 p-값을 보정합니다.
#'
#' @param results_df `pseudobulk_linear_fit`에서 반환된 데이터 프레임.
#'   `group`, `gene`, `slope`, `slope_se`, `n_samples_in_model` 열을 포함해야 합니다.
#' @param gene_col `results_df`에서 유전자 이름을 포함하는 열의 이름. 기본값 "gene".
#' @param group_col `results_df`에서 그룹 식별자를 포함하는 열의 이름. 기본값 "group".
#' @param slope_col `results_df`에서 기울기 값을 포함하는 열의 이름. 기본값 "slope".
#' @param se_col `results_df`에서 기울기의 표준 오차를 포함하는 열의 이름. 기본값 "slope_se".
#' @param n_samples_col `results_df`에서 모델 피팅에 사용된 샘플 수를 포함하는 열의 이름. 기본값 "n_samples_in_model".
#' @param p_adjust_method 다중 비교를 위한 p-값 보정 방법. 기본값 "BH".
#'
#' @return 각 유전자 내 그룹 쌍 간의 기울기 비교 결과를 담은 데이터 프레임.
#'   (gene, group1, group2, slope1, slope2, se1, se2, n1, n2, t_statistic, df, p_value, adj_p_value).
#'
#' @importFrom dplyr group_by do filter arrange
#' @importFrom rlang sym !!
#' @importFrom utils combn
#' @importFrom stats pt
#'
#' @export
post_hoc_slope_comparison_legacy <- function(results_df,
                                             gene_col = "gene",
                                             group_col = "group",
                                             slope_col = "slope",
                                             se_col = "slope_se",
                                             n_samples_col = "n_samples_in_model",
                                             p_adjust_method = "BH") {
  
  if (!all(c(gene_col, group_col, slope_col, se_col, n_samples_col) %in% names(results_df))) {
    stop("필수 열 중 일부가 results_df에 없습니다. (gene, group, slope, slope_se, n_samples_in_model)")
  }
  
  # NA 값이 있는 행 제거 (기울기, SE, 샘플 수 중 하나라도 NA면 비교 불가)
  results_df_complete <- results_df %>%
    filter(!is.na(.data[[slope_col]]), !is.na(.data[[se_col]]), !is.na(.data[[n_samples_col]]))
  
  comparison_results <- results_df_complete %>%
    group_by(!!sym(gene_col)) %>%
    do({
      gene_data <- .
      if (nrow(gene_data) < 2) { # 비교할 그룹이 2개 미만이면 건너뜀
        return(data.frame())
      }
      
      pairs <- if (nrow(gene_data) == 2) { # 정확히 2개 그룹
        matrix(gene_data[[group_col]], ncol=1) # combn 형식을 맞추기 위해 list of vectors
      } else { # 3개 이상 그룹
        combn(gene_data[[group_col]], 2, simplify = FALSE) # 각 쌍을 list의 element로
      }
      
      pairwise_tests <- lapply(pairs, function(pair_groups) {
        # combn(simplify=FALSE)는 list of vectors, simplify=TRUE (default)는 matrix
        # nrow(gene_data) == 2 일 때의 처리를 위해 pair_groups 구조 통일
        if (is.matrix(pair_groups)) { # combn(simplify=TRUE) for >2 groups
          g1_name <- pair_groups[1,1] # Assuming matrix output from combn for >2 groups handling
          g2_name <- pair_groups[2,1]
        } else { # list output or special handling for 2 groups
          g1_name <- pair_groups[1]
          g2_name <- pair_groups[2]
        }
        
        
        d1 <- gene_data[gene_data[[group_col]] == g1_name, ]
        d2 <- gene_data[gene_data[[group_col]] == g2_name, ]
        
        b1 <- d1[[slope_col]]; se1 <- d1[[se_col]]; n1 <- d1[[n_samples_col]]
        b2 <- d2[[slope_col]]; se2 <- d2[[se_col]]; n2 <- d2[[n_samples_col]]
        
        # 각 회귀 모델의 잔차 자유도 (df = n - p, p=2 (절편,기울기))
        # n_samples_in_model이 2 미만이면 slope_se가 NA일 것이므로 위에서 필터링됨.
        # 그러나 n_samples_in_model이 정확히 2이면 df_reg = 0이 되어 Welch-Satterthwaite에서 0으로 나누기 오류 발생.
        # slope_se가 계산되려면 n_samples_in_model >= 3 이어야 함. 이 조건은 pseudobulk_linear_fit에서 관리.
        df_reg1 <- n1 - 2
        df_reg2 <- n2 - 2
        
        if (df_reg1 <= 0 || df_reg2 <= 0) { # 자유도가 0 이하면 t-분포 정의 안됨
          t_stat <- NA_real_
          df_ws <- NA_real_
          p_val <- NA_real_
        } else {
          t_stat <- (b1 - b2) / sqrt(se1^2 + se2^2)
          
          # Welch-Satterthwaite degrees of freedom
          numerator_ws <- (se1^2 + se2^2)^2
          denominator_ws <- (se1^4 / df_reg1) + (se2^4 / df_reg2)
          
          if (denominator_ws == 0) { # 분모가 0이면 df 계산 불가 (se가 모두 0인 극단적 경우)
            df_ws <- NA_real_
            p_val <- NA_real_
          } else {
            df_ws <- numerator_ws / denominator_ws
            p_val <- 2 * pt(abs(t_stat), df = df_ws, lower.tail = FALSE) # 양측 검정
          }
        }
        
        data.frame(
          group1 = g1_name, group2 = g2_name,
          slope1 = b1, slope2 = b2,
          se1 = se1, se2 = se2,
          n1 = n1, n2 = n2,
          t_statistic = t_stat, df = df_ws, p_value_raw = p_val,
          stringsAsFactors = FALSE
        )
      }) # End lapply(pairs, ...)
      
      # Pairwise t-test 핸들링 수정: combn의 결과가 matrix일 때와 list일 때 구분
      if (nrow(gene_data) > 2) { # 3개 이상 그룹의 경우 combn은 matrix를 반환 (simplify=TRUE 기본값)
        # combn은 열별로 쌍을 반환
        pairs_matrix <- combn(gene_data[[group_col]], 2) # simplify=TRUE
        
        pairwise_tests_list <- list()
        for(j in 1:ncol(pairs_matrix)){
          g1_name <- pairs_matrix[1,j]
          g2_name <- pairs_matrix[2,j]
          
          d1 <- gene_data[gene_data[[group_col]] == g1_name, ]
          d2 <- gene_data[gene_data[[group_col]] == g2_name, ]
          
          b1 <- d1[[slope_col]]; se1 <- d1[[se_col]]; n1 <- d1[[n_samples_col]]
          b2 <- d2[[slope_col]]; se2 <- d2[[se_col]]; n2 <- d2[[n_samples_col]]
          
          df_reg1 <- n1 - 2
          df_reg2 <- n2 - 2
          
          t_stat <- NA_real_; df_ws <- NA_real_; p_val <- NA_real_
          if (df_reg1 > 0 && df_reg2 > 0 && se1 > 0 && se2 > 0) { # SE가 0이면 분모 0 가능성
            t_stat <- (b1 - b2) / sqrt(se1^2 + se2^2)
            numerator_ws <- (se1^2 + se2^2)^2
            denominator_ws <- (se1^4 / df_reg1) + (se2^4 / df_reg2)
            if(denominator_ws > 0) { # SE가 모두 0이 아니어야 함
              df_ws <- numerator_ws / denominator_ws
              p_val <- 2 * pt(abs(t_stat), df = df_ws, lower.tail = FALSE)
            }
          }
          pairwise_tests_list[[j]] <- data.frame(
            group1 = g1_name, group2 = g2_name, slope1 = b1, slope2 = b2,
            se1 = se1, se2 = se2, n1 = n1, n2 = n2,
            t_statistic = t_stat, df = df_ws, p_value_raw = p_val,
            stringsAsFactors = FALSE
          )
        }
        do.call(rbind, pairwise_tests_list)
        
      } else if (nrow(gene_data) == 2) { # 정확히 2개 그룹
        g1_name <- gene_data[[group_col]][1]
        g2_name <- gene_data[[group_col]][2]
        
        d1 <- gene_data[gene_data[[group_col]] == g1_name, ]
        d2 <- gene_data[gene_data[[group_col]] == g2_name, ]
        
        b1 <- d1[[slope_col]]; se1 <- d1[[se_col]]; n1 <- d1[[n_samples_col]]
        b2 <- d2[[slope_col]]; se2 <- d2[[se_col]]; n2 <- d2[[n_samples_col]]
        
        df_reg1 <- n1 - 2
        df_reg2 <- n2 - 2
        
        t_stat <- NA_real_; df_ws <- NA_real_; p_val <- NA_real_
        if (df_reg1 > 0 && df_reg2 > 0 && se1 > 0 && se2 > 0) {
          t_stat <- (b1 - b2) / sqrt(se1^2 + se2^2)
          numerator_ws <- (se1^2 + se2^2)^2
          denominator_ws <- (se1^4 / df_reg1) + (se2^4 / df_reg2)
          if(denominator_ws > 0) {
            df_ws <- numerator_ws / denominator_ws
            p_val <- 2 * pt(abs(t_stat), df = df_ws, lower.tail = FALSE)
          }
        }
        data.frame(
          group1 = g1_name, group2 = g2_name, slope1 = b1, slope2 = b2,
          se1 = se1, se2 = se2, n1 = n1, n2 = n2,
          t_statistic = t_stat, df = df_ws, p_value_raw = p_val,
          stringsAsFactors = FALSE
        )
      } else { # 그룹이 1개 이하면 비교 불가 (위에서 이미 처리)
        data.frame()
      }
      
    }) %>% # End group_by %>% do
    filter(nrow(.) > 0) # 비어있는 데이터프레임 (비교 못한 유전자) 제거
  
  if (nrow(comparison_results) == 0) {
    warning("비교할 수 있는 그룹 쌍이 없거나, 모든 유전자의 그룹 수가 2개 미만입니다.")
    return(data.frame(gene=character(), group1=character(), group2=character(), 
                      slope1=numeric(), slope2=numeric(), se1=numeric(), se2=numeric(),
                      n1=integer(), n2=integer(), t_statistic=numeric(), df=numeric(),
                      p_value_raw=numeric(), adj_p_value=numeric()))
  }
  
  # 전체 pairwise p-값에 대해 보정 (유전자별로 보정할 수도 있으나, 전체 비교 수에 대해 한번에 하는 것이 일반적)
  # 또는 유전자별로 보정하려면, 위 do 블록 내에서 adj_p_value를 계산해야 함.
  # 여기서는 유전자별로 그룹 쌍이 여러 개 있을 수 있으므로, 각 유전자 내에서 보정.
  comparison_results <- comparison_results %>%
    group_by(!!sym(gene_col)) %>%
    mutate(adj_p_value = stats::p.adjust(.data$p_value_raw, method = p_adjust_method)) %>%
    ungroup() %>%
    arrange(!!sym(gene_col), .data$adj_p_value)
  
  return(comparison_results)
}

#' @title 그룹 간 선형 회귀 기울기 사후 비교
#' @description
#' `pseudobulk_linear_fit` 함수의 결과를 사용하여 각 유전자에 대해 그룹 간 기울기를 비교합니다.
#' 두 그룹의 경우 t-검정을 수행하고, 세 그룹 이상인 경우 모든 쌍에 대해 t-검정을 수행하고 p-값을 보정합니다.
#'
#' @param results_df `pseudobulk_linear_fit`에서 반환된 데이터 프레임.
#'   `group`, `gene`, `slope`, `slope_se`, `n_samples_in_model` 열을 포함해야 합니다.
#' @param gene_col `results_df`에서 유전자 이름을 포함하는 열의 이름. 기본값 "gene".
#' @param group_col `results_df`에서 그룹 식별자를 포함하는 열의 이름. 기본값 "group".
#' @param slope_col `results_df`에서 기울기 값을 포함하는 열의 이름. 기본값 "slope".
#' @param se_col `results_df`에서 기울기의 표준 오차를 포함하는 열의 이름. 기본값 "slope_se".
#' @param n_samples_col `results_df`에서 모델 피팅에 사용된 샘플 수를 포함하는 열의 이름. 기본값 "n_samples_in_model".
#' @param p_adjust_method 다중 비교를 위한 p-값 보정 방법. 기본값 "BH".
#'
#' @return 각 유전자 내 그룹 쌍 간의 기울기 비교 결과를 담은 데이터 프레임.
#'   (gene, group1, group2, slope1, slope2, se1, se2, n1, n2, t_statistic, df, p_value_raw, adj_p_value).
#'
#' @importFrom dplyr group_by do filter arrange mutate ungroup slice
#' @importFrom rlang sym !! .data
#' @importFrom utils combn
#' @importFrom stats pt p.adjust
#'
#' @export
post_hoc_slope_comparison_legacy2 <- function(results_df,
                                              gene_col = "gene",
                                              group_col = "group",
                                              slope_col = "slope",
                                              se_col = "slope_se",
                                              n_samples_col = "n_samples_in_model",
                                              p_adjust_method = "BH") {
  
  # --- 입력값 및 필수 열 확인 ---
  required_cols <- c(gene_col, group_col, slope_col, se_col, n_samples_col)
  if (!all(required_cols %in% names(results_df))) {
    missing_cols <- required_cols[!required_cols %in% names(results_df)]
    stop(paste0("필수 열 중 일부가 results_df에 없습니다: ", paste(missing_cols, collapse=", ")))
  }
  
  # --- 중복된 (유전자-그룹) 항목 처리 ---
  # results_df가 pseudobulk_linear_fit의 정상적인 출력이면 중복이 없어야 함.
  # 하지만, 만약을 위해 중복 제거 및 경고
  results_df_deduplicated <- results_df %>%
    group_by(!!sym(gene_col), !!sym(group_col)) %>%
    slice(1) %>% # 각 유전자-그룹 조합에 대해 첫 번째 행만 유지
    ungroup()
  
  if(nrow(results_df_deduplicated) < nrow(results_df)){
    warning("입력 `results_df`에 중복된 유전자-그룹 항목이 발견되었습니다. 각 조합의 첫 번째 항목을 사용합니다.")
  }
  
  # --- NA 값 및 유효하지 않은 샘플 수/SE 처리 ---
  # slope_se > 0 이어야 Welch-Satterthwaite에서 분모가 0이 되는 것을 방지.
  # n_samples_in_model >=3 이어야 df_reg > 0 이 됨 (pseudobulk_linear_fit에서 이미 처리)
  results_df_complete <- results_df_deduplicated %>%
    filter(
      !is.na(.data[[slope_col]]), 
      !is.na(.data[[se_col]]),
      .data[[se_col]] > 0, # 표준오차가 0보다 커야 함
      !is.na(.data[[n_samples_col]]),
      .data[[n_samples_col]] >= 3 # 최소 샘플 수 (df_reg > 0 보장)
    )
  
  if (nrow(results_df_complete) == 0) {
    warning("NA 값, 유효하지 않은 SE 또는 샘플 수 필터링 후 비교할 데이터가 없습니다.")
    return(data.frame(gene=character(), group1=character(), group2=character(), 
                      slope1=numeric(), slope2=numeric(), se1=numeric(), se2=numeric(),
                      n1=integer(), n2=integer(), t_statistic=numeric(), df=numeric(),
                      p_value_raw=numeric(), adj_p_value=numeric()))
  }
  
  # --- 유전자별 그룹 간 비교 수행 ---
  comparison_results_list <- results_df_complete %>%
    group_by(!!sym(gene_col)) %>%
    do({
      gene_data_for_current_gene <- .
      
      if (nrow(gene_data_for_current_gene) < 2) {
        # 이 유전자에 대해 비교할 그룹이 2개 미만이면 빈 데이터프레임 반환
        return(data.frame()) 
      }
      
      # combn은 그룹 이름 벡터를 받아 조합 매트릭스를 반환 (열별로 쌍)
      group_combinations_matrix <- combn(gene_data_for_current_gene[[group_col]], 2)
      
      individual_pair_results_list_for_gene <- lapply(1:ncol(group_combinations_matrix), function(col_idx) {
        g1_name <- group_combinations_matrix[1, col_idx]
        g2_name <- group_combinations_matrix[2, col_idx]
        
        # 각 그룹에 대한 데이터 추출 (이제 반드시 한 행이어야 함)
        d1 <- gene_data_for_current_gene[gene_data_for_current_gene[[group_col]] == g1_name, ]
        d2 <- gene_data_for_current_gene[gene_data_for_current_gene[[group_col]] == g2_name, ]
        
        b1 <- d1[[slope_col]]; se1 <- d1[[se_col]]; n1 <- d1[[n_samples_col]]
        b2 <- d2[[slope_col]]; se2 <- d2[[se_col]]; n2 <- d2[[n_samples_col]]
        
        df_reg1 <- n1 - 2 # 각 회귀 모델의 잔차 자유도
        df_reg2 <- n2 - 2
        
        t_stat <- NA_real_; df_ws <- NA_real_; p_val <- NA_real_
        
        # df_reg1, df_reg2 > 0은 results_df_complete 필터링에서 n_samples_col >= 3으로 보장됨
        # se1, se2 > 0도 위에서 필터링됨
        t_stat <- (b1 - b2) / sqrt(se1^2 + se2^2)
        
        numerator_ws <- (se1^2 + se2^2)^2
        denominator_ws <- (se1^4 / df_reg1) + (se2^4 / df_reg2) # 이 값은 se > 0 이므로 0보다 큼
        
        df_ws <- numerator_ws / denominator_ws
        p_val <- 2 * stats::pt(abs(t_stat), df = df_ws, lower.tail = FALSE) # 양측 검정
        
        data.frame(
          # gene 열은 나중에 추가 (group_by 변수를 결과에 유지하는 옵션 사용)
          group1 = g1_name, group2 = g2_name,
          slope1 = b1, slope2 = b2,
          se1 = se1, se2 = se2,
          n1 = n1, n2 = n2,
          t_statistic = t_stat, df = df_ws, p_value_raw = p_val,
          stringsAsFactors = FALSE
        )
      }) # End lapply over columns of group_combinations_matrix
      
      # 현재 유전자에 대한 모든 쌍의 결과를 rbind
      if (length(individual_pair_results_list_for_gene) > 0) {
        do.call(rbind, individual_pair_results_list_for_gene)
      } else {
        data.frame() # 쌍이 없는 경우 (이론상 위 nrow < 2 에서 걸러짐)
      }
    }) %>% # End group_by %>% do
    ungroup() %>% # dplyr::do는 자동으로 ungroup하지 않으므로 명시적 ungroup
    filter(nrow(.) > 0) # 혹시 빈 데이터프레임이 반환된 경우 제거 (실제로는 do 내부에서 처리)
  
  if (nrow(comparison_results_list) == 0) {
    warning("최종적으로 비교 결과가 생성되지 않았습니다.")
    return(data.frame(gene=character(), group1=character(), group2=character(), 
                      slope1=numeric(), slope2=numeric(), se1=numeric(), se2=numeric(),
                      n1=integer(), n2=integer(), t_statistic=numeric(), df=numeric(),
                      p_value_raw=numeric(), adj_p_value=numeric()))
  }
  
  # p-값 보정 (각 유전자 내에서 수행된 모든 pairwise 비교에 대해)
  final_comparison_results <- comparison_results_list %>%
    group_by(!!sym(gene_col)) %>%
    mutate(adj_p_value = stats::p.adjust(.data$p_value_raw, method = p_adjust_method)) %>%
    ungroup() %>%
    arrange(!!sym(gene_col), .data$adj_p_value, .data$p_value_raw)
  
  return(final_comparison_results)
}