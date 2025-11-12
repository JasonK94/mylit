#' Input Validation Utilities
#'
#' This module provides common validation functions used across the package.
#' It standardizes error messages and validation logic.
#'
#' @name validation
NULL

#' @param sobj Seurat object
#' @param patient Column name for patient identifier
#' @param treatment Column name for treatment/drug identifier
#' @param timepoint Column name for timepoint identifier
#' @param tissue Column name for tissue identifier
#' @param check_key Column name for the key to check
#' @param check_values Values to check for parity
#' @return A list containing:
#'   \itemize{
#'     \item `detailed`: The original metadata with added diagnostic columns.
#'     \item `diag_table`: A one-row-per-group summary table.
#'     \item `message`: A summary message.
#'     \item `paired_ids`: A vector of `SegmentDisplayName`s for samples that passed the checks.
#'   }
#' @examples
#' \dontrun{
#' diagnosis_parity(sobj, patient="EMRID", treatment="drug", timepoint="pre_post", tissue="ck_str", check_key="pre_post", check_values=c("pre", "post"))
#' }
#' @export
diagnosis_parity <- function(sobj,
                                      patient=NULL, treatment=NULL, timepoint=NULL, tissue=NULL,   # 문자열 컬럼명 허용
                                      check_key, check_values,
                                      id_col="SegmentDisplayName") {
  # 동적 진단 열 이름: "pre","post" -> "n_pre","n_post"
  safe_counts_names <- paste0("n_", make.names(check_values, unique = TRUE))
  
  make_unique_id <- function(meta, patient=NULL, treatment=NULL, timepoint=NULL, tissue=NULL, check_key=NULL) {
    # unique_id를 만들 때 쓸 후보
    cols <- c(patient, treatment, timepoint, tissue)
  
    # check_key가 들어있으면 제외
    cols <- cols[!cols %in% check_key]
  
    # NULL 제거
    cols <- cols[!sapply(cols, is.null)]
  
    # 조합
    if (length(cols) == 0) {
      meta$unique_id <- seq_len(nrow(meta))  # 전부 NULL이면 그냥 row별 고유 번호
    } else {
      meta <- meta %>%
        mutate(across(all_of(cols), ~ tidyr::replace_na(as.character(.x), ""))) %>%
        mutate(unique_id = do.call(paste0, .[cols]))
    }
    meta
  }
  
  meta <- sobj@meta.data %>%
    make_unique_id(patient=patient,
                 treatment=treatment,
                 timepoint=timepoint,
                 tissue=tissue,
                 check_key=check_key) %>%
    mutate(
      .val_chr  = as.character(.data[[check_key]])
    )

  # unique_id별 전체 개수
  totals <- meta %>% count(unique_id, name = "total_n")

  # 선택값 카운트 (없어도 0으로 생성되도록 complete + pivot_wider)
  counts_sel <- meta %>%
    mutate(.val_chr = factor(.val_chr, levels = check_values)) %>%
    count(unique_id, .val_chr) %>%
    complete(unique_id, .val_chr, fill = list(n = 0)) %>%
    pivot_wider(names_from = .val_chr, values_from = n, values_fill = 0)

  # pivot에서 특정 값 컬럼 자체가 안 생긴 경우 대비해서 강제 추가
  for (v in check_values) {
    if (!hasName(counts_sel, v)) counts_sel[[v]] <- 0L
  }

  # 읽기 쉬운 이름으로 통일: pre -> n_pre, post -> n_post ...
  counts_sel <- counts_sel %>%
    rename_with(~ paste0("n_", make.names(.x, unique = TRUE)),
                .cols = all_of(check_values))

  # 진단 테이블 계산 (rowwise 없이)
  diag_table <- totals %>%
    left_join(counts_sel, by = "unique_id") %>%
    mutate(across(c(total_n, all_of(safe_counts_names)),
                  ~ replace_na(.x, 0))) %>%
    mutate(
      # 선택값 합계
      n_selected = purrr::reduce(across(all_of(safe_counts_names)), `+`),
      # 모든 선택값이 최소 1개 이상인가?
      has_all    = if_all(all_of(safe_counts_names), ~ .x > 0),
      # 선택값들의 개수가 모두 동일한가? (pre==post==... 형태)
      parity_eq  = pmap_lgl(across(all_of(safe_counts_names)),
                            ~ { v <- c(...); length(unique(v)) == 1 }),
      # 기타 값 개수
      n_other    = total_n - n_selected,
      # 최종 통과 조건(원래 의도 유지)
      ok         = has_all & parity_eq & n_other == 0
    ) %>%
    dplyr::select(unique_id, total_n, all_of(safe_counts_names),
           n_selected, n_other, has_all, parity_eq, ok)

  # 원본 롱형 테이블에 진단 열 조인
  detailed <- meta %>% left_join(diag_table, by = "unique_id")

  # 실패 id 메시지
  no_pass_ids <- diag_table %>% filter(!ok) %>% pull(unique_id)
  msg <- if (length(no_pass_ids) == 0) {
    "✅ All groups passed parity & exclusivity checks."
  } else {
    sprintf("checkout %s", paste(no_pass_ids, collapse = ", "))
  }
  message(msg)
  message("View($detailed)")
  message("enter View(sobj$diag_table) to check which row has the problem")
  ok_id=detailed%>%
    filter(ok)%>%
    pull(.data[[id_col]])
  
  list(
    detailed    = detailed,    # 원본 행 + n_pre/n_post/... + ok 등
    diag_table  = diag_table,  # id별 한 줄 요약
    message     = msg,
    paired_ids = ok_id
  )
}

#' @title 분석 입력값 검증 헬퍼 함수
#' @description 통계 분석/모델링에 사용될 데이터, 메타데이터, 변수, 포뮬러, 방법론을 검증합니다.
#'
#' @param data 기본 데이터 객체 (data.frame, matrix, Seurat 등).
#' @param metadata 메타데이터 (data.frame). NULL일 경우 `data`에서 추출을 시도합니다.
#' @param target_vars 분석 대상 변수 (캐릭터 벡터).
#' @param control_vars 통제 변수 (캐릭터 벡터).
#' @param formula 분석에 사용될 포뮬러 (formula 또는 character).
#' @param method 분석 방법론 (character, 예: "lm", "wilcoxon", "random forest").
#' @param auto_correct_vars (logical) TRUE일 경우, 문자형(character) 변수를 팩터(factor)로 자동 변환합니다.
#' @param auto_correct_formula (logical) TRUE일 경우, `method`에 맞게 포뮬러 수정을 시도합니다 (예: ML에 랜덤 이펙트 제거).
#' @param na_action (character) 'report'(기본값)는 NA를 보고만 하고, 'remove_rows'는 `target_vars`, `control_vars`에 
#'                  지정된 변수들에 NA가 있는 행을 제거한 `metadata_clean`을 반환합니다.
#'
#' @return 검증이 완료되고 (필요시) 수정된 `metadata`, `formula` 및 (na_action='remove_rows'인 경우) 
#'         NA가 제거된 `metadata_clean`을 포함하는 리스트(list).
#' @examples
#' \dontrun{
#' # ... (이전 예제 코드와 동일) ...
#'
#' # 3. NA 제거 옵션 테스트
#' validated_na_removed <- validate_analysis_inputs(
#'   data = sim_data,
#'   target_vars = "target",
#'   control_vars = c("group", "age", "sex"), # 'age'에 NA가 있음
#'   formula = "target ~ group + age + sex",
#'   method = "lm",
#'   na_action = "remove_rows"
#' )
#'
#' # 원본 메타데이터 행 개수
#' # nrow(validated_na_removed$metadata) 
#' # NA 제거된 메타데이터 행 개수
#' # nrow(validated_na_removed$metadata_clean)
#' }
validate_analysis_inputs <- function(data, metadata = NULL,
                                     target_vars = NULL, control_vars = NULL,
                                     formula = NULL, method = NULL,
                                     auto_correct_vars = FALSE,
                                     auto_correct_formula = FALSE,
                                     na_action = c("report", "remove_rows")) {
  
  na_action <- match.arg(na_action) # 파라미터 기본값 매칭

  message("--- 시작: 입력값 검증 리포트 ---")

  # --- 0. 초기 설정 및 반환 객체 준비 ---
  validated_output <- list(
    # data = data, # [수정] 원본 데이터 반환 제거 (낭비)
    metadata = NULL, # 수정된 메타데이터가 여기에 할당됨
    formula = NULL,  # 수정된 포뮬러가 여기에 할당됨
    metadata_clean = NULL, # NA가 제거된 메타데이터가 여기에 할당됨
    report = list() # 검증 메시지를 저장할 리스트
  )
  
  is_seurat <- inherits(data, "Seurat")
  data_obs <- 0 # 관측치 개수 (행)
  
  # --- 1. Data 객체 검증 ---
  message("\n1. Data 객체:")
  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      message("  - ERROR: Seurat 객체를 다루려면 'SeuratObject' 패키지가 필요합니다.")
      return(invisible(NULL))
    }
    data_obs <- nrow(data@meta.data)
    message(paste0("  - Type: Seurat Object (", data_obs, " cells, ", nrow(data), " features)"))
  } else if (is.data.frame(data)) {
    data_obs <- nrow(data)
    message(paste0("  - Type: data.frame (", data_obs, " obs, ", ncol(data), " vars)"))
  } else if (is.matrix(data)) {
    data_obs <- nrow(data)
    message(paste0("  - Type: matrix (", data_obs, " rows, ", ncol(data), " cols)"))
    message("  - Warning: matrix의 경우, metadata를 명시적으로 제공해야 합니다.")
  } else {
    message(paste0("  - Type: ", class(data)[1], " (지원되지 않거나 확인 필요)"))
    try(data_obs <- nrow(data), silent = TRUE)
  }
  
  if (data_obs == 0) {
    message("  - ERROR: 데이터의 관측치 개수(행)를 확인할 수 없습니다.")
    return(invisible(NULL))
  }

  # --- 2. Metadata 객체 준비 및 검증 ---
  message("\n2. Metadata 객체:")
  
  internal_metadata <- NULL # 함수 내부에서 사용할 메타데이터
  
  if (is.null(metadata)) {
    message("  - `metadata`가 NULL입니다. `data`에서 추출을 시도합니다.")
    if (is_seurat) {
      internal_metadata <- data@meta.data
      message("  - `data@meta.data`를 메타데이터로 사용합니다.")
    } else if (is.data.frame(data)) {
      internal_metadata <- data
      message("  - `data` 자체를 메타데이터로 사용합니다.")
    } else {
      message("  - ERROR: `data`가 data.frame이 아니므로 메타데이터를 자동 추출할 수 없습니다. `metadata`를 제공해주세요.")
      return(invisible(NULL))
    }
  } else {
    if (!is.data.frame(metadata)) {
      message("  - ERROR: `metadata`는 data.frame이어야 합니다.")
      return(invisible(NULL))
    }
    internal_metadata <- metadata
    message("  - 제공된 `metadata`를 사용합니다.")
  }

  # 메타데이터 기본 검증
  message("  - Metadata Dim: ", nrow(internal_metadata), " rows, ", ncol(internal_metadata), " cols")
  
  # 데이터-메타데이터 관측치 일치 여부
  if (nrow(internal_metadata) != data_obs) {
    message(paste0(
      "  - FATAL ERROR: 데이터 관측치 개수 (", data_obs, 
      ")와 메타데이터 행 개수 (", nrow(internal_metadata), ")가 일치하지 않습니다."
    ))
    return(invisible(NULL))
  }
  
  # 메타데이터 개요 (str)
  str_capture <- utils::capture.output(utils::str(internal_metadata))
  message("  - Structure (상위 6줄): \n", paste("    ", str_capture[1:min(6, length(str_capture))], collapse = "\n"))
  
  # 결측치 (NA) 확인 (변수별)
  # [수정] NA 리포트 형식 변경
  meta_na_counts <- colSums(is.na(internal_metadata))
  na_cols_with_counts <- meta_na_counts[meta_na_counts > 0]
  
  if (length(na_cols_with_counts) > 0) {
    message("  - WARNING: 결측치(NA)가 발견된 컬럼:")
    # Create a formatted string for each column
    na_report_lines <- paste0(
      "    ... ", names(na_cols_with_counts), " : ", 
      na_cols_with_counts, "개"
    )
    # Print all lines
    message(paste(na_report_lines, collapse = "\n"))
  } else {
    message("  - 결측치(NA): 없음 (Good)")
  }
  
  # --- 3. 변수 (target_vars, control_vars) 검증 ---
  message("\n3. 변수 검증 (in Metadata):")
  all_vars <- unique(c(target_vars, control_vars))
  
  if (length(all_vars) == 0) {
    message("  - 검증할 변수가 지정되지 않았습니다.")
  } else {
    meta_names <- names(internal_metadata)
    
    for (var in all_vars) {
      if (!var %in% meta_names) {
        message(paste0("  - ERROR: '", var, "' 변수가 메타데이터에 존재하지 않습니다."))
      } else {
        var_class <- class(internal_metadata[[var]])[1]
        var_na_count <- sum(is.na(internal_metadata[[var]]))
        
        message(paste0("  - '", var, "': 발견 (Type: ", var_class, ", NAs: ", var_na_count, ")"))
        
        # 문자형 변수 처리
        if (var_class == "character") {
          if (auto_correct_vars) {
            internal_metadata[[var]] <- as.factor(internal_metadata[[var]])
            message(paste0("    ... (auto_correct) '", var, "'를 factor로 변환했습니다."))
          } else {
            message(paste0("    ... WARNING: '", var, "'는 character입니다. 모델링 전 factor 변환을 권장합니다."))
          }
        }
        
        # 팩터 레벨 개수 확인
        if (is.factor(internal_metadata[[var]])) {
          level_count <- length(levels(internal_metadata[[var]]))
          if (level_count > 50 && level_count > data_obs / 10) {
             message(paste0("    ... WARNING: '", var, "' 팩터의 레벨 개수 (", level_count, ")가 관측치에 비해 많습니다. (Overfitting 위험)"))
          }
        }
      }
    }
  }
  
  # --- [신규] 3b. NA 처리 (na_action) ---
  if (na_action == "remove_rows") {
    if (length(all_vars) > 0) {
      # all_vars에 지정된 컬럼들에서 NA가 없는 행만 선택
      rows_to_keep <- stats::complete.cases(internal_metadata[, all_vars, drop = FALSE])
      cleaned_metadata <- internal_metadata[rows_to_keep, ]
      
      removed_count <- nrow(internal_metadata) - nrow(cleaned_metadata)
      if (removed_count > 0) {
         message(paste0(
           "\n  - (na_action) '", paste(all_vars, collapse=","), "' 변수 기준 NA ", 
           removed_count, "개 행 제거됨. (metadata_clean에 저장)"
         ))
      }
      validated_output$metadata_clean <- cleaned_metadata
    } else {
       message("\n  - (na_action) 'remove_rows'를 선택했으나, `target_vars` 또는 `control_vars`가 지정되지 않아 NA 제거를 스킵합니다.")
       validated_output$metadata_clean <- internal_metadata
    }
  } else {
    # na_action == "report" (default)
    # NA 제거 안 함, 원본(수정됐을 수 있는) 메타데이터를 clean으로 간주
    validated_output$metadata_clean <- internal_metadata
  }
  
  # --- 4. 포뮬러 (Formula) 검증 ---
  message("\n4. 포뮬러 검증:")
  internal_formula <- formula
  
  if (is.null(internal_formula)) {
    message("  - 포뮬러가 제공되지 않았습니다.")
  } else {
    # 문자열일 경우 포뮬러로 변환
    if (is.character(internal_formula)) {
      tryCatch(
        { internal_formula <- as.formula(internal_formula) },
        error = function(e) {
          message(paste0("  - ERROR: '", internal_formula, "'는 유효한 포뮬러가 아닙니다. (", e$message, ")"))
          internal_formula <<- NULL
        }
      )
    }
    
    if (inherits(internal_formula, "formula")) {
      message(paste0("  - Formula: ", deparse1(internal_formula)))
      
      # 포뮬러 변수들이 메타데이터에 있는지 확인
      formula_vars <- all.vars(internal_formula)
      missing_vars <- setdiff(formula_vars, names(internal_metadata))
      if (length(missing_vars) > 0) {
        message(paste0("  - ERROR: 포뮬러 변수 '", paste(missing_vars, collapse = ", "), "'가 메타데이터에 없습니다."))
      }
      
      # 포뮬러 복잡성 (항 개수 vs 관측치)
      term_labels <- attr(terms(internal_formula), "term.labels")
      num_terms <- length(term_labels)
      message(paste0("  - Terms (", num_terms, "개): ", paste(term_labels, collapse = ", ")))
      
      # [수정] NA 제거시 관측치 개수 재계산
      effective_obs <- ifelse(na_action == "remove_rows", nrow(validated_output$metadata_clean), data_obs)
      if (effective_obs < num_terms * 10) {
         message(paste0(
            "  - WARNING: 유효 관측치 개수 (", effective_obs, ") 대비 포뮬러 항 (", num_terms, 
            ")이 너무 많습니다. (10:1 비율 권장, 현재: ", round(effective_obs / num_terms, 1), ":1) Overfitting 위험."
         ))
      }
      
      # 포뮬러 구조 분석 (랜덤 이펙트, 상호작용 등)
      formula_str <- deparse1(internal_formula)
      has_random <- grepl("\\|", formula_str)
      has_nested <- grepl("/", formula_str)
      has_interaction <- grepl(":", formula_str) || grepl("\\*", formula_str)
      
      if (has_random) {
        message("  - 구조: Random effects (예: (1|var)) 발견. lme4, glmmTMB, GAM 등 Mixed model에서 지원됩니다.")
      }
      if (has_nested) {
        message("  - 구조: Nested random effects (예: (1|var1/var2)) 발견. lme4 등에서 지원됩니다.")
      }
      if (has_interaction) {
        message("  - 구조: Interaction effects (예: var1:var2) 발견.")
      }
    }
  }

  # --- 5. 방법론 (Method) 호환성 검증 ---
  message("\n5. 방법론(Method) 호환성:")
  
  # [수정] 다중 method 지원
  if (is.null(method) || length(method) == 0) {
    message("  - `method`가 지정되지 않았습니다.")
  } else {
    
    # 여러 메소드를 루프로 순회
    for (m in method) {
      current_method <- tolower(m)
      message(paste0("\n  --- 검증 (Method: '", current_method, "') ---"))
      
      # 5a. ML 계열 (Random Forest, XGBoost, SVM 등)
      if (current_method %in% c("random forest", "rf", "xgboost", "svm", "lasso")) {
        message("  - 타입: Machine Learning / Regularized Regression")
        if (has_random) {
          message("  - WARNING: ML/Lasso 모델은 (1|var) 같은 랜덤 이펙트 포뮬러를 직접 지원하지 않습니다.")
          if (auto_correct_formula) {
            # lme4::findbars를 사용하여 랜덤 이펙트만 제거 시도
            if (requireNamespace("lme4", quietly = TRUE)) {
              # lme4::nobars()를 사용하여 랜덤 이펙트 항 제거
              fixed_formula <- lme4::nobars(internal_formula)
              if (!is.null(fixed_formula)) {
                internal_formula <- fixed_formula # internal_formula를 직접 수정
                message(paste0("    ... (auto_correct) 포뮬러에서 랜덤 이펙트 항을 제거했습니다: ", deparse1(internal_formula)))
                # 포뮬러가 변경되었으므로 has_random 플래그도 업데이트
                has_random <- FALSE 
              } else {
                 message("    ... (auto_correct) 랜덤 이펙트 항 제거에 실패했습니다. (포뮬러가 너무 복잡할 수 있음)")
              }
            } else {
              message("    ... (auto_correct) 포뮬러 수정을 위해 'lme4' 패키지가 필요합니다.")
            }
          } else {
            message("    ... `auto_correct_formula = TRUE`로 설정하여 랜덤 이펙트 항을 자동으로 제거할 수 있습니다.")
          }
        }
      }
      
      # 5b. 전통적 통계 검정 (wilcoxon, t.test, limma)
      else if (current_method %in% c("wilcoxon", "t.test", "limma")) {
        message("  - 타입: 통계적 가설 검정")
        if (has_random && current_method != "limma") {
          message(paste0("  - WARNING: '", current_method, "'는 랜덤 이펙트를 지원하지 않습니다."))
        }
        if (current_method == "limma" && has_random) {
           message("  - INFO: 'limma'는 `duplicateCorrelation()`을 통해 랜덤 이펙트(반복측정)를 처리할 수 있습니다.")
        }
         if (num_terms > 3) { # num_terms는 포뮬러에서 계산됨
           message("  - WARNING: '", current_method, "'는 보통 'target ~ group' 형태의 간단한 포뮬러를 사용합니다.")
         }
      }
      
      # 5c. 기본 선형 모델 (lm, glm)
      else if (current_method %in% c("lm", "glm")) {
         message("  - 타입: (일반) 선형 모델")
         if (has_random) {
           message("  - ERROR: 'lm'/'glm'은 랜덤 이펙트를 지원하지 않습니다. 'lme4' 패키지의 `lmer`/`glmer`를 사용해야 합니다.")
         }
      }
      
      # 5d. 혼합 효과 모델 (lmer, glmer, lme, gam)
      else if (current_method %in% c("lmer", "glmer", "lme", "lme4", "gam", "gamm")) {
         message("  - 타입: 혼합 효과 / 일반화 가법 모델 (Mixed / GAM)")
         if (has_random) {
           message("  - INFO: 랜덤 이펙트를 지원하는 모델입니다. (Good)")
         } else {
           message("  - INFO: 랜덤 이펙트가 없지만, 'lm'/'glm'/'gam'으로도 분석 가능할 수 있습니다.")
         }
      }
      
      # 5e. 차원 축소 (PCA)
      else if (current_method %in% c("pca", "pca_loadings")) {
         message("  - 타입: 차원 축소")
         message("  - INFO: 'pca'는 보통 포뮬러를 사용하지 않으며 (혹은 `~.` 사용), 숫자형 데이터(matrix)가 필요합니다.")
      }
      
      else {
        message(paste0("  - '", current_method, "'는 알 수 없거나 검증 로직에 없는 방법론입니다."))
      }
    } # end for method loop
  }

  message("\n--- 종료: 검증 리포트 ---")

  # --- 6. 결과 반환 ---
  # auto_correct_vars가 적용된 internal_metadata를 반환 리스트에 할당
  validated_output$metadata <- internal_metadata
  # auto_correct_formula가 적용된 internal_formula를 반환 리스트에 할당
  validated_output$formula <- internal_formula
  
  # metadata_clean은 위에서 (Section 3b) 이미 할당됨
  
  # 사용자가 이 리스트를 받아 후속 작업에 사용해야 함
  return(invisible(validated_output))
}