



# LMM analysis suite -----
# PART 1: 데이터 준비 및 설정

#' 분석 설정 객체 생성
#'
#' 분석 파이프라인 전반에서 사용될 메타데이터의 컬럼명(변수명)을
#' 표준화된 리스트 객체로 생성합니다.
#'
#' @param patient 환자 ID 컬럼명. (기본값: "patient_id")
#' @param drug 약물 정보 컬럼명. (기본값: "drug")
#' @param timepoint 시간점(e.g., "pre", "post") 컬럼명. (기본값: "timepoint")
#' @param ck CK status (e.g., "CK+", "CK-") 컬럼명. (기본값: "ck_status")
#' @param response 반응 여부(e.g., "R", "NR") 컬럼명. (기본값: "response")
#' @param aoi AOI (Area of Interest) ID 컬럼명. (기본값: "aoi_id")
#'
#' @return 컬럼명 설정이 담긴 명명된 리스트(named list).
#' @export
#' @examples
#' config <- create_analysis_config(patient = "Subject", drug = "Treatment")
create_analysis_config <- function(
  patient = "patient_id",
  drug = "drug", 
  timepoint = "timepoint",
  ck = "ck_status",
  response = "response",
  aoi = "aoi_id"
) {
  list(
    patient = patient,
    drug = drug,
    timepoint = timepoint,
    ck = ck,
    response = response,
    aoi = aoi
  )
}

#' GeoMx 데이터 준비 및 Seurat 객체 생성
#'
#' Excel/CSV 파일 또는 R 객체(matrix, data.frame)로부터 count 및 metadata를 읽어
#' 분석을 위한 Seurat 객체를 생성합니다.
#'
#' @param count_file Count matrix 파일 경로 (Excel 또는 CSV).
#' @param metadata_file Metadata 파일 경로 (Excel 또는 CSV).
#' @param count_matrix (선택 사항) 파일 대신 사용할 count matrix 객체.
#' @param metadata (선택 사항) 파일 대신 사용할 metadata data.frame 객체.
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param q3_normalize Q3 정규화가 이미 수행되었는지 여부.
#'        `TRUE` (기본값)이면 별도 정규화를 생략하고, `FALSE`이면 LogNormalize를 수행합니다.
#'
#' @return 전처리된 Seurat 객체. `condition_full`과 `condition_simple` 메타데이터가 추가됩니다.
#' @importFrom Seurat CreateSeuratObject NormalizeData
#' @importFrom openxlsx read.xlsx
#' @export
prepare_geomx_data <- function(count_file = NULL, 
                             metadata_file = NULL,
                             count_matrix = NULL,
                             metadata = NULL,
                             config = create_analysis_config(),
                             q3_normalize = TRUE) {
  
  # 파일에서 읽기 또는 직접 입력
  if (!is.null(count_file)) {
    if (grepl("\\.xlsx$", count_file)) {
      count_matrix <- openxlsx::read.xlsx(count_file, sheet = 1, rowNames = TRUE)
    } else {
      count_matrix <- read.csv(count_file, row.names = 1)
    }
  }
  
  if (!is.null(metadata_file)) {
    if (grepl("\\.xlsx$", metadata_file)) {
      metadata <- openxlsx::read.xlsx(metadata_file, sheet = 2)
    } else {
      metadata <- read.csv(metadata_file)
    }
  }
  
  # AOI 이름 매칭
  rownames(metadata) <- metadata[[config$aoi]]
  
  # Seurat 객체 생성
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(count_matrix),
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
  
  # Q3 정규화는 이미 되어있다고 가정
  if (!q3_normalize) {
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  }
  
  # 조합 변수 생성 (빠른 스크리닝용)
  seurat_obj$condition_full <- paste0(
    seurat_obj@meta.data[[config$timepoint]], "_",
    seurat_obj@meta.data[[config$drug]], "_", 
    seurat_obj@meta.data[[config$response]]
  )
  
  seurat_obj$condition_simple <- paste0(
    seurat_obj@meta.data[[config$timepoint]], "_",
    seurat_obj@meta.data[[config$drug]]
  )
  
  return(seurat_obj)
}

# PART 2: 빠른 스크리닝 (Wilcoxon)

#' 빠른 유전자 스크리닝 (Wilcoxon Test)
#'
#' `FindAllMarkers` (Wilcoxon rank sum test)를 사용하여 그룹 간 차등 발현 유전자를
#' 신속하게 스크리닝합니다.
#'
#' @param seurat_obj Seurat 객체.
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param comparison_type 비교할 그룹 기준. "pre_post", "drug", "response", "combined" 중 선택.
#' @param ck_subset (선택 사항) "CK+" 또는 "CK-"를 지정하여 특정 CK status로 필터링.
#' @param min_pct `FindAllMarkers`의 `min.pct` 파라미터. (기본값: 0.1)
#' @param logfc_threshold `FindAllMarkers`의 `logfc.threshold` 파라미터. (기본값: 0.25)
#' @param top_n 반환할 상위 유전자 개수. (기본값: 1000)
#' @param use_adj `TRUE`이면 `p_val_adj` 기준, `FALSE` (기본값)이면 `p_val` 기준으로 정렬.
#'
#' @return 리스트:
#' \item{markers}{`FindAllMarkers` 결과 (data.frame)}
#' \item{top_genes}{상위 유전자 벡터}
#' \item{seurat_obj}{필터링이 적용된 Seurat 객체 (if `ck_subset` is used)}
#' @importFrom Seurat Idents FindAllMarkers
#' @importFrom dplyr arrange
#' @importFrom rlang sym
#' @export
quick_screen_genes <- function(seurat_obj,
                               config = create_analysis_config(),
                               comparison_type = "pre_post",
                               ck_subset = NULL,
                               min_pct = 0.1,
                               logfc_threshold = 0.25,
                               top_n = 1000,
                               use_adj = FALSE) {
  
  # CK 서브셋 필터링
  if (!is.null(ck_subset)) {
    seurat_obj <- subset(seurat_obj, 
                         subset = !!sym(config$ck) == ck_subset)
  }
  
  # 비교 그룹 설정
  if (comparison_type == "pre_post") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$timepoint]]
  } else if (comparison_type == "drug") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$drug]]
  } else if (comparison_type == "response") {
    Idents(seurat_obj) <- seurat_obj@meta.data[[config$response]]
  } else if (comparison_type == "combined") {
    Idents(seurat_obj) <- seurat_obj$condition_full
  }
  
  # FindAllMarkers 실행
  markers <- FindAllMarkers(
    seurat_obj,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = FALSE,
    return.thresh = 1  # 모든 유전자 반환
  )
  
  # p-value 선택 및 상위 유전자 선택
  if (use_adj) {
    top_genes <- markers %>%
      arrange(p_val_adj) %>%
      head(top_n) %>%
      pull(gene) %>%
      unique()
  } else {
    top_genes <- markers %>%
      arrange(p_val) %>%
      head(top_n) %>%
      pull(gene) %>%
      unique()
  }
  
  return(list(
    markers = markers,
    top_genes = top_genes,
    seurat_obj = seurat_obj
  ))
}

# PART 3: Linear Mixed Model 분석

#' 단일 유전자에 대한 LMM (Linear Mixed Model) 적합
#'
#' (내부 헬퍼 함수) `run_lmm_multiple_genes` 내에서 단일 유전자에 대해
#' `lmer`를 사용하여 LMM을 적합합니다.
#'
#' @param gene_expr 단일 유전자의 발현값 벡터 (numeric).
#' @param metadata Seurat 객체의 `meta.data` 슬롯 (data.frame).
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param formula_components (선택 사항) LMM 공식을 동적으로 생성하기 위한 리스트.
#'        (e.g., `list(fixed = c("drug", "timepoint"), random = "(1|patient)")`)
#' @param use_config_names `formula_components` 사용 시, "patient"와 같은 일반 용어를
#'        `config$patient` (`"patient_id"`)로 자동 매핑할지 여부. (기본값: `TRUE`)
#'
#' @return 모델 적합 결과 리스트.
#' \item{model}{`lmerMod` 객체}
#' \item{effects}{계수 요약 (data.frame)}
#' \item{anova}{ANOVA 테이블}
#' \item{converged}{모델 수렴 여부 (logical)}
#' \item{formula}{모델에 사용된 공식 (string)}
#' \item{all_drugs}{데이터에 있는 모든 약물}
#' \item{ref_drug}{모델의 참조(reference) 약물}
#' @importFrom lmerTest anova lmer
#' @importFrom stats as.formula relevel
fit_lmm_single_gene <- function(gene_expr,
                                metadata,
                                config = create_analysis_config(),
                                formula_components = NULL,
                                use_config_names = TRUE) {
  
  # 데이터프레임 생성
  df <- cbind(
    data.frame(Expression = gene_expr),
    metadata
  )
  
  # Formula 생성
  if (is.null(formula_components)) {
    # 기본값: config 변수명 사용
    fixed_effects <- c(config$drug, config$timepoint, config$response)
    interactions <- c(
      paste(config$drug, config$timepoint, sep = ":"),
      paste(config$drug, config$response, sep = ":"),
      paste(config$timepoint, config$response, sep = ":"),
      paste(config$drug, config$timepoint, config$response, sep = ":")
    )
    random <- paste0("(1|", config$patient, ")")
    
    formula_str <- paste0(
      "Expression ~ ",
      paste(c(fixed_effects, interactions), collapse = " + "),
      " + ", random
    )
  } else {
    # formula_components 제공 시
    if (use_config_names) {
      # config 매핑 사용
      mapping <- list(
        "patient" = config$patient,
        "drug" = config$drug,
        "timepoint" = config$timepoint,
        "response" = config$response,
        "ck" = config$ck
      )
      
      # fixed effects 매핑
      fixed_mapped <- sapply(formula_components$fixed, function(x) {
        if (x %in% names(mapping)) mapping[[x]] else x
      })
      
      # interactions 매핑
      if (!is.null(formula_components$interactions)) {
        interactions_mapped <- sapply(formula_components$interactions, function(x) {
          parts <- strsplit(x, ":")[[1]]
          mapped_parts <- sapply(parts, function(p) {
            if (p %in% names(mapping)) mapping[[p]] else p
          })
          paste(mapped_parts, collapse = ":")
        })
      } else {
        interactions_mapped <- NULL
      }
      
      # random effects 매핑
      random_mapped <- formula_components$random
      for (key in names(mapping)) {
        random_mapped <- gsub(key, mapping[[key]], random_mapped)
      }
      
      formula_str <- paste0(
        "Expression ~ ",
        paste(c(fixed_mapped, interactions_mapped), collapse = " + "),
        " + ", random_mapped
      )
    } else {
      # 직접 사용 (매핑 없이)
      fixed_effects <- formula_components$fixed
      interactions <- formula_components$interactions
      formula_str <- paste0(
        "Expression ~ ",
        paste(c(fixed_effects, interactions), collapse = " + "),
        " + ", formula_components$random
      )
    }
  }
  
  # 모델 적합
  tryCatch({
    # Factor 레벨 재정렬 (reference level 설정)
    if (config$drug %in% names(df)) {
      df[[config$drug]] <- relevel(as.factor(df[[config$drug]]), ref = levels(as.factor(df[[config$drug]]))[1])
    }
    if (config$response %in% names(df)) {
      df[[config$response]] <- relevel(as.factor(df[[config$response]]), ref = levels(as.factor(df[[config$response]]))[1])
    }
    
    model <- lmer(as.formula(formula_str), data = df, REML = FALSE)
    
    # 결과 추출
    coef_summary <- summary(model)$coefficients
    anova_result <- anova(model)
    
    # 모든 약물 대비 추출 (reference level 포함)
    all_drugs <- unique(df[[config$drug]])
    ref_drug <- levels(as.factor(df[[config$drug]]))[1]
    
    # 주요 효과 계산
    effects <- data.frame(
      term = rownames(coef_summary),
      estimate = coef_summary[, "Estimate"],
      std_error = coef_summary[, "Std. Error"],
      t_value = coef_summary[, "t value"],
      p_value = coef_summary[, "Pr(>|t|)"],
      ref_drug = ref_drug
    )
    
    return(list(
      model = model,
      effects = effects,
      anova = anova_result,
      converged = TRUE,
      formula = formula_str,
      all_drugs = all_drugs,
      ref_drug = ref_drug
    ))
    
  }, error = function(e) {
    return(list(
      converged = FALSE,
      error = e$message,
      formula = formula_str
    ))
  })
}

#' 다중 유전자에 대한 LMM 병렬 수행
#'
#' 지정된 유전자 목록에 대해 `fit_lmm_single_gene` 함수를 병렬로 실행하여
#' LMM 분석을 수행합니다.
#'
#' @param seurat_obj Seurat 객체.
#' @param genes (선택 사항) 분석할 유전자 이름 벡터.
#'        `NULL` (기본값)인 경우, 경고와 함께 첫 100개 유전자를 사용합니다.
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#' @param formula_components (선택 사항) `fit_lmm_single_gene`로 전달될 LMM 공식 리스트.
#' @param use_config_names `formula_components` 사용 시, `config` 매핑 사용 여부. (기본값: `TRUE`)
#' @param n_cores 병렬 처리에 사용할 CPU 코어 수. (기본값: 16)
#' @param verbose 진행 상황 메시지를 출력할지 여부. (기본값: `TRUE`)
#'
#' @return 리스트:
#' \item{raw_results}{각 유전자에 대한 `fit_lmm_single_gene`의 원시 결과 리스트}
#' \item{summary}{`summarize_lmm_results`로 요약된 전체 결과 (data.frame)}
#' \item{converged_genes}{수렴에 성공한 유전자 수}
#' \item{total_genes}{시도한 총 유전자 수}
#' @importFrom Seurat GetAssayData
#' @importFrom parallel makeCluster clusterEvalQ clusterExport parLapply stopCluster
#' @export
run_lmm_multiple_genes <- function(seurat_obj,
                                 genes = NULL,
                                 config = create_analysis_config(),
                                 formula_components = NULL,
                                 use_config_names = TRUE,
                                 n_cores = 16,
                                 verbose = TRUE) {
  
  if (is.null(genes)) {
    genes <- rownames(seurat_obj)[1:100]  # 기본값: 상위 100개
    warning("No genes specified. Using first 100 genes.")
  }
  
  # 발현 매트릭스와 메타데이터 추출
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  metadata <- seurat_obj@meta.data
  
  # 진행상황 표시
  if (verbose) {
    message(sprintf("Running LMM for %d genes using %d cores...", 
                    length(genes), n_cores))
  }
  
  # 병렬 처리
  if (n_cores > 1) {
    cl <- makeCluster(n_cores, type="PSOCK")
    clusterEvalQ(cl, {
      library(lme4)
      library(lmerTest)
    })
    clusterExport(cl, c("fit_lmm_single_gene", "config", "metadata", 
                        "expr_matrix", "formula_components", "use_config_names"),
                  envir = environment())
    
    results <- parLapply(cl, genes, function(gene) {
      if (gene %in% rownames(expr_matrix)) {
        result <- fit_lmm_single_gene(
          gene_expr = as.numeric(expr_matrix[gene, ]),
          metadata = metadata,
          config = config,
          formula_components = formula_components,
          use_config_names = use_config_names
        )
        result$gene <- gene
        return(result)
      } else {
        return(list(gene = gene, converged = FALSE, 
                    error = "Gene not found"))
      }
    })
    
    stopCluster(cl)
  } else {
    # 순차 처리
    results <- lapply(genes, function(gene) {
      if (verbose && which(genes == gene) %% 100 == 0) {
        message(sprintf("  Processing gene %d/%d", 
                        which(genes == gene), length(genes)))
      }
      
      if (gene %in% rownames(expr_matrix)) {
        result <- fit_lmm_single_gene(
          gene_expr = as.numeric(expr_matrix[gene, ]),
          metadata = metadata,
          config = config,
          formula_components = formula_components
        )
        result$gene <- gene
        return(result)
      } else {
        return(list(gene = gene, converged = FALSE, 
                    error = "Gene not found"))
      }
    })
  }
  
  names(results) <- genes
  
  # 결과 요약
  summary_df <- summarize_lmm_results(results, config)
  
  return(list(
    raw_results = results,
    summary = summary_df,
    converged_genes = sum(sapply(results, function(x) x$converged)),
    total_genes = length(genes)
  ))
}

#' LMM 결과 요약 및 p-value 보정
#'
#' (내부 헬퍼 함수) `run_lmm_multiple_genes`에서 생성된 원시 결과 리스트를
#' 하나의 요약 data.frame으로 통합하고, 용어(term)별로 p-value를 보정(BH)합니다.
#'
#' @param lmm_results `run_lmm_multiple_genes`의 `$raw_results` 객체 (리스트).
#' @param config `create_analysis_config()`로 생성된 변수명 설정 리스트.
#'
#' @return 모든 유전자와 모든 효과(term)를 포함하는 요약 data.frame.
#'         `p_adj` 및 `significant` 컬럼이 추가됩니다.
#' @importFrom dplyr bind_rows group_by mutate ungroup
#' @importFrom stats p.adjust
summarize_lmm_results <- function(lmm_results, config) {
  # 수렴한 모델만 처리
  converged <- lmm_results[sapply(lmm_results, function(x) x$converged)]
  
  if (length(converged) == 0) {
    warning("No models converged!")
    return(NULL)
  }
  
  # 모든 효과 수집
  all_effects <- do.call(rbind, lapply(names(converged), function(gene) {
    effects <- converged[[gene]]$effects
    effects$gene <- gene
    return(effects)
  }))
  
  # p-value 보정
  all_effects <- all_effects %>%
    group_by(term) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant = p_adj < 0.05
    ) %>%
    ungroup()
  
  # Reference drug 정보 추가
  if ("ref_drug" %in% names(converged[[1]]$effects)) {
    ref_drugs <- unique(sapply(converged, function(x) x$ref_drug))
    attr(all_effects, "ref_drug") <- ref_drugs[1]
  }
  
  return(all_effects)
}
