
# --------- analyses: LISI, PERMANOVA ---------

#' LISI 점수에 대해 의사 반복을 피해 통계 검정 수행
#'
#' @param lisi_scores LISI 점수 (compute_lisi의 결과. Rowname이 cell barcode)
#' @param metadata Seurat 객체의 메타데이터 (is@meta.data)
#' @param lisi_col_name LISI 점수 컬럼명 (e.g., "GEM", "batch")
#' @param patient_var 환자/샘플 ID가 있는 메타데이터 컬럼명 (e.g., "patient_ID")
#' @param celltype_var 세포 유형 정보가 있는 메타데이터 컬럼명 (e.g., "cell_type")
#' @param biological_group_var 비교하려는 그룹 (e.g., "disease_status", "set")
#'
#' @return cell_type별 통계 검정 결과 (adj_p_value로 정렬됨)
#' @export
TestLISI = function(lisi_scores, 
                    metadata, 
                    lisi_col_name, 
                    patient_var, 
                    celltype_var, 
                    biological_group_var) {
  
  # 1. 데이터 준비 (메타데이터와 LISI 점수 결합)
  
  # LISI 점수 컬럼이 있는지 확인
  if (!lisi_col_name %in% colnames(lisi_scores)) {
    stop(paste0("Error: '", lisi_col_name, "' 컬럼이 lisi_scores에 없습니다."))
  }
  
  # LISI 점수 컬럼을 'LISI_score'라는 공통 이름으로 변경
  lisi_scores_renamed <- lisi_scores %>%
    select(LISI_score = !!sym(lisi_col_name))
  
  # 메타데이터와 LISI 점수 결합 (rownames 기준)
  meta_data_with_lisi <- metadata %>%
    rownames_to_column("cell_barcode_temp") %>%
    left_join(
      lisi_scores_renamed %>% rownames_to_column("cell_barcode_temp"),
      by = "cell_barcode_temp"
    ) %>%
    select(-cell_barcode_temp) # 임시 바코드 컬럼 제거

  # 2. 세포 유형별 루프 시작
  all_cell_types <- unique(meta_data_with_lisi[[celltype_var]])
  
  message(paste("Starting LISI statistical test for", length(all_cell_types), "cell types..."))
  
  stats_results <- purrr::map_dfr(all_cell_types, function(current_cell_type) {
    
    # 3. 환자 레벨로 LISI 점수 집계 (Aggregation)
    # (핵심) 각 환자별 LISI 점수의 중앙값을 계산
    median_lisi_per_patient <- meta_data_with_lisi %>%
      filter(!!sym(celltype_var) == current_cell_type) %>%
      filter(!is.na(!!sym(biological_group_var))) %>% # 비교 그룹이 NA인 셀 제외
      group_by(!!sym(patient_var), !!sym(biological_group_var)) %>%
      summarise(
        median_LISI = median(LISI_score, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      # 한 환자가 여러 메타데이터 행을 가질 경우(e.g. 여러 배치) 
      # group_by/summarise가 이미 환자당 1줄로 만들지만,
      # 만일을 대비해 환자 ID로 중복 제거
      distinct(!!sym(patient_var), .keep_all = TRUE) 
      
    # 4. 통계 검정을 위한 데이터 유효성 검사
    if (nrow(median_lisi_per_patient) < 3 || 
        length(unique(median_lisi_per_patient[[biological_group_var]])) < 2) {
      
      message(paste("  [Skipping]", current_cell_type, "| Not enough samples or < 2 groups"))
      return(tibble(
        cell_type = current_cell_type,
        statistic = NA,
        p_value = NA,
        error = "Not enough samples or only one group"
      ))
    }
    
    # 5. Wilcoxon Rank-Sum Test 수행
    test_formula <- as.formula(paste("median_LISI ~", biological_group_var))
    
    tryCatch({
      wilcox_res <- wilcox.test(test_formula, data = median_lisi_per_patient)
      
      return(tibble(
        cell_type = current_cell_type,
        statistic_W = wilcox_res$statistic,
        p_value = wilcox_res$p.value,
        error = NA
      ))
    }, error = function(e) {
      # (e.g., 한 그룹의 모든 값이 동일할 때 오류 발생 가능)
      message(paste("  [Error]", current_cell_type, ":", e$message))
      return(tibble(
        cell_type = current_cell_type,
        statistic_W = NA,
        p_value = NA,
        error = e$message
      ))
    })
  })
  
  # 6. P-value 보정 (BH) 및 정렬
  stats_results <- stats_results %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    arrange(adj_p_value)
    
  message("LISI statistical test finished.")
  print(stats_results)
  return(stats_results)
}

#' @export
PERMANOVA=function(sobj, patient_var, celltype_var, group_var, assay="RNA", na.drop=TRUE){
  # 환자 ID가 들어있는 컬럼명 (예: "p1", "p2"...)
  # patient_var
  # 세포 유형 정보가 있는 컬럼명
  # celltype_var <- "seurat_clusters"
  # 테스트하려는 생물학적 그룹 컬럼명 (예: "set1", "set2", "set3")
  # group_var <- "set"
  # 사용할 Assay
  # assay <- "SCT"
  # na.drop: group_var에 NA가 있는 셀을 제거할지 여부
  
  # --- 분석 시작 ---
  
  # 1. 분석할 모든 세포 유형 목록 가져오기
  all_cell_types <- unique(sobj@meta.data[[celltype_var]])
  
  # 2. 각 세포 유형에 대해 반복적으로 PERMANOVA 수행
  # map_dfr은 결과를 data.frame으로 깔끔하게 묶어줍니다.
  pb_permanova_results <- purrr::map_dfr(all_cell_types, function(current_cell_type) {
    message(paste("Processing:", current_cell_type))
    # 3. 현재 세포 유형으로만 Seurat 객체 필터링
    seurat_subset <- subset(sobj, !!sym(celltype_var) == current_cell_type)
    
    # --- [수정] group_var에 NA가 있으면 해당 셀 제거 (na.drop=TRUE일 경우) ---
    if (na.drop) {
      seurat_subset <- subset(seurat_subset, !is.na(!!sym(group_var)))
    }
    
    # --- [추가] 필터링 후 셀이 남아있는지 확인 ---
    if (nrow(seurat_subset@meta.data) == 0) {
      message(paste("  [Skipping]", current_cell_type, "| No cells remaining after NA drop"))
      return(data.frame(
        cell_type = current_cell_type,
        R2 = NA,
        p_value = NA,
        error = "No cells remaining after NA drop"
      ))
    }
    # --- [추가 끝] ---
    
    # 4. 환자(patient_var)별로 Pseudo-bulk 수행
    pb_subset <- AggregateExpression(
      seurat_subset,
      group.by = patient_var,
      assays = assay,
      slot = "data" # SCT/RNA의 정규화된 데이터를 사용
    )
    
    # 5. 매트릭스 준비 (샘플 x 유전자)
    pb_matrix_t <- t(pb_subset[[assay]])
    
    # 6. Pseudo-bulk 샘플(환자)에 대한 메타데이터 생성
    # 환자 ID와 그 환자가 속한 'set' 그룹을 매칭
    meta_pb_subset <- seurat_subset@meta.data %>%
      select(!!sym(patient_var), !!sym(group_var)) %>%
      filter(!is.na(!!sym(patient_var))) %>% #NA 값을 가진 행을 미리 제거합니다.
      distinct() # 환자-그룹 매핑 정보만 남김

    # 매트릭스 row 순서와 메타데이터 순서 맞추기
    # AggregateExpression은 숫자로 시작하는 ID 앞에 'g'를 붙일 수 있으므로,
    # join_key를 만들어 두 경우 모두 처리합니다. (숫자이거나, 캐릭터이거나.)
    
    # --- [수정] dplyr의 NSE(Non-Standard Evaluation) 오류를 피하기 위해 base R로 변경 ---
    patient_ids <- as.character(meta_pb_subset[[patient_var]])
    needs_g_prefix <- grepl("^[0-9]", patient_ids)
    join_keys <- ifelse(needs_g_prefix, paste0("g", patient_ids), patient_ids)
    
    meta_pb_subset$join_key <- join_keys
    # column_to_rownames를 사용하기 전에 기존의 행 이름(cell barcodes)을 제거합니다.
    rownames(meta_pb_subset) <- NULL
    meta_pb_subset <- tibble::column_to_rownames(meta_pb_subset, "join_key")
    
    # AggregateExpression이 생성한 실제 샘플 이름으로 메타데이터 순서를 맞춥니다.
    meta_pb_subset <- meta_pb_subset[rownames(pb_matrix_t), , drop = FALSE]


    # # 매트릭스 row 순서와 메타데이터 순서 맞추기
    # rownames(meta_pb_subset) <- paste0("g",meta_pb_subset[[patient_var]])
    # meta_pb_subset <- meta_pb_subset[rownames(pb_matrix_t), ]
    
    # (오류 방지) 해당 세포 유형에 샘플이 너무 적거나 그룹이 1개면 스킵
    if (nrow(meta_pb_subset) < 3 || length(unique(meta_pb_subset[[group_var]])) < 2) {
      # --- [확인용 코드 추가] ---
      message(paste(" [Skipping]", current_cell_type,
                    "| Samples:", nrow(meta_pb_subset),
                    "| Unique Groups:", length(unique(meta_pb_subset[[group_var]]))))
      # --- [추가 끝] ---
      return(data.frame(
        cell_type = current_cell_type,
        R2 = NA,
        p_value = NA,
        error = "Not enough samples or only one group"
      ))
    }
    
    # 7. 유클리드 거리 계산
    dist_pb <- dist(pb_matrix_t, method = "euclidean")
    
    # 8. PERMANOVA 수행 (formula: 거리 ~ 그룹)
    set_formula <- as.formula(paste("dist_pb ~", group_var))
    
    res <- adonis2(
      set_formula,
      data = meta_pb_subset,
      permutations = 999)
    
    # 9. 결과 추출
    data.frame(
      cell_type = current_cell_type,
      R2 = res$R2[1], # 'set' 그룹이 설명하는 분산의 양
      p_value = res$`Pr(>F)`[1], # 'set' 그룹 간 차이의 p-value
      error = NA
    )}
  )
  
  # --- 10. P-value 보정 (Benjamini-Hochberg) ---
  pb_permanova_results$adj_p_value <- p.adjust(pb_permanova_results$p_value, method = "BH")
  
  # --- 11. adj_p_value 기준으로 정렬 ---
  pb_permanova_results <- pb_permanova_results %>%
    dplyr::arrange(adj_p_value)
  
  # --- 최종 결과 확인 ---
  print(pb_permanova_results)
  
  return(pb_permanova_results)
}


# --------- analyses: CLMM, Spearman ---------

#' Seurat 객체에서 유전자별 CLMM 분석 수행 (v2: 기능 추가)
#'
#' @param sobj Seurat 객체
#' @param ordinal_var 메타데이터의 순서형 반응 변수 (예: "heal")
#' @param patient_col 메타데이터의 환자 ID(random effect) (예: "emrid")
#' @param assay 사용할 Seurat assay (기본값: "RNA")
#' @param layer 사용할 Seurat layer (v4/v3의 'slot') (기본값: "data")
#' @param top_n 분석할 상위 유전자 개수 (기본값: NULL = 모든 유전자)
#' @param sort_by 'top_n'을 사용할 경우 정렬 기준 ("variance" 또는 "mean")
#' @param decreasing 'top_n' 정렬 시 내림차순 (기본값: TRUE)
#' @param ... clmm() 함수에 전달할 추가 인자 (예: control = clmm.control(maxIter = 200))
#'
#' @return 'gene' 및 CLMM 결과 (Estimate, p-value 등) data.frame
#'
run_clmm_analysis <- function(sobj, ordinal_var, patient_col, assay = "RNA", layer = "data", 
                              top_n = NULL, sort_by = "variance", decreasing = TRUE, ...) {
  
  # 0. 필수 패키지 확인
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("'ordinal' 패키지가 필요합니다: install.packages('ordinal')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("'dplyr' 패키지가 필요합니다: install.packages('dplyr')")
  }
  
  # 1. 메타데이터 준비
  meta <- sobj@meta.data
  unique_levels <- sort(unique(meta[[ordinal_var]]))
  meta[[ordinal_var]] <- factor(meta[[ordinal_var]], levels = unique_levels, ordered = TRUE)
  cat("Note: '", ordinal_var, "'를 ordered factor로 변환했습니다: ", paste(unique_levels, collapse = "<"), "\n")
  meta[[patient_col]] <- as.factor(meta[[patient_col]])
  
  # 2. 유전자 발현 데이터 및 유전자 목록 준비
  expr_matrix <- GetAssayData(sobj, assay = assay, slot = layer)
  
  if (!is.null(top_n)) {
    cat(paste("'top_n'이 설정되었습니다.", top_n, "개 유전자를", sort_by, "기준으로 정렬합니다...\n"))
    
    if (sort_by == "variance") {
      metrics <- apply(expr_matrix, 1, var)
    } else if (sort_by == "mean") {
      metrics <- rowMeans(expr_matrix)
    } else {
      stop("'sort_by'는 'variance' 또는 'mean'이어야 합니다.")
    }
    
    ordered_indices <- order(metrics, decreasing = decreasing)
    genes_to_test <- rownames(expr_matrix)[ordered_indices][1:min(top_n, length(metrics))]
    
  } else {
    genes_to_test <- rownames(expr_matrix)
  }
  
  total_genes <- length(genes_to_test)
  cat("총", total_genes, "개의 유전자에 대해 CLMM 모델을 실행합니다...\n")
  
  # 3. clmm 모델 실행 (for loop로 변경: 진행 상황 및 시간 측정)
  clmm_results_list <- list()
  overall_start_time <- Sys.time()
  
  for (i in seq_along(genes_to_test)) {
    gene <- genes_to_test[i]
    
    model_data <- data.frame(
      ord_var = meta[[ordinal_var]],
      pat_var = meta[[patient_col]],
      gene_expr = as.numeric(expr_matrix[gene, ])
    )
    
    # 모델 피팅 (오류 발생 시 NULL 반환)
    model_fit <- tryCatch({
      # '...' 인자를 통해 clmm.control(maxIter=200) 등 전달 가능
      ordinal::clmm(ord_var ~ gene_expr + (1|pat_var), data = model_data, ...)
    }, error = function(e) {
      NULL
    })
    
    # 결과 추출
    if (!is.null(model_fit)) {
      coef_summary <- summary(model_fit)$coefficients
      if ("gene_expr" %in% rownames(coef_summary)) {
        clmm_results_list[[gene]] <- coef_summary["gene_expr", ]
      } else {
        clmm_results_list[[gene]] <- rep(NA, 4)
      }
    } else {
      clmm_results_list[[gene]] <- rep(NA, 4)
    }
    
    # 4. 진행 상황 및 예상 시간 리포트
    if (i == 100) {
      first_100_time <- Sys.time()
      time_per_100 <- as.numeric(difftime(first_100_time, overall_start_time, units = "secs"))
      total_estimated_secs <- (time_per_100 / 100) * total_genes
      total_estimated_mins <- total_estimated_secs / 60
      
      cat(paste0("[진행 상황] 첫 100개 유전자 완료 (", round(time_per_100, 1), "초 소요)\n"))
      cat(paste0("   > 총 예상 소요 시간: 약 ", round(total_estimated_mins, 1), " 분\n"))
      
    } else if (i > 100 && i %% 100 == 0) {
      current_time <- Sys.time()
      elapsed_mins <- as.numeric(difftime(current_time, overall_start_time, units = "mins"))
      cat(paste0("[진행 상황] ", i, " / ", total_genes, " (", round(i/total_genes*100), "%) 완료. (", round(elapsed_mins, 1), "분 경과)\n"))
    }
  } # end for loop
  
  # 5. 결과 data.frame으로 변환
  clmm_results_df <- do.call(rbind, clmm_results_list)
  rownames(clmm_results_df) <- names(clmm_results_list) # genes_to_test 순서
  colnames(clmm_results_df) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  clmm_results_df <- as.data.frame(clmm_results_df)
  clmm_results_df$gene <- rownames(clmm_results_df)
  
  clmm_results_df <- clmm_results_df[!is.na(clmm_results_df$`Pr(>|z|)`), ]
  clmm_results_df <- clmm_results_df %>% 
    dplyr::arrange(`Pr(>|z|)`) %>%
    dplyr::select(gene, dplyr::everything())
  
  cat("CLMM 분석 완료.\n")
  return(clmm_results_df)
}

#' Seurat 객체에서 유전자별 스피어만 상관관계 분석 (Pseudobulk) 수행 (v2: 버그 수정 및 기능 추가)
#'
#' @param sobj Seurat 객체
#' @param ordinal_var 메타데이터의 순서형 변수 (예: "heal")
#' @param patient_col 메타데이터의 환자 ID (예: "emrid")
#' @param assay 사용할 Seurat assay (기본값: "RNA")
#' @param layer 사용할 Seurat layer (v4/v3의 'slot') (기본값: "data")
#' @param top_n 분석할 상위 유전자 개수 (기본값: NULL = 모든 유전자)
#' @param sort_by 'top_n'을 사용할 경우 정렬 기준 ("variance" 또는 "mean")
#' @param decreasing 'top_n' 정렬 시 내림차순 (기본값: TRUE)
#'
#' @return 'gene', 'rho', 'p.value'를 포함하는 data.frame
#'
run_spearman_pseudobulk <- function(sobj, ordinal_var, patient_col, assay = "RNA", layer = "data",
                                    top_n = NULL, sort_by = "variance", decreasing = TRUE) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("'dplyr' 패키지가 필요합니다: install.packages('dplyr')")
  }

  # 1. (오류 수정) 숫자로 시작하는 ID에 'g'를 덧붙인 임시 컬럼 생성
  #    AggregateExpression이 내부적으로 하는 작업을 미리 수행하여 ID 불일치 방지
  temp_group_by_col <- "temp_spearman_group_by"
  sobj@meta.data[[temp_group_by_col]] <- ifelse(
    grepl("^[0-9]", sobj@meta.data[[patient_col]]),
    paste0("g", sobj@meta.data[[patient_col]]),
    as.character(sobj@meta.data[[patient_col]])
  )

  # 2. Pseudobulk 발현 매트릭스 생성
  cat("환자별 Pseudobulk 프로필을 생성합니다 (임시 ID 사용)...\n")
  avg_expr <- AggregateExpression(
    sobj,
    group.by = temp_group_by_col, # 수정된 임시 컬럼 사용
    assays = assay,
    slot = layer, 
    return.seurat = FALSE
  )
  
  pseudobulk_matrix <- avg_expr[[assay]]
  
  # (참고) AggregateExpression의 'Inf' 경고 메시지:
  # 'layer = "data"' (보통 log-normalized)를 사용할 때 발생할 수 있습니다.
  # 이는 'AggregateExpression'이 내부적으로 'exp(x)-1'을 시도할 때 발생할 수 있으나,
  # 'slot = "data"'의 경우 평균(mean)을 계산하므로 결과(평균 발현 값)는 유효합니다.

  # 3. 환자별 메타데이터(ordinal_var) 준비
  meta <- sobj@meta.data
  
  # 임시 ID를 기준으로 유니크한 메타데이터 추출
  patient_meta <- meta %>%
    dplyr::select(dplyr::all_of(c(temp_group_by_col, ordinal_var))) %>%
    dplyr::distinct(!!dplyr::sym(temp_group_by_col), .keep_all = TRUE)
  
  # pseudobulk_matrix의 컬럼 이름(환자 ID) 순서대로 patient_meta 정렬
  pseudobulk_patients <- colnames(pseudobulk_matrix)
  match_indices <- match(pseudobulk_patients, patient_meta[[temp_group_by_col]])
  
  # 오류 검사 (NA가 없는지)
  if (any(is.na(match_indices))) {
      stop(paste("Patient ID 매칭 실패. 임시 ID 생성 후에도 NA가 발생:",
                 paste(pseudobulk_patients[is.na(match_indices)], collapse = ",")))
  }
  
  ordered_patient_meta <- patient_meta[match_indices, ]
  
  # (오류 수정) 원본 코드의 if문 검사를 더 강력하게 수정
  if (!all(pseudobulk_patients == ordered_patient_meta[[temp_group_by_col]])) {
    stop("환자 ID 정렬 오류: Pseudobulk 순서와 메타데이터 순서가 일치하지 않습니다.")
  }
  
  heal_values <- ordered_patient_meta[[ordinal_var]]
  
  # 4. 유전자 필터링 (Pseudobulk Matrix 기준)
  if (!is.null(top_n)) {
    cat(paste("'top_n'이 설정되었습니다.", top_n, "개 유전자를", sort_by, "기준으로 정렬합니다 (Pseudobulk 기준)...\n"))
    
    if (sort_by == "variance") {
      metrics <- apply(pseudobulk_matrix, 1, var)
    } else if (sort_by == "mean") {
      metrics <- rowMeans(pseudobulk_matrix)
    } else {
      stop("'sort_by'는 'variance' 또는 'mean'이어야 합니다.")
    }
    
    ordered_indices <- order(metrics, decreasing = decreasing)
    genes_to_test <- rownames(pseudobulk_matrix)[ordered_indices][1:min(top_n, length(metrics))]
    
    # 분석할 매트릭스 서브셋
    pseudobulk_matrix_subset <- pseudobulk_matrix[genes_to_test, ]
    
  } else {
    genes_to_test <- rownames(pseudobulk_matrix)
    pseudobulk_matrix_subset <- pseudobulk_matrix
  }

  # 5. 스피어만 상관관계 계산 (for loop 사용)
  total_genes <- length(genes_to_test)
  cat("총", total_genes, "개의 유전자에 대해 스피어만 상관관계를 계산합니다...\n")
  
  spearman_results_list <- list()
  overall_start_time <- Sys.time()
  
  for (i in seq_along(genes_to_test)) {
    gene <- genes_to_test[i]
    gene_avg_expr <- pseudobulk_matrix_subset[gene, ]
    
    # 발현에 분산이 있는지 확인
    if (var(gene_avg_expr, na.rm=TRUE) == 0) {
      spearman_results_list[[gene]] <- c(rho = NA, p.value = NA)
      next # 다음 유전자로
    }
    
    cor_test_result <- tryCatch({
      stats::cor.test(gene_avg_expr, heal_values, method = "spearman")
    }, error = function(e) { NULL })
    
    if (!is.null(cor_test_result)) {
      spearman_results_list[[gene]] <- c(
        rho = cor_test_result$estimate,
        p.value = cor_test_result$p.value
      )
    } else {
      spearman_results_list[[gene]] <- c(rho = NA, p.value = NA)
    }
    
    # 진행 상황 리포트
    if (i == 100) {
      first_100_time <- Sys.time()
      time_per_100 <- as.numeric(difftime(first_100_time, overall_start_time, units = "secs"))
      total_estimated_secs <- (time_per_100 / 100) * total_genes
      total_estimated_mins <- total_estimated_secs / 60
      
      cat(paste0("[진행 상황] 첫 100개 유전자 완료 (", round(time_per_100, 1), "초 소요)\n"))
      cat(paste0("   > 총 예상 소요 시간: 약 ", round(total_estimated_mins, 1), " 분\n"))
      
    } else if (i > 100 && i %% 100 == 0) {
      current_time <- Sys.time()
      elapsed_mins <- as.numeric(difftime(current_time, overall_start_time, units = "mins"))
      cat(paste0("[진행 상황] ", i, " / ", total_genes, " (", round(i/total_genes*100), "%) 완료. (", round(elapsed_mins, 1), "분 경과)\n"))
    }
  } # end for loop

  # 6. 결과 data.frame으로 변환
  spearman_results_df <- do.call(rbind, spearman_results_list)
  spearman_results_df <- as.data.frame(spearman_results_df)
  spearman_results_df$gene <- rownames(spearman_results_df)
  
  spearman_results_df <- spearman_results_df[!is.na(spearman_results_df$p.value), ]
  spearman_results_df <- spearman_results_df %>%
    dplyr::arrange(p.value) %>%
    dplyr::select(gene, rho, p.value)
  
  cat("스피어만 상관관계 분석 완료.\n")
  return(spearman_results_df)
}
