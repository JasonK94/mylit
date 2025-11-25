
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

# --------- Differential Expression Pipelines (MAST / MUSCAT / NEBULA) ---------

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
runMAST <- function(sobj,
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
runMAST_v1 <- function(...) {
  .Deprecated("runMAST", package = "myR")
  runMAST(...)
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
runNEBULA_v1 <- function(...) {
  .Deprecated("runNEBULA", package = "myR")
  runNEBULA(...)
}



#' @keywords internal
#' @noRd
#' Legacy wrapper retained for backwards compatibility.
runMUSCAT_v5 <- function(...) {
  .Deprecated("runMUSCAT", package = "myR")
  runMUSCAT(..., remove_na_groups = FALSE)
}

runMUSCAT2_v1 <- function(...) {
  .Deprecated("runMUSCAT", package = "myR")
  runMUSCAT(..., remove_na_groups = TRUE)
}



#' @export


# ============================================================================

#' Run MUSCAT Analysis
#'
#' @description
#' Unified MUSCAT helper that performs pseudobulking by cluster/sample and
#' then applies muscat's `pbDS` (edgeR/DESeq2/limma).  Missing-value handling
#' is controlled via `remove_na_groups`, allowing backwards-compatible behaviour
#' with the legacy v5 workflow.
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
#' @param remove_na_groups Remove cells with NA in group_id before analysis (default: TRUE)
#'
#' @return Data frame with differential expression results per cluster
#'
#' @note
#' This function automatically performs pseudobulking via
#' \code{muscat::aggregateData()}. For single-cell level analysis, use
#' \code{runNEBULA2_v1}.
#'
#' @export
runMUSCAT <- function(
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
  cluster_label_map = NULL,
  remove_na_groups = TRUE
){
  if (is.null(contrast)) stop("'contrast'를 지정하세요. 예: 'IS - SAH'")
  filter_genes <- match.arg(filter_genes)

  # deps
  req <- c("Seurat","muscat","SingleCellExperiment","SummarizedExperiment","S4Vectors","limma","dplyr")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("필요 패키지 설치: ", paste(miss, collapse=", "))

  # --- 0. NA 값 처리 (g3 등 결측치 제거) ---
  message("0/7: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  
  # 필수 컬럼 확인
  required_cols <- c(cluster_id, sample_id, group_id)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  # NA 값이 있는 셀 확인 (R의 NA와 character "NA" 모두 제거)
  if (remove_na_groups) {
    # R의 NA 값 확인
    na_mask <- is.na(meta[[group_id]]) | 
               is.na(meta[[cluster_id]]) | 
               is.na(meta[[sample_id]])
    if (!is.null(batch_id) && batch_id %in% colnames(meta)) {
      na_mask <- na_mask | is.na(meta[[batch_id]])
    }
    
    # character "NA" 문자열도 제거 (group_id, cluster_id, sample_id, batch_id)
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
  
  # 최소 그룹 수 확인
  if (length(unique(meta[[group_id]])) < 2) {
    stop(sprintf("group_id ('%s')에 최소 2개의 그룹이 필요합니다. 현재: %s",
                 group_id, paste(unique(meta[[group_id]]), collapse=", ")))
  }

  # 1) Seurat -> SCE, prepSCE
  message("1/7: Seurat -> SCE 변환 중...")
  sce <- Seurat::as.SingleCellExperiment(sobj)
  sce <- muscat::prepSCE(sce, kid = cluster_id, sid = sample_id, gid = group_id)

  # factor 보장 및 NA 제거
  message("2/7: 메타데이터 정리 중...")
  sce$cluster_id <- droplevels(factor(SummarizedExperiment::colData(sce)$cluster_id))
  sce$sample_id  <- droplevels(factor(SummarizedExperiment::colData(sce)$sample_id))
  sce$group_id   <- droplevels(factor(SummarizedExperiment::colData(sce)$group_id))
  if (!is.null(batch_id) && batch_id %in% colnames(SummarizedExperiment::colData(sce))) {
    sce[[batch_id]] <- droplevels(factor(SummarizedExperiment::colData(sce)[[batch_id]]))
  }

  # 2) Pseudobulk
  message("3/7: Pseudobulking 중...")
  pb <- muscat::aggregateData(sce, assay = "counts", by = c("cluster_id","sample_id"))

  # (선택) 특정 클러스터만
  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    pb <- pb[names(SummarizedExperiment::assays(pb)) %in% keep_clusters]
    if (length(SummarizedExperiment::assays(pb)) == 0L) stop("keep_clusters에 해당하는 클러스터가 없습니다.")
  }

  # 2-1) pb 메타 보강 (sample_id / group_id / batch)
  message("4/7: Pseudobulk 메타데이터 보강 중...")
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
  
  # NA 제거
  sce_map <- sce_map[complete.cases(sce_map), ]

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
  message("6/7: Design matrix 생성 및 DE 분석 실행 중...")
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
  message("7/7: 결과 정리 중...")
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

  message("MUSCAT2_v1 분석 완료.")
  return(combined)
}

#' @keywords internal
#' @noRd
#' Run NEBULA Analysis (v2 - Improved with missing value handling)
#'
#' @description
#' Performs differential expression analysis using NEBULA (Negative Binomial
#' mixed-effects model). NEBULA accounts for patient/sample-level random effects
#' and is suitable for multi-level experimental designs.
#' 
#' This version (v2) improves upon v1 by:
#' - Better handling of missing values (especially g3) in metadata
#' - More robust NA filtering before analysis
#' - Better error messages and validation
#'
#' @param sobj Seurat object
#' @param layer Assay layer to use (default: "counts")
#' @param fixed_effects Character vector of fixed effect variables
#'   (e.g., c("g3", "celltype_col"))
#' @param covar_effects Character vector of covariate variables
#'   (e.g., c("batch_col"))
#' @param patient_col Column name for patient/sample ID (default: "hos_no")
#' @param offset Column name for offset variable (default: "nCount_RNA")
#' @param min_count Minimum number of cells expressing a gene (default: 10)
#' @param remove_na_cells Remove cells with NA in any required variable (default: TRUE)
#'
#' @return NEBULA result object from \code{nebula::nebula()}
#'
#' @note
#' This function operates on single-cell level data. Pseudobulking can be
#' applied before calling this function using \code{runMUSCAT} or
#' \code{runNEBULA2_v1_with_pseudobulk}.
#'
._nebula_fit_with_retry <- function(grouped_data, fallback_gene_cap = 1000, min_gene_floor = 10) {
  call_nebula <- function(grouped, hl_mode = TRUE) {
    nebula::nebula(
      count = grouped$count,
      id = grouped$id,
      pred = grouped$pred,
      offset = grouped$offset,
      model = "NBLMM",
      method = if (hl_mode) "HL" else "LN"
    )
  }
  fit_once <- function(grouped) {
    tryCatch(
      call_nebula(grouped, hl_mode = TRUE),
      error = function(e) {
        warning("NEBULA with HL method failed, trying with default settings...")
        tryCatch(call_nebula(grouped, hl_mode = FALSE), error = function(e2) e2)
      }
    )
  }
  current <- grouped_data
  if (nrow(current$count) > fallback_gene_cap) {
    message(sprintf("... 유전자 수가 많아 %d개로 제한합니다.", fallback_gene_cap))
    current$count <- current$count[seq_len(fallback_gene_cap), , drop = FALSE]
  }
  repeat {
    result <- fit_once(current)
    if (!inherits(result, "error")) {
      return(result)
    }
    msg <- conditionMessage(result)
    if (!grepl("objective in x0 returns NA", msg, fixed = TRUE)) {
      stop(sprintf("NEBULA failed: %s. Try reducing number of genes or checking data quality.", msg))
    }
    new_gene_count <- max(min_gene_floor, floor(nrow(current$count) / 2))
    if (new_gene_count >= nrow(current$count)) {
      stop(sprintf("NEBULA failed: %s. Try reducing number of genes or checking data quality.", msg))
    }
    message(sprintf("... NEBULA objective returned NA. Retrying with %d genes.", new_gene_count))
    current$count <- current$count[seq_len(new_gene_count), , drop = FALSE]
  }
}

._nebula_fixed_impl <- function(sobj,
                                layer = "counts",
                                fixed_effects = c("g3"),
                                covar_effects = NULL,
                                patient_col = "hos_no",
                                offset = "nCount_RNA",
                                min_count = 10,
                                remove_na_cells = TRUE) {
  
  # --- 0. 데이터 추출 ---
  meta <- sobj@meta.data
  counts <- GetAssayData(sobj, layer = layer) # dgCMatrix (희소 행렬)
  
  # --- 1. 유전자 필터링 ---
  message(sprintf("1/7: 유전자 필터링 (min %d cells)...", min_count))
  keep_genes <- rowSums(counts > 0) >= min_count
  counts_filtered <- counts[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(counts)))
  
  # --- 2. NA 값 확인 및 제거 ---
  # NEBULA는 'covar_effects'도 고정 효과로 처리해야 함
  all_fixed_vars <- c(fixed_effects, covar_effects)
  all_fixed_vars <- all_fixed_vars[!is.null(all_fixed_vars)]
  
  # NA 체크할 모든 변수 목록
  vars_to_check <- c(all_fixed_vars, patient_col, offset)
  vars_to_check <- vars_to_check[vars_to_check %in% colnames(meta)]
  
  message("2/7: 모델 변수에서 NA 값 확인 중...")
  message(paste("... 확인 대상:", paste(vars_to_check, collapse = ", ")))
  
  # 각 변수별 NA 개수 확인
  if (remove_na_cells) {
    na_counts <- vapply(vars_to_check, function(v) sum(is.na(meta[[v]])), integer(1))
    if (any(na_counts > 0)) {
      message("... 변수별 NA 개수:")
      for (v in names(na_counts[na_counts > 0])) {
        message(sprintf("    %s: %d 개", v, na_counts[v]))
      }
    }
    
    keep_cells_idx <- complete.cases(meta[, vars_to_check])
    n_removed <- sum(!keep_cells_idx)
    
    if (n_removed > 0) {
      message(sprintf("... NA 값으로 인해 %d 개의 세포를 제거합니다.", n_removed))
    } else {
      message("... NA 값이 없습니다.")
    }
  } else {
    keep_cells_idx <- rep(TRUE, nrow(meta))
    n_removed <- 0
    # NA가 있으면 경고
    if (any(!complete.cases(meta[, vars_to_check]))) {
      warning("일부 변수에 NA가 있지만 제거하지 않습니다. 분석 결과에 영향을 줄 수 있습니다.")
    }
  }
  
  # NA가 없는 '깨끗한' 데이터 생성
  meta_clean <- meta[keep_cells_idx, ]
  counts_clean <- counts_filtered[, keep_cells_idx]
  
  message(sprintf("... 최종 분석 대상 세포: %d 개", nrow(meta_clean)))
  
  # --- 3. 디자인 행렬 및 벡터 생성 ---
  message("3/7: 디자인 행렬 및 벡터 생성 중...")
  
  # 범주형 변수와 숫자형 변수 분리
  # 'offset'은 숫자로 유지해야 함!
  factor_vars <- c(all_fixed_vars, patient_col)
  factor_vars <- factor_vars[factor_vars %in% colnames(meta_clean)]
  
  # factor 변환
  meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], as.factor)
  
  # 각 factor 변수의 레벨 수 확인
  for (v in factor_vars) {
    n_levels <- length(levels(meta_clean[[v]]))
    message(sprintf("... %s: %d 레벨 (%s)", v, n_levels, 
                    paste(levels(meta_clean[[v]]), collapse=", ")))
  }
  
  # 'offset'은 numeric으로 강제 변환 (안전장치)
  if (offset %in% colnames(meta_clean)) {
    meta_clean[[offset]] <- as.numeric(meta_clean[[offset]])
    if (any(is.na(meta_clean[[offset]]))) {
        stop(sprintf("'%s' 컬럼에 NA가 아닌 숫자형 값만 있어야 합니다.", offset))
    }
    if (any(meta_clean[[offset]] <= 0)) {
      warning(sprintf("'%s' 컬럼에 0 이하 값이 있습니다. offset은 양수여야 합니다.", offset))
    }
  } else {
    stop(sprintf("'%s' 컬럼이 메타데이터에 없습니다.", offset))
  }
  
  # formula 생성
  formula_str <- paste("~", paste(all_fixed_vars, collapse = " + "))
  message(sprintf("... 사용할 고정 효과 포뮬러: %s", formula_str))
  
  # 설계 행렬 생성 및 특이성(singularity) 확인
  design_matrix <- model.matrix(as.formula(formula_str), data = meta_clean)
  zero_var_cols <- which(apply(design_matrix, 2, function(x) length(unique(x)) <= 1))
  intercept_idx <- which(colnames(design_matrix) == "(Intercept)")
  zero_var_cols <- setdiff(zero_var_cols, intercept_idx)
  if (length(zero_var_cols) > 0) {
    message("... zero-variance predictors removed: ",
            paste(colnames(design_matrix)[zero_var_cols], collapse = ", "))
    design_matrix <- design_matrix[, -zero_var_cols, drop = FALSE]
  }
  
  # 설계 행렬의 특이성 확인 (rank deficiency) - qr() 사용
  design_qr <- qr(design_matrix)
  design_rank <- design_qr$rank
  n_cols <- ncol(design_matrix)
  
  # 완전 분리된 조합 확인 (설계 행렬이 특이하기 전에 미리 확인)
  if (length(all_fixed_vars) >= 2) {
    message("... 변수 간 완전 분리(complete separation) 확인 중...")
    separation_issues <- FALSE
    for (i in 1:(length(all_fixed_vars)-1)) {
      for (j in (i+1):length(all_fixed_vars)) {
        var1 <- all_fixed_vars[i]
        var2 <- all_fixed_vars[j]
        contingency <- table(meta_clean[[var1]], meta_clean[[var2]])
        zero_cells <- sum(contingency == 0)
        if (zero_cells > 0) {
          separation_issues <- TRUE
          warning(sprintf("%s와 %s 사이에 완전 분리된 조합이 있습니다 (0인 셀: %d개).", 
                         var1, var2, zero_cells))
          message(sprintf("  %s x %s contingency table:", var1, var2))
          print(contingency)
        }
      }
    }
    if (!separation_issues) {
      message("... 완전 분리 문제 없음")
    }
  }
  
  if (design_rank < n_cols) {
    warning(sprintf("설계 행렬이 특이(singular)합니다. rank=%d < columns=%d. 완전 분리(complete separation) 문제로 인해 분석이 실패할 수 있습니다.", 
                    design_rank, n_cols))
    dependent_idx <- design_qr$pivot[seq.int(design_rank + 1, n_cols)]
    drop_cols <- colnames(design_matrix)[dependent_idx]
    message("... Removing linearly dependent predictors: ", paste(drop_cols, collapse = ", "))
    keep_idx <- design_qr$pivot[seq_len(design_rank)]
    design_matrix <- design_matrix[, keep_idx, drop = FALSE]
    design_qr <- qr(design_matrix)
    design_rank <- design_qr$rank
    n_cols <- ncol(design_matrix)
    message("설계 행렬이 특이하지만 분석을 계속 진행합니다. 오류가 발생하면:")
    message("  1. covar_effects를 제거하거나 다른 변수를 사용하세요.")
    message("  2. 완전 분리된 조합을 제거하거나 데이터를 필터링하세요.")
    message("  3. min_count를 높여서 더 적은 유전자로 분석하세요.")
  }
  
  message(sprintf("... 설계 행렬: %d 행 x %d 열 (rank=%d)", 
                  nrow(design_matrix), ncol(design_matrix), design_rank))
  
  # id 및 offset 벡터 추출
  id_vector <- meta_clean[[patient_col]]
  offset_vector <- meta_clean[[offset]]
  
  # --- 4. group_cell()로 데이터 정렬 ---
  message(sprintf("4/7: NEBULA를 위해 id (%s) 기준으로 데이터 정렬 중...", patient_col))
  data_grouped <- nebula::group_cell(
    count = counts_clean,
    id = id_vector,
    pred = design_matrix,
    offset = offset_vector
  )
  
  message(sprintf("... %d 개의 유전자, %d 개의 샘플", 
                  nrow(data_grouped$count), length(unique(data_grouped$id))))
  
  # --- 5. NEBULA 실행 ---
  message("5/7: NEBULA 실행 중 (NBLMM)...")
  
  # 유전자 수가 너무 많으면 일부만 테스트
  n_genes <- nrow(data_grouped$count)
  n_samples <- length(unique(data_grouped$id))
  
  message(sprintf("... %d 개의 유전자, %d 개의 샘플 분석 시작", n_genes, n_samples))
  
  # 샘플 수가 적거나 설계 행렬이 특이하면 유전자 수 제한
  if (n_samples < 10 || design_rank < n_cols) {
    if (n_genes > 1000) {
      message(sprintf("... 샘플 수가 적거나 설계 행렬이 특이하여 유전자 수를 1000개로 제한합니다."))
      gene_subset <- sample(1:n_genes, min(1000, n_genes))
      data_grouped$count <- data_grouped$count[gene_subset, , drop = FALSE]
      message(sprintf("... %d 개의 유전자로 분석 진행", nrow(data_grouped$count)))
    }
  }
  
  re_nebula <- ._nebula_fit_with_retry(data_grouped)
  
  # --- 6. 결과 정리 ---
  message("6/7: 결과 정리 중...")
  # NEBULA 결과는 리스트로 반환됨
  # summary는 자동으로 포함되어 있음
  
  # --- 7. 완료 ---
  message("7/7: 분석 완료.")
  
  return(re_nebula)
}


#' Run NEBULA differential expression analysis
#'
#' @description
#' Unified interface for NEBULA-based differential expression.  Use `formula` to
#' activate the formula parser or `pseudobulk`/`pseudobulk_args` to trigger the
#' pseudobulk workflow.
#'
#' @param sobj Seurat object
#' @param layer Assay layer to use (default: "counts")
#' @param fixed_effects Character vector of fixed effects used in the classic
#'   interface.
#' @param covar_effects Character vector of covariates.
#' @param patient_col Column name for patient/sample ID.
#' @param offset Column name for offset variable.
#' @param min_count Minimum number of cells expressing a gene.
#' @param remove_na_cells Remove cells with NA in model variables.
#' @param formula Optional lme4-style formula string or object.
#' @param pseudobulk Logical or list to enable the pseudobulk workflow.
#' @param pseudobulk_args Named list overriding pseudobulk defaults
#'   (e.g., cluster/sample identifiers, offsets, filtering thresholds).
#' @param separation_min_cells Minimum cells required for each factor level when
#'   using the formula mode. Levels below the threshold are pruned to avoid
#'   complete separation (default: 20).
#' @param simplify_interactions When TRUE (default) problematic interaction
#'   terms detected during formula analysis are dropped automatically to avoid
#'   singular design matrices.
#' @param separation_group_var Optional metadata column used to enforce gene-level
#'   balance (e.g., treatment group). When supplied, genes must have non-zero
#'   expression in at least \code{separation_min_cells_per_group} cells for every
#'   level of this column.
#' @param separation_min_cells_per_group Minimum number of cells per level for
#'   the gene-level balance filter (default: 3).
#'
#' @return NEBULA result object. When pseudobulk mode is enabled a list is
#' returned containing the NEBULA result plus the pseudobulk counts/metadata.
#'
#' @export
runNEBULA <- function(
  sobj,
  layer = "counts",
  fixed_effects = c("g3"),
  covar_effects = NULL,
  patient_col = "hos_no",
  offset = "nCount_RNA",
  min_count = 10,
  remove_na_cells = TRUE,
  formula = NULL,
  pseudobulk = FALSE,
  pseudobulk_args = list(),
  separation_min_cells = 20,
  simplify_interactions = TRUE,
  separation_group_var = NULL,
  separation_min_cells_per_group = 3
) {
  if (!is.null(formula)) {
    return(._nebula_formula_impl(
      sobj = sobj,
      formula = formula,
      layer = layer,
      patient_col = patient_col,
      offset = offset,
      min_count = min_count,
      remove_na_cells = remove_na_cells,
      separation_min_cells = separation_min_cells,
      simplify_interactions = simplify_interactions,
      separation_group_var = separation_group_var,
      separation_min_cells_per_group = separation_min_cells_per_group
    ))
  }
  use_pb <- isTRUE(pseudobulk) ||
    (is.list(pseudobulk) && length(pseudobulk) > 0) ||
    length(pseudobulk_args) > 0
  if (use_pb) {
    pb_defaults <- list(
      layer = layer,
      cluster_id = "seurat_clusters",
      sample_id = "hos_no",
      group_id = "type",
      fixed_effects = fixed_effects,
      covar_effects = covar_effects,
      patient_col = patient_col,
      offset_method = "sum",
      min_count = min_count,
      min_cells_per_pb = 3,
      remove_na_cells = remove_na_cells,
      keep_clusters = NULL
    )
    pb_override <- if (is.list(pseudobulk)) pseudobulk else list()
    pb_args <- utils::modifyList(pb_defaults,
                                 utils::modifyList(pb_override, pseudobulk_args))
    return(do.call(._nebula_pseudobulk_impl, c(list(sobj = sobj), pb_args)))
  }
  ._nebula_fixed_impl(
    sobj = sobj,
    layer = layer,
    fixed_effects = fixed_effects,
    covar_effects = covar_effects,
    patient_col = patient_col,
    offset = offset,
    min_count = min_count,
    remove_na_cells = remove_na_cells
  )
}

runNEBULA2_v1 <- function(...) {
  .Deprecated("runNEBULA", package = "myR")
  runNEBULA(...)
}

runNEBULA2_v1_with_formula <- function(sobj,
                                       formula,
                                       layer = "counts",
                                       patient_col = NULL,
                                       offset = "nCount_RNA",
                                       min_count = 10,
                                       remove_na_cells = TRUE) {
  .Deprecated("runNEBULA", package = "myR")
  runNEBULA(
    sobj = sobj,
    layer = layer,
    patient_col = if (is.null(patient_col)) "hos_no" else patient_col,
    offset = offset,
    min_count = min_count,
    remove_na_cells = remove_na_cells,
    formula = formula
  )
}

runNEBULA2_v1_with_pseudobulk <- function(sobj,
                                          layer = "counts",
                                          cluster_id = "seurat_clusters",
                                          sample_id  = "hos_no",
                                          group_id   = "type",
                                          fixed_effects = c("g3"),
                                          covar_effects = NULL,
                                          patient_col = "hos_no",
                                          offset_method = c("sum", "mean", "n_cells"),
                                          min_count = 10,
                                          min_cells_per_pb = 3,
                                          remove_na_cells = TRUE,
                                          keep_clusters = NULL) {
  .Deprecated("runNEBULA", package = "myR")
  runNEBULA(
    sobj = sobj,
    layer = layer,
    fixed_effects = fixed_effects,
    covar_effects = covar_effects,
    patient_col = patient_col,
    min_count = min_count,
    remove_na_cells = remove_na_cells,
    pseudobulk = list(
      cluster_id = cluster_id,
      sample_id = sample_id,
      group_id = group_id,
      offset_method = match.arg(offset_method),
      min_cells_per_pb = min_cells_per_pb,
      keep_clusters = keep_clusters
    )
  )
}


#' @keywords internal
#' @noRd
#' Run NEBULA Analysis with Pseudobulk (v2)
#'
#' @description
#' Performs NEBULA analysis after pseudobulking by cluster and sample.
#' This function first aggregates cells within each cluster and sample to create
#' pseudobulk samples, then runs NEBULA analysis on the pseudobulked data.
#' 
#' This allows combining the benefits of pseudobulking (reduced computational
#' burden, better signal-to-noise ratio) with NEBULA's mixed-effects modeling
#' (accounting for patient-level random effects).
#'
#' @param sobj Seurat object
#' @param layer Assay layer to use (default: "counts")
#' @param cluster_id Column name for cell type/cluster (default: "seurat_clusters")
#' @param sample_id Column name for sample ID (default: "hos_no")
#' @param group_id Column name for group/condition (default: "type")
#' @param fixed_effects Character vector of fixed effect variables
#'   (e.g., c("g3", "cluster_id"))
#' @param covar_effects Character vector of covariate variables
#'   (e.g., c("batch_col"))
#' @param patient_col Column name for patient/sample ID (default: "hos_no")
#' @param offset_method Method for calculating offset: "sum" (sum of counts),
#'   "mean" (mean of counts), or "n_cells" (number of cells) (default: "sum")
#' @param min_count Minimum number of cells expressing a gene (default: 10)
#' @param min_cells_per_pb Minimum cells per pseudobulk sample (default: 3)
#' @param remove_na_cells Remove cells with NA in any required variable (default: TRUE)
#' @param keep_clusters Optional vector of cluster IDs to keep
#'
#' @return List with:
#'   \itemize{
#'     \item nebula_result: NEBULA result object
#'     \item pseudobulk_meta: Metadata for pseudobulk samples
#'     \item pseudobulk_counts: Pseudobulk count matrix
#'   }
#'
#' @note
#' This function performs pseudobulking by cluster and sample, then runs NEBULA.
#' The offset is calculated from the pseudobulk counts (sum by default).
#'
._nebula_pseudobulk_impl <- function(sobj,
                                layer = "counts",
                                cluster_id = "seurat_clusters",
                                sample_id  = "hos_no",
                                group_id   = "type",
                                fixed_effects = c("g3"),
                                covar_effects = NULL,
                                patient_col = "hos_no",
                                offset_method = c("sum", "mean", "n_cells"),
                                min_count = 10,
                                min_cells_per_pb = 3,
                                remove_na_cells = TRUE,
                                keep_clusters = NULL) {
  
  offset_method <- match.arg(offset_method)
  
  # --- 0. 데이터 추출 및 NA 처리 ---
  message("0/8: 메타데이터에서 NA 값 확인 중...")
  meta <- sobj@meta.data
  
  # 필수 컬럼 확인
  required_cols <- c(cluster_id, sample_id, group_id, patient_col)
  missing_cols <- required_cols[!required_cols %in% colnames(meta)]
  if (length(missing_cols) > 0) {
    stop(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
  }
  
  # NA 값이 있는 셀 확인
  if (remove_na_cells) {
    na_mask <- is.na(meta[[group_id]]) | 
               is.na(meta[[cluster_id]]) | 
               is.na(meta[[sample_id]]) |
               is.na(meta[[patient_col]])
    
    # fixed_effects와 covar_effects도 확인
    all_vars <- c(fixed_effects, covar_effects)
    all_vars <- all_vars[!is.null(all_vars) & all_vars %in% colnames(meta)]
    for (v in all_vars) {
      na_mask <- na_mask | is.na(meta[[v]])
    }
    
    n_na_cells <- sum(na_mask)
    if (n_na_cells > 0) {
      message(sprintf("... NA 값이 있는 %d 개의 세포를 제거합니다.", n_na_cells))
      sobj <- sobj[, !na_mask]
      meta <- sobj@meta.data
    } else {
      message("... NA 값이 없습니다.")
    }
  }
  

  # --- 1. Pseudobulking ---
  message("1/8: Pseudobulking 중...")
  
  # 특정 클러스터만 필터링
  if (!is.null(keep_clusters)) {
    keep_clusters <- as.character(keep_clusters)
    cells_to_keep <- rownames(meta)[meta[[cluster_id]] %in% keep_clusters]
    sobj <- sobj[, cells_to_keep]
    meta <- sobj@meta.data
    message(sprintf("... 클러스터 %s만 사용 (%d 세포)", 
                    paste(keep_clusters, collapse=", "), ncol(sobj)))
  }
  
  # Pseudobulk 생성: cluster_id와 sample_id로 집계
  counts <- GetAssayData(sobj, layer = layer)
  
  # 각 cluster-sample 조합에 대해 집계
  pb_list <- list()
  pb_meta_list <- list()
  
  clusters <- unique(meta[[cluster_id]])
  samples <- unique(meta[[sample_id]])
  
  for (clust in clusters) {
    for (samp in samples) {
      # 해당 cluster-sample 조합의 세포 선택
      cell_mask <- meta[[cluster_id]] == clust & meta[[sample_id]] == samp
      n_cells <- sum(cell_mask)
      
      if (n_cells >= min_cells_per_pb) {
        # Pseudobulk ID 생성
        pb_id <- paste(clust, samp, sep = "_")
        
        # Counts 합계
        if (n_cells == 1) {
          pb_counts <- counts[, cell_mask, drop = FALSE]
          colnames(pb_counts) <- pb_id
        } else {
          pb_counts <- Matrix::rowSums(counts[, cell_mask, drop = FALSE])
          pb_counts <- as.matrix(pb_counts)
          colnames(pb_counts) <- pb_id
        }
        
        pb_list[[pb_id]] <- pb_counts
        
        # Metadata 생성
        pb_meta_row <- data.frame(
          pb_id = pb_id,
          cluster_id = clust,
          sample_id = samp,
          patient_id = unique(meta[cell_mask, patient_col])[1],
          n_cells = n_cells,
          stringsAsFactors = FALSE
        )
        
        # group_id 추가
        pb_meta_row[[group_id]] <- unique(meta[cell_mask, group_id])[1]
        
        # fixed_effects와 covar_effects 추가 (sample-level이면 첫 번째 값 사용)
        all_vars <- c(fixed_effects, covar_effects)
        all_vars <- all_vars[!is.null(all_vars) & all_vars %in% colnames(meta)]
        for (v in all_vars) {
          pb_meta_row[[v]] <- unique(meta[cell_mask, v])[1]
        }
        
        pb_meta_list[[pb_id]] <- pb_meta_row
      }
    }
  }
  
  if (length(pb_list) == 0) {
    stop("Pseudobulk 샘플이 생성되지 않았습니다. min_cells_per_pb를 낮추거나 데이터를 확인하세요.")
  }
  
  # Pseudobulk count matrix 생성
  message("2/8: Pseudobulk count matrix 생성 중...")
  pb_counts <- do.call(cbind, pb_list)
  pb_meta <- do.call(rbind, pb_meta_list)
  rownames(pb_meta) <- pb_meta$pb_id
  
  # 컬럼 순서 맞추기
  pb_counts <- pb_counts[, pb_meta$pb_id, drop = FALSE]
  
  message(sprintf("... %d 개의 pseudobulk 샘플 생성 (%d 클러스터, %d 샘플)",
                  ncol(pb_counts), length(clusters), length(samples)))
  
  # --- 3. Offset 계산 ---
  message("3/8: Offset 계산 중...")
  if (offset_method == "sum") {
    pb_meta$offset <- colSums(pb_counts)
  } else if (offset_method == "mean") {
    pb_meta$offset <- colMeans(pb_counts)
  } else if (offset_method == "n_cells") {
    pb_meta$offset <- pb_meta$n_cells
  }
  
  message(sprintf("... Offset 방법: %s (범위: %.2f - %.2f)",
                  offset_method, min(pb_meta$offset), max(pb_meta$offset)))
  
  # --- 4. 유전자 필터링 ---
  message("4/8: 유전자 필터링 중...")
  keep_genes <- rowSums(pb_counts > 0) >= min_count
  pb_counts_filtered <- pb_counts[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(pb_counts)))
  
  # --- 5. 디자인 행렬 생성 ---
  message("5/8: 디자인 행렬 생성 중...")
  
  # fixed_effects에 cluster_id가 없으면 추가
  if (!cluster_id %in% fixed_effects) {
    fixed_effects <- c(fixed_effects, "cluster_id")
    message(sprintf("... cluster_id를 fixed_effects에 추가"))
  }
  
  all_fixed_vars <- c(fixed_effects, covar_effects)
  all_fixed_vars <- all_fixed_vars[!is.null(all_fixed_vars)]
  all_fixed_vars <- all_fixed_vars[all_fixed_vars %in% colnames(pb_meta)]
  
  # factor 변환
  factor_vars <- c(all_fixed_vars, "patient_id")
  pb_meta[factor_vars] <- lapply(pb_meta[factor_vars], as.factor)
  
  # formula 생성
  formula_str <- paste("~", paste(all_fixed_vars, collapse = " + "))
  message(sprintf("... 사용할 고정 효과 포뮬러: %s", formula_str))
  
  design_matrix <- model.matrix(as.formula(formula_str), data = pb_meta)
  
  # id 및 offset 벡터 추출
  id_vector <- pb_meta[["patient_id"]]
  offset_vector <- pb_meta[["offset"]]
  
  # --- 6. group_cell()로 데이터 정렬 ---
  message("6/8: NEBULA를 위해 id 기준으로 데이터 정렬 중...")
  data_grouped <- nebula::group_cell(
    count = pb_counts_filtered,
    id = id_vector,
    pred = design_matrix,
    offset = offset_vector
  )
  
  message(sprintf("... %d 개의 유전자, %d 개의 샘플", 
                  nrow(data_grouped$count), length(unique(data_grouped$id))))
  
  # --- 7. NEBULA 실행 ---
  message("7/8: NEBULA 실행 중 (NBLMM)...")
  re_nebula <- tryCatch({
    nebula::nebula(
      count = data_grouped$count,
      id = data_grouped$id,
      pred = data_grouped$pred,
      offset = data_grouped$offset,
      model = "NBLMM",
      method = "HL"
    )
  }, error = function(e) {
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
  
  # --- 8. 결과 반환 ---
  message("8/8: 분석 완료.")
  
  return(list(
    nebula_result = re_nebula,
    pseudobulk_meta = pb_meta,
    pseudobulk_counts = pb_counts_filtered
  ))
}

#' @keywords internal
#' @noRd
#' Remove sparse factor levels to mitigate complete separation
._nebula_prune_sparse_levels <- function(meta,
                                         factor_vars,
                                         min_cells_per_level = 10,
                                         max_passes = 3) {
  if (!length(factor_vars) || min_cells_per_level <= 1) {
    return(list(
      meta = meta,
      keep_idx = rep(TRUE, nrow(meta)),
      removed_cells = 0L,
      passes = 0L,
      removal_log = list()
    ))
  }
  keep_idx <- rep(TRUE, nrow(meta))
  removal_log <- list()
  passes <- 0L
  repeat {
    meta_view <- meta[keep_idx, , drop = FALSE]
    if (!nrow(meta_view)) {
      break
    }
    rows_to_drop <- rep(FALSE, nrow(meta_view))
    for (var in factor_vars) {
      if (!var %in% colnames(meta_view)) {
        next
      }
      freq <- table(meta_view[[var]], useNA = "no")
      low_levels <- names(freq)[freq < min_cells_per_level]
      if (!length(low_levels)) {
        next
      }
      removal_log[[var]] <- sort(unique(c(removal_log[[var]], low_levels)))
      rows_to_drop <- rows_to_drop | meta_view[[var]] %in% low_levels
    }
    if (!any(rows_to_drop) || passes >= max_passes) {
      break
    }
    orig_idx <- which(keep_idx)[rows_to_drop]
    keep_idx[orig_idx] <- FALSE
    passes <- passes + 1L
  }
  list(
    meta = meta[keep_idx, , drop = FALSE],
    keep_idx = keep_idx,
    removed_cells = sum(!keep_idx),
    passes = passes,
    removal_log = removal_log
  )
}


#' @keywords internal
#' @noRd
#' Detect complete separation issues between factor pairs
._nebula_detect_complete_separation <- function(meta, vars) {
  result <- list(
    has_issue = FALSE,
    pairs = list(),
    problem_vars = character(0)
  )
  if (length(vars) < 2) {
    return(result)
  }
  combs <- utils::combn(vars, 2, simplify = FALSE)
  for (pair in combs) {
    var1 <- pair[1]
    var2 <- pair[2]
    if (!var1 %in% colnames(meta) || !var2 %in% colnames(meta)) {
      next
    }
    contingency <- table(meta[[var1]], meta[[var2]])
    zero_cells <- sum(contingency == 0)
    if (zero_cells > 0) {
      result$has_issue <- TRUE
      result$pairs[[length(result$pairs) + 1]] <- list(
        var1 = var1,
        var2 = var2,
        zero_cells = zero_cells,
        table = contingency
      )
    }
  }
  if (length(result$pairs)) {
    result$problem_vars <- unique(unlist(lapply(result$pairs, function(x) c(x$var1, x$var2))))
  }
  result
}


#' @keywords internal
#' @noRd
#' Emit warnings/messages for separation diagnostics
._nebula_log_separation <- function(separation_info) {
  if (!is.list(separation_info) || !isTRUE(separation_info$has_issue)) {
    message("... 완전 분리 문제 없음")
    return(invisible(NULL))
  }
  for (entry in separation_info$pairs) {
    warning(sprintf("%s와 %s 사이에 완전 분리된 조합이 있습니다 (0인 셀: %d개).",
                    entry$var1, entry$var2, entry$zero_cells))
    message(sprintf("  %s x %s contingency table:", entry$var1, entry$var2))
    print(entry$table)
  }
  invisible(NULL)
}


#' @keywords internal
#' @noRd
#' Remove factor levels that lack group balance
._nebula_prune_unbalanced_levels <- function(meta,
                                             factor_vars,
                                             group_var,
                                             min_cells_per_group = 1,
                                             max_passes = 3) {
  if (is.null(group_var) ||
      !group_var %in% colnames(meta) ||
      min_cells_per_group <= 0 ||
      !length(factor_vars)) {
    return(list(
      meta = meta,
      keep_idx = rep(TRUE, nrow(meta)),
      removed_cells = 0L,
      passes = 0L,
      removal_log = list()
    ))
  }
  group_levels <- levels(droplevels(as.factor(meta[[group_var]])))
  if (length(group_levels) < 2) {
    return(list(
      meta = meta,
      keep_idx = rep(TRUE, nrow(meta)),
      removed_cells = 0L,
      passes = 0L,
      removal_log = list()
    ))
  }
  keep_idx <- rep(TRUE, nrow(meta))
  removal_log <- list()
  passes <- 0L
  repeat {
    meta_view <- meta[keep_idx, , drop = FALSE]
    if (!nrow(meta_view)) {
      break
    }
    group_view <- droplevels(as.factor(meta_view[[group_var]]))
    if (length(levels(group_view)) < 2) {
      break
    }
    rows_to_drop <- rep(FALSE, nrow(meta_view))
    for (var in factor_vars) {
      if (!var %in% colnames(meta_view)) {
        next
      }
      tbl <- table(meta_view[[var]], group_view)
      if (!nrow(tbl)) {
        next
      }
      low_levels <- rownames(tbl)[apply(tbl < min_cells_per_group, 1, any)]
      if (!length(low_levels)) {
        next
      }
      removal_log[[var]] <- sort(unique(c(removal_log[[var]], low_levels)))
      rows_to_drop <- rows_to_drop | meta_view[[var]] %in% low_levels
    }
    if (!any(rows_to_drop) || passes >= max_passes) {
      break
    }
    orig_idx <- which(keep_idx)[rows_to_drop]
    keep_idx[orig_idx] <- FALSE
    passes <- passes + 1L
  }
  list(
    meta = meta[keep_idx, , drop = FALSE],
    keep_idx = keep_idx,
    removed_cells = sum(!keep_idx),
    passes = passes,
    removal_log = removal_log
  )
}

# Formula 1: ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/patient)

#' @keywords internal
#' @noRd
#' Run NEBULA Analysis with Formula (v2)
#'
#' @description
#' Performs differential expression analysis using NEBULA with a formula interface.
#' Supports lme4-style formulas including nested random effects.
#' Note: NEBULA only supports single-level random effects, so nested random effects
#' like (1|GEM/patient) will be converted to (1|patient) with GEM as a fixed effect.
#'
#' @param sobj Seurat object
#' @param formula Character string or formula object (e.g., "~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/patient)")
#' @param layer Assay layer to use (default: "counts")
#' @param patient_col Column name for patient/sample ID. If NULL, extracted from formula's random effects
#' @param offset Column name for offset variable (default: "nCount_RNA")
#' @param min_count Minimum number of cells expressing a gene (default: 10)
#' @param remove_na_cells Remove cells with NA in model variables (default: TRUE)
#' @param separation_min_cells Minimum cells required for each factor level.
#'   Levels below the threshold are pruned iteratively to reduce complete
#'   separation (default: 20).
#' @param simplify_interactions Automatically drop interaction terms that involve
#'   variables flagged during separation diagnostics (default: TRUE).
#' @param separation_group_var Optional metadata column name used to enforce gene
#'   balance across groups (e.g., "g3"). When provided, each retained gene must
#'   have non-zero expression in at least
#'   \code{separation_min_cells_per_group} cells for every level of this column.
#' @param separation_min_cells_per_group Minimum number of non-zero cells per
#'   level for the gene-level group balance filter (default: 3).
#'
#' @return NEBULA result object from \code{nebula::nebula()}
#'
._nebula_formula_impl <- function(sobj,
                                       formula,
                                       layer = "counts",
                                       patient_col = NULL,
                                       offset = "nCount_RNA",
                                       min_count = 10,
                                       remove_na_cells = TRUE,
                                       separation_min_cells = 20,
                                       simplify_interactions = TRUE,
                                       separation_group_var = NULL,
                                       separation_min_cells_per_group = 3) {
  
  # --- 0. Formula 파싱 ---
  if (is.character(formula)) {
    formula_obj <- as.formula(formula)
  } else if (inherits(formula, "formula")) {
    formula_obj <- formula
  } else {
    stop("'formula'는 문자열 또는 formula 객체여야 합니다.")
  }
  
  message("0/8: Formula 파싱 중...")
  formula_deparse <- deparse(formula_obj)
  formula_str <- paste(formula_deparse, collapse = " ")  # 여러 줄을 하나로 합치기
  message(sprintf("... Formula: %s", formula_str))
  
  # Random effects 추출: (1|...) 또는 (1 | ...) 패턴 찾기 (공백 허용)
  random_pattern <- "\\(1\\s*\\|[^)]+\\)"
  random_matches <- regmatches(formula_str, gregexpr(random_pattern, formula_str))[[1]]
  
  if (length(random_matches) == 0) {
    stop("Formula에 random effects가 없습니다. 예: (1|patient) 또는 (1|GEM/patient)")
  }
  
  # Nested random effects 처리: (1|GEM/patient) -> patient만 추출
  # random_matches는 "(1 | GEM/hos_no)" 형태
  # 괄호와 "1 |" 제거하여 실제 변수만 추출
  random_terms <- gsub("^\\(1\\s*\\|\\s*", "", random_matches)  # 앞부분 제거
  random_terms <- gsub("\\)$", "", random_terms)  # 뒷부분 괄호 제거
  
  # nested 처리: GEM/patient -> patient만 사용
  nested_terms <- strsplit(random_terms, "/")
  actual_random_terms <- sapply(nested_terms, function(x) {
    if (length(x) > 1) {
      # nested인 경우: 마지막 요소 (patient)만 사용
      # 앞의 요소들 (GEM)은 fixed effect로 추가
      message(sprintf("... Nested random effect 발견: %s", paste(x, collapse="/")))
      message(sprintf("   → Random effect로 사용: %s", x[length(x)]))
      if (length(x) > 2) {
        warning("3단계 이상의 nested random effects는 지원되지 않습니다. 마지막 레벨만 사용합니다.")
      }
      x[length(x)]
    } else {
      x[1]
    }
  })
  
  # patient_col 추출
  if (is.null(patient_col)) {
    patient_col <- actual_random_terms[1]  # 첫 번째 random effect를 patient로 사용
    message(sprintf("... patient_col 자동 설정: %s", patient_col))
  }
  
  # Nested 구조에서 상위 레벨들 (예: GEM)을 fixed effect로 추가
  fixed_effects_from_nested <- character(0)
  for (i in seq_along(nested_terms)) {
    if (length(nested_terms[[i]]) > 1) {
      # GEM/patient -> GEM을 fixed effect로 추가
      upper_levels <- nested_terms[[i]][-length(nested_terms[[i]])]
      fixed_effects_from_nested <- c(fixed_effects_from_nested, upper_levels)
    }
  }
  
  # Random effects 제거한 formula에서 fixed effects 추출
  formula_no_random <- formula_str
  for (rm in random_matches) {
    # 정규식 이스케이프 처리하여 정확히 제거
    rm_escaped <- gsub("([()|])", "\\\\\\1", rm)  # 특수문자 이스케이프
    formula_no_random <- gsub(rm_escaped, "", formula_no_random)
  }
  formula_no_random <- gsub("\\+\\s*\\+", "+", formula_no_random)  # 연속된 + 제거
  formula_no_random <- gsub("~\\s*\\+", "~", formula_no_random)    # ~ 뒤의 + 제거
  formula_no_random <- gsub("\\s*\\+\\s*$", "", formula_no_random) # 끝의 + 제거
  formula_no_random <- gsub("\\s+", " ", formula_no_random)  # 연속된 공백 제거
  formula_no_random <- trimws(formula_no_random)
  
  # 빈 formula 확인
  if (formula_no_random == "~" || nchar(trimws(gsub("~", "", formula_no_random))) == 0) {
    stop("Random effects 제거 후 fixed effects가 없습니다.")
  }
  
  # Formula를 파싱하여 terms 추출
  tryCatch({
    temp_formula <- as.formula(formula_no_random)
    formula_terms <- attr(terms(temp_formula), "term.labels")
  }, error = function(e) {
    stop(sprintf("Formula 파싱 실패: %s", conditionMessage(e)))
  })
  
  # 고정 효과 목록 (중복 제거)
  fixed_effects_all <- unique(c(fixed_effects_from_nested, formula_terms))
  
  # 교호작용 항 처리 (g3:anno3.scvi -> 그대로 유지)
  # main effects 추출 (교호작용에 포함된 변수)
  main_effects_from_interactions <- character(0)
  interaction_terms <- fixed_effects_all[grepl(":", fixed_effects_all)]
  for (it in interaction_terms) {
    main_effects_from_interactions <- c(main_effects_from_interactions, 
                                       strsplit(it, ":")[[1]])
  }
  
  # 모든 변수 추출 (메타데이터 컬럼명과 매칭)
  all_vars <- unique(c(
    setdiff(fixed_effects_all, interaction_terms),  # 교호작용 제외한 main effects
    main_effects_from_interactions,  # 교호작용에서 추출한 main effects
    patient_col
  ))
  
  message(sprintf("... Fixed effects: %s", paste(fixed_effects_all, collapse=", ")))
  message(sprintf("... Random effects: %s", paste(actual_random_terms, collapse=", ")))
  message(sprintf("... Patient column: %s", patient_col))
  
  # --- 1. 데이터 추출 ---
  meta <- sobj@meta.data
  counts <- GetAssayData(sobj, layer = layer) # dgCMatrix (희소 행렬)
  
  # --- 2. 유전자 필터링 ---
  message(sprintf("1/8: 유전자 필터링 (min %d cells)...", min_count))
  keep_genes <- rowSums(counts > 0) >= min_count
  counts_filtered <- counts[keep_genes, ]
  message(sprintf("... %d / %d 유전자 통과", sum(keep_genes), nrow(counts)))
  
  # --- 3. NA 값 확인 및 제거 ---
  vars_to_check <- c(all_vars, offset)
  vars_to_check <- vars_to_check[vars_to_check %in% colnames(meta)]
  
  message("2/8: 모델 변수에서 NA 값 확인 중...")
  message(paste("... 확인 대상:", paste(vars_to_check, collapse = ", ")))
  
  # 각 변수별 NA 개수 확인
  if (remove_na_cells) {
    na_counts <- vapply(vars_to_check, function(v) sum(is.na(meta[[v]])), integer(1))
    if (any(na_counts > 0)) {
      message("... 변수별 NA 개수:")
      for (v in names(na_counts[na_counts > 0])) {
        message(sprintf("    %s: %d 개", v, na_counts[v]))
      }
    }
    
    keep_cells_idx <- complete.cases(meta[, vars_to_check])
    n_removed <- sum(!keep_cells_idx)
    
    if (n_removed > 0) {
      message(sprintf("... NA 값으로 인해 %d 개의 세포를 제거합니다.", n_removed))
    } else {
      message("... NA 값이 없습니다.")
    }
  } else {
    keep_cells_idx <- rep(TRUE, nrow(meta))
    n_removed <- 0
    if (any(!complete.cases(meta[, vars_to_check]))) {
      warning("일부 변수에 NA가 있지만 제거하지 않습니다. 분석 결과에 영향을 줄 수 있습니다.")
    }
  }
  
  # NA가 없는 '깨끗한' 데이터 생성
  meta_clean <- meta[keep_cells_idx, ]
  counts_clean <- counts_filtered[, keep_cells_idx]
  
  message(sprintf("... 최종 분석 대상 세포: %d 개", nrow(meta_clean)))
  
  separation_guard <- list(
    min_cells = separation_min_cells,
    min_cells_per_group = separation_min_cells_per_group,
    group_var = separation_group_var,
    removed_cells = 0L,
    removal_log = list(),
    balance_removal_log = list(),
    simplified_interactions = character(0),
    gene_balance = NULL
  )
  
  # --- 4. 디자인 행렬 생성 ---
  message("3/8: 디자인 행렬 생성 중...")
  
  # 범주형 변수 변환
  factor_vars <- all_vars[all_vars %in% colnames(meta_clean)]
  prunable_vars <- factor_vars[vapply(
    factor_vars,
    function(v) {
      col <- meta_clean[[v]]
      is.factor(col) || is.character(col) || length(unique(col)) <= 50
    },
    logical(1)
  )]
  meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], as.factor)

  prune_result <- ._nebula_prune_sparse_levels(
    meta = meta_clean,
    factor_vars = prunable_vars,
    min_cells_per_level = separation_min_cells
  )
  if (prune_result$removed_cells > 0) {
    meta_clean <- prune_result$meta
    counts_clean <- counts_clean[, prune_result$keep_idx, drop = FALSE]
    separation_guard$removed_cells <- separation_guard$removed_cells + prune_result$removed_cells
    separation_guard$removal_log <- prune_result$removal_log
    message(sprintf("... 희소 레벨 제거: %d 개의 세포 제외 (threshold=%d)", 
                    prune_result$removed_cells, separation_min_cells))
    for (var in names(prune_result$removal_log)) {
      message(sprintf("    %s → 제거된 레벨: %s",
                      var, paste(prune_result$removal_log[[var]], collapse = ", ")))
    }
  }
  if (!nrow(meta_clean)) {
    stop("희소 레벨 제거 후 남은 세포가 없습니다. threshold를 낮추거나 데이터를 확인하세요.")
  }
  
  balance_factor_candidates <- intersect(
    prunable_vars,
    unique(fixed_effects_from_nested)
  )
  if (!is.null(separation_group_var) &&
      separation_group_var %in% colnames(meta_clean) &&
      length(balance_factor_candidates) > 0) {
    balance_result <- ._nebula_prune_unbalanced_levels(
      meta = meta_clean,
      factor_vars = balance_factor_candidates,
      group_var = separation_group_var,
      min_cells_per_group = separation_min_cells_per_group
    )
    if (balance_result$removed_cells > 0) {
      meta_clean <- balance_result$meta
      counts_clean <- counts_clean[, balance_result$keep_idx, drop = FALSE]
      separation_guard$removed_cells <- separation_guard$removed_cells + balance_result$removed_cells
      separation_guard$balance_removal_log <- balance_result$removal_log
      message(sprintf("... 그룹 균형 기준으로 %d 개의 세포 제거 (threshold=%d per level)", 
                      balance_result$removed_cells, separation_min_cells_per_group))
      for (var in names(balance_result$removal_log)) {
        message(sprintf("    %s → 제거된 레벨 (group coverage 부족): %s",
                        var, paste(balance_result$removal_log[[var]], collapse = ", ")))
      }
    }
    if (!nrow(meta_clean)) {
      stop("그룹 균형 필터 적용 후 남은 세포가 없습니다. threshold를 낮추거나 데이터를 확인하세요.")
    }
  }
  meta_clean[factor_vars] <- lapply(meta_clean[factor_vars], droplevels)
  message(sprintf("... 희소/불균형 레벨 정리 후 세포: %d 개", nrow(meta_clean)))
  
  if (!is.null(separation_group_var) &&
      separation_group_var %in% colnames(meta_clean) &&
      nrow(counts_clean) > 0) {
    group_vec <- droplevels(meta_clean[[separation_group_var]])
    group_levels <- levels(group_vec)
    if (length(group_levels) >= 2) {
      message(sprintf("... %s 기준 유전자 균형 필터 적용 (레벨당 최소 %d 셀)",
                      separation_group_var, separation_min_cells_per_group))
      group_counts <- sapply(group_levels, function(gl) {
        mask <- group_vec == gl
        if (!any(mask)) {
          return(rep(0, nrow(counts_clean)))
        }
        Matrix::rowSums(counts_clean[, mask, drop = FALSE] > 0)
      })
      keep_gene_mask <- apply(group_counts >= separation_min_cells_per_group, 1, all)
      removed_genes <- sum(!keep_gene_mask)
      if (removed_genes > 0) {
        counts_clean <- counts_clean[keep_gene_mask, , drop = FALSE]
        separation_guard$gene_balance <- list(
          group_var = separation_group_var,
          threshold = separation_min_cells_per_group,
          removed_genes = removed_genes
        )
        message(sprintf("... 그룹 균형으로 %d 개 유전자 제거 (잔여 %d)",
                        removed_genes, nrow(counts_clean)))
      } else {
        message("... 모든 유전자가 그룹 균형 기준을 통과했습니다.")
      }
    } else {
      message(sprintf("... %s 컬럼에 단일 레벨만 있어 유전자 균형 필터를 건너뜁니다.",
                      separation_group_var))
    }
  }
  if (nrow(counts_clean) == 0) {
    stop("그룹 균형 필터 적용 후 남은 유전자가 없습니다. threshold를 낮추거나 데이터를 확인하세요.")
  }
  
  # 각 factor 변수의 레벨 수 확인
  for (v in factor_vars) {
    if (!v %in% colnames(meta_clean)) {
      next
    }
    n_levels <- length(levels(meta_clean[[v]]))
    levels_str <- paste(levels(meta_clean[[v]])[1:min(5, n_levels)], collapse=", ")
    if(n_levels > 5) {
      levels_str <- paste0(levels_str, "...")
    }
    message(sprintf("... %s: %d 레벨 (%s)", v, n_levels, levels_str))
  }
  
  # offset 처리
  if (offset %in% colnames(meta_clean)) {
    meta_clean[[offset]] <- as.numeric(meta_clean[[offset]])
    if (any(is.na(meta_clean[[offset]]))) {
        stop(sprintf("'%s' 컬럼에 NA가 아닌 숫자형 값만 있어야 합니다.", offset))
    }
    if (any(meta_clean[[offset]] <= 0)) {
      warning(sprintf("'%s' 컬럼에 0 이하 값이 있습니다. offset은 양수여야 합니다.", offset))
    }
  } else {
    stop(sprintf("'%s' 컬럼이 메타데이터에 없습니다.", offset))
  }
  
  # Formula에서 random effects 제거한 후 design matrix 생성
  interaction_terms <- fixed_effects_all[grepl(":", fixed_effects_all)]
  main_effects_only <- setdiff(fixed_effects_all, interaction_terms)
  
  build_design <- function(effect_terms) {
    if (!length(effect_terms)) {
      stop("Formula에서 사용할 고정 효과가 없습니다. 교호작용을 모두 제거하지 않았는지 확인하세요.")
    }
    formula_obj <- as.formula(paste("~", paste(effect_terms, collapse = " + ")))
    message(sprintf("... Design formula: %s", paste(deparse(formula_obj), collapse = " ")))
    matrix <- model.matrix(formula_obj, data = meta_clean)
    zero_var_cols <- which(apply(matrix, 2, function(x) length(unique(x)) <= 1))
    intercept_idx <- which(colnames(matrix) == "(Intercept)")
    zero_var_cols <- setdiff(zero_var_cols, intercept_idx)
    if (length(zero_var_cols) > 0) {
      message("... zero-variance predictors removed: ",
              paste(colnames(matrix)[zero_var_cols], collapse = ", "))
      matrix <- matrix[, -zero_var_cols, drop = FALSE]
    }
    list(matrix = matrix, formula = formula_obj)
  }
  
  build_result <- build_design(fixed_effects_all)
  design_matrix <- build_result$matrix
  design_formula <- build_result$formula
  design_qr <- qr(design_matrix)
  design_rank <- design_qr$rank
  n_cols <- ncol(design_matrix)
  
  message("... 변수 간 완전 분리(complete separation) 확인 중...")
  separation_info <- ._nebula_detect_complete_separation(meta_clean, main_effects_only)
  ._nebula_log_separation(separation_info)
  dropped_interactions <- character(0)
  if (simplify_interactions && separation_info$has_issue && length(interaction_terms) > 0) {
    problem_vars <- separation_info$problem_vars
    interactions_to_drop <- interaction_terms[vapply(
      interaction_terms,
      function(term) {
        term_vars <- strsplit(term, ":")[[1]]
        any(term_vars %in% problem_vars)
      },
      logical(1)
    )]
    if (length(interactions_to_drop)) {
      message("... 완전 분리 해소를 위해 교호작용 항 제거: ",
              paste(interactions_to_drop, collapse = ", "))
      fixed_effects_all <- setdiff(fixed_effects_all, interactions_to_drop)
      dropped_interactions <- interactions_to_drop
      interaction_terms <- setdiff(interaction_terms, interactions_to_drop)
      main_effects_only <- setdiff(fixed_effects_all, interaction_terms)
      build_result <- build_design(fixed_effects_all)
      design_matrix <- build_result$matrix
      design_formula <- build_result$formula
      design_qr <- qr(design_matrix)
      design_rank <- design_qr$rank
      n_cols <- ncol(design_matrix)
      if (length(main_effects_only) >= 2) {
        message("... 교호작용 제거 후 완전 분리 재확인...")
        separation_info <- ._nebula_detect_complete_separation(meta_clean, main_effects_only)
        ._nebula_log_separation(separation_info)
      }
    }
  }
  separation_guard$simplified_interactions <- dropped_interactions
  separation_guard$diagnostics <- separation_info
  
  if (design_rank < n_cols) {
    warning(sprintf("설계 행렬이 특이(singular)합니다. rank=%d < columns=%d.", 
                   design_rank, n_cols))
    dependent_idx <- design_qr$pivot[seq.int(design_rank + 1, n_cols)]
    drop_cols <- colnames(design_matrix)[dependent_idx]
    message("... Removing linearly dependent predictors: ", paste(drop_cols, collapse = ", "))
    keep_idx <- design_qr$pivot[seq_len(design_rank)]
    design_matrix <- design_matrix[, keep_idx, drop = FALSE]
    design_qr <- qr(design_matrix)
    design_rank <- design_qr$rank
    n_cols <- ncol(design_matrix)
    message("오류가 발생하면:")
    message("  1. 교호작용 항을 제거하거나 단순화하세요.")
    message("  2. 완전 분리된 조합을 제거하거나 데이터를 필터링하세요.")
    message("  3. min_count를 높여서 더 적은 유전자로 분석하세요.")
  }
  
  message(sprintf("... 설계 행렬: %d 행 x %d 열 (rank=%d)", 
                 nrow(design_matrix), ncol(design_matrix), design_rank))
  
  # id 및 offset 벡터 추출
  id_vector <- meta_clean[[patient_col]]
  offset_vector <- meta_clean[[offset]]
  
  # --- 5. group_cell()로 데이터 정렬 ---
  message(sprintf("4/8: NEBULA를 위해 id (%s) 기준으로 데이터 정렬 중...", patient_col))
  data_grouped <- nebula::group_cell(
    count = counts_clean,
    id = id_vector,
    pred = design_matrix,
    offset = offset_vector
  )
  
  # --- 6. NEBULA 실행 ---
  message("5/8: NEBULA 실행 중 (NBLMM)...")
  re_nebula <- ._nebula_fit_with_retry(data_grouped)
  
  # --- 7. 결과에 formula 정보 추가 ---
  re_nebula$formula <- deparse(formula_obj)
  re_nebula$design_formula <- deparse(design_formula)
  re_nebula$patient_col <- patient_col
  re_nebula$fixed_effects <- fixed_effects_all
  re_nebula$separation_guard <- separation_guard
  
  # --- 8. 완료 ---
  message("6/8: 분석 완료.")
  
  return(re_nebula)
}

