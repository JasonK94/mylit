#' 클러스터 비율 비교 박스플롯 생성 (Pseudobulk 방식)
#'
#' @param sobj Seurat 객체
#' @param patient_col 메타데이터 내 환자 또는 샘플 ID 컬럼명
#' @param group.by 비교할 클러스터 컬럼명
#' @param split.by 비교할 그룹 컬럼명
#' @param test.use 사용할 통계 검정 ("wilcox.test", "t.test", "kruskal.test", "anova")
#' @param p.adj p-value 조정 방법 ("BH", "bonferroni", "none" 등)
#' @param split.na.use 'split.by' 컬럼의 NA 값 처리 방법 (FALSE: 제외, TRUE: 별도 그룹으로 취급)
#' @param colors 그룹별로 사용할 색상 팔레트 (선택 사항)
#'
#' @return ggplot 객체와 통계 결과가 포함된 리스트
#'
#' @export
PlotClusterProportions <- function(sobj, 
                                   patient_col, 
                                   group.by, 
                                   split.by, 
                                   test.use = "wilcox.test",
                                   p.adj = "BH",
                                   split.na.use = FALSE,
                                   colors = NULL) {
  
  # 1. 메타데이터 추출
  meta_data <- sobj@meta.data
  
  # 2. 필수 컬럼 확인
  required_cols <- c(patient_col, group.by, split.by)
  if (!all(required_cols %in% colnames(meta_data))) {
    stop("제공된 컬럼명 중 일부가 메타데이터에 없습니다: ", 
         paste(required_cols[!required_cols %in% colnames(meta_data)], collapse = ", "))
  }
  
  # 3. 환자(샘플)별 그룹 정보 추출 (dplyr:: 명시적 호출)
  patient_group_info <- meta_data %>%
    dplyr::select(all_of(c(patient_col, split.by))) %>%
    dplyr::distinct()
  
  # 4. Pseudobulk 비율 계산 (dplyr:: 명시적 호출)
  prop_data <- meta_data %>%
    dplyr::group_by(!!sym(patient_col), !!sym(group.by)) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(!!sym(patient_col)) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(patient_group_info, by = patient_col)
  
  # 5. NA 처리 로직 (dplyr:: 명시적 호출)
  any_na <- any(is.na(prop_data[[split.by]]))
  if (any_na) {
    if (split.na.use == FALSE) {
      message("NA values found in '", split.by, "'. Ignoring samples with NA.")
      prop_data <- prop_data %>% dplyr::filter(!is.na(!!sym(split.by)))
    } else { 
      message("NA values found in '", split.by, "'. Treating as a separate group: 'NA_as_Group'.")
      prop_data <- prop_data %>%
        dplyr::mutate(!!sym(split.by) := as.character(!!sym(split.by))) %>%
        dplyr::mutate(!!sym(split.by) := dplyr::if_else(is.na(!!sym(split.by)), "NA_as_Group", !!sym(split.by))) %>%
        dplyr::mutate(!!sym(split.by) := factor(!!sym(split.by)))
    }
  }
  
  # 5.5. 팩터 레벨 정리
  prop_data <- droplevels(prop_data)
  
  # 5.6. split.by 컬럼을 factor로 변환 (V8 수정)
  prop_data[[split.by]] <- as.factor(prop_data[[split.by]])
  
  # *** 5.7. group.by 컬럼을 factor로 변환 (V10 수정) ***
  # 이 부분이 .subset2 오류의 최종 원인일 가능성이 높습니다.
  prop_data[[group.by]] <- as.factor(prop_data[[group.by]])
  
  # 6. 통계 검정을 위한 데이터 필터링
  n_groups <- length(unique(prop_data[[split.by]]))
  stat_test <- data.frame() # 통계 결과 초기화
  
  if (n_groups < 2) {
    warning("NA 제거 후 비교할 그룹이 2개 미만입니다. 통계 검정을 생략합니다.")
  } else {
    # 7. 비교 가능한 클러스터 탐색 (dplyr:: 명시적 호출)
    groups_per_cluster <- prop_data %>%
      dplyr::group_by(!!sym(group.by)) %>%
      dplyr::summarise(n_groups_in_cluster = dplyr::n_distinct(!!sym(split.by)))
    
    clusters_to_test <- groups_per_cluster %>%
      dplyr::filter(n_groups_in_cluster == n_groups) %>%
      dplyr::pull(!!sym(group.by))
    
    # 8. 통계 검정 수행 (dplyr:: 명시적 호출)
    if (length(clusters_to_test) > 0) {
      
      prop_data_for_stat_test <- prop_data %>%
        dplyr::filter(!!sym(group.by) %in% clusters_to_test)
      
      if (n_groups > 2 && test.use == "t.test") {
        message("그룹이 3개 이상이므로 'anova'를 사용합니다.")
        test.use <- "anova"
      }
      if (n_groups > 2 && test.use == "wilcox.test") {
        message("그룹이 3개 이상이므로 'kruskal.test'를 사용합니다.")
        test.use <- "kruskal.test"
      }
      
      stat_test <- prop_data_for_stat_test %>%
        dplyr::group_by(!!sym(group.by)) %>%
        dplyr::do(
          if(test.use == "wilcox.test") {
            rstatix::wilcox_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else if (test.use == "t.test") {
            rstatix::t_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else if (test.use == "kruskal.test") {
            rstatix::kruskal_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else if (test.use == "anova") {
            rstatix::anova_test(data = ., as.formula(paste("proportion ~", split.by)))
          } else {
            stop("지원하지 않는 test.use입니다.")
          }
        ) %>%
        dplyr::ungroup()
      
      # 9. P-value 조정
      if (nrow(stat_test) > 0 && p.adj != "none") {
        stat_test <- stat_test %>%
          rstatix::adjust_pvalue(method = p.adj) %>%
          rstatix::add_significance(p.col = "p.adj")
        stat_test$p_label_col <- "p.adj.signif"
      } else if (nrow(stat_test) > 0) {
        stat_test <- stat_test %>%
          rstatix::add_significance(p.col = "p")
        stat_test$p_label_col <- "p.signif"
      }
      
      # 10. P-value 브라켓 위치 계산
      if (nrow(stat_test) > 0 && n_groups == 2) {
        stat_test <- stat_test %>%
          rstatix::add_xy_position(
            data = prop_data, 
            formula = as.formula(paste("proportion ~", split.by)), 
            x = group.by, 
            dodge = 0.8
          )
      }
      
    } else {
      warning("NA 제거 후, 모든 그룹 간 비교가 가능한 클러스터가 없습니다. 통계 검정을 생략합니다.")
    }
  }
  
  # 11. 박스플롯 생성
  p <- ggplot2::ggplot(prop_data, ggplot2::aes(x = !!sym(group.by), y = proportion, fill = !!sym(split.by))) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                        alpha = 0.7, size = 1.5) +
    ggplot2::labs(
      title = "Cluster Proportion by Group (Pseudobulk)",
      x = group.by,
      y = "Proportion per Sample",
      fill = split.by
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent)
  
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  
  # 12. P-value 브라켓 추가
  if (nrow(stat_test) > 0) {
    if (n_groups > 2) {
      p <- p + ggpubr::stat_pvalue_manual(
        data = stat_test,
        aes(x = !!sym(group.by)), 
        label = "{p_label_col}",
        y.position = max(prop_data$proportion, na.rm = TRUE) * 1.05,
        hide.ns = TRUE,
        inherit.aes = FALSE 
      )
    } else {
      p <- p + ggpubr::stat_pvalue_manual(
        stat_test, 
        label = "{p_label_col}",
        tip.length = 0.01,
        hide.ns = TRUE
      )
    }
  }
  
  # 13. 반환
  return(list(plot = p, statistics = stat_test))
}

#' FindAllMarkers 결과에서 Top N 유전자를 LLM용으로 출력
#'
#' @param markers_df FindAllMarkers()의 결과인 데이터프레임
#' @param top_n 각 클러스터별로 추출할 상위 유전자 개수
#' @param cluster_col 'cluster' 정보가 있는 컬럼명 (기본값: "cluster")
#' @param gene_col 'gene' 정보가 있는 컬럼명 (기본값: "gene")
#' @param p_val_col p-value 컬럼명 (기본값: "p_val")
#' @param p_val_adj_col p-value (adjusted) 컬럼명 (기본값: "p_val_adj")
#' @param avg_log2FC_col log-fold change 컬럼명 (기본값: "avg_log2FC")
#'
#' @return 콘솔에 클러스터별 Top N 유전자 목록을 출력합니다.
#'
#' @export
PrintTopMarkersForLLM <- function(markers_df,
                                  top_n = 50,
                                  cluster_col = "cluster",
                                  gene_col = "gene",
                                  p_val_col = "p_val",
                                  p_val_adj_col = "p_val_adj",
                                  avg_log2FC_col = "avg_log2FC") {
  
  # 1. 필수 컬럼 존재 여부 확인
  required_cols <- c(cluster_col, gene_col, avg_log2FC_col)
  
  if (!all(required_cols %in% colnames(markers_df))) {
    missing_cols <- required_cols[!required_cols %in% colnames(markers_df)]
    stop("필수 컬럼이 markers_df에 없습니다: ", paste(missing_cols, collapse = ", "),
         "\n'cluster_col', 'gene_col', 'avg_log2FC_col' 인자를 확인하세요.")
  }
  
  # 2. 정렬을 위한 데이터 전처리 (p-value 0 처리)
  markers_processed <- markers_df
  
  use_pval <- p_val_col %in% colnames(markers_processed)
  use_padj <- p_val_adj_col %in% colnames(markers_processed)
  
  # p_val_mod 컬럼을 안전하게 추가
  if (use_pval) {
    markers_processed$p_val_mod <- ifelse(markers_processed[[p_val_col]] == 0, 
                                          1e-300, # 0 대신 매우 작은 값
                                          markers_processed[[p_val_col]])
  }
  
  # 3. 정렬 기준에 따라 Top N 마커 필터링
  markers_grouped <- markers_processed %>%
    group_by(!!sym(cluster_col))
  
  # 정렬 로직
  if (use_padj && use_pval) {
    message(paste("Sorting by", p_val_adj_col, ",", p_val_col, ", and", avg_log2FC_col))
    sorted_markers <- markers_grouped %>%
      arrange(!!sym(p_val_adj_col), p_val_mod, desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
    
  } else if (use_pval) {
    message(paste("Sorting by", p_val_col, "and", avg_log2FC_col))
    sorted_markers <- markers_grouped %>%
      arrange(p_val_mod, desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
    
  } else {
    message(paste("Sorting by", avg_log2FC_col, "only."))
    sorted_markers <- markers_grouped %>%
      arrange(desc(!!sym(avg_log2FC_col)), .by_group = TRUE)
  }
  
  final_markers <- sorted_markers %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # 4. LLM에 복사/붙여넣기 편한 형태로 콘솔에 출력
  
  # 클러스터 목록 (숫자형/문자형 모두 정렬되도록)
  clusters <- unique(final_markers[[cluster_col]])
  tryCatch({
    clusters_sorted <- sort(as.numeric(as.character(clusters)))
  }, warning = function(w) {
    clusters_sorted <<- sort(as.character(clusters))
  })
  
  
  cat("--- Top", top_n, "Markers per Cluster ---\n\n")
  
  for (cl in clusters_sorted) {
    cat("--- Cluster:", cl, "---\n")
    
    genes <- final_markers %>%
      filter(!!sym(cluster_col) == cl) %>%
      pull(!!sym(gene_col))
    
    # 쉼표로 구분된 한 줄의 문자열로 출력
    cat(paste(genes, collapse = ", "))
    cat("\n\n") # 클러스터 간 구분을 위한 공백
  }
  
  message("--- End of Marker List ---")
}

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
  
  stats_results <- map_dfr(all_cell_types, function(current_cell_type) {
    
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
  pb_permanova_results <- map_dfr(all_cell_types, function(current_cell_type) {
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
    rownames(meta_pb_subset) <- paste0("g",meta_pb_subset[[patient_var]])
    meta_pb_subset <- meta_pb_subset[rownames(pb_matrix_t), ]
    
    # (오류 방지) 해당 세포 유형에 샘플이 너무 적거나 그룹이 1개면 스킵
    if (nrow(meta_pb_subset) < 3 || length(unique(meta_pb_subset[[group_var]])) < 2) {
      # --- [확인용 코드 추가] ---
      message(paste("  [Skipping]", current_cell_type, 
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
      permutations = 999
      S  )
    
    # 9. 결과 추출
    data.frame(
      cell_type = current_cell_type,
      R2 = res$R2[1],        # 'set' 그룹이 설명하는 분산의 양
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