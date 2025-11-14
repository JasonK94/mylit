# ============================================================================
# LDS 단계 8: SVA 상관관계 Heatmap 생성
# ============================================================================

#' @title LDS 단계 8: SVA 상관관계 Heatmap 생성
#' @description SVA와 메타데이터 간 상관관계를 3가지 방식으로 시각화
#' @param svs_final 사용된 SV 매트릭스
#' @param meta.data 메타데이터 (SV 포함 가능)
#' @param output_dir 출력 디렉터리
#' @param prefix 파일명 접두사
#' @param method 상관관계 계산 방법 ("spearman" 또는 "pearson")
#' @param top_n 상위 n개 메타데이터 변수만 표시 (p-value 기준). 기본값 10
#' @param p_threshold p-value 임계값. 기본값 0.05
#' @param save_intermediate 중간 결과 저장 여부
#' @return list(heatmap_files, correlation_matrices)
#' @export
lds_08_create_heatmaps <- function(svs_final,
                                    meta.data,
                                    output_dir = tempdir(),
                                    prefix = "lds",
                                    method = c("spearman", "pearson"),
                                    top_n = 10,
                                    p_threshold = 0.05,
                                    save_intermediate = FALSE) {
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("'pheatmap' 패키지가 필요합니다: install.packages('pheatmap')")
  }
  
  method <- match.arg(method)
  
  if (is.null(svs_final) || ncol(svs_final) == 0) {
    message("8/8: SV가 없어 heatmap 생성을 스킵합니다.")
    return(NULL)
  }
  
  message("8/8: SVA 상관관계 Heatmap 생성 중...")
  
  # 메타데이터에서 숫자형/팩터형/문자형 변수만 선택
  meta_vars <- colnames(meta.data)
  # SV 변수 제외
  meta_vars <- meta_vars[!grepl("^SV\\d+$", meta_vars)]
  # 숫자형, 팩터형, 또는 문자형 변수만 (문자형은 나중에 factor로 변환)
  numeric_vars <- sapply(meta_vars, function(v) {
    var_data <- meta.data[[v]]
    is.numeric(var_data) || is.factor(var_data) || is.character(var_data)
  })
  meta_vars_numeric <- meta_vars[numeric_vars]
  
  if (length(meta_vars_numeric) == 0) {
    message("... 숫자형/팩터형/문자형 메타데이터 변수가 없어 heatmap 생성을 스킵합니다.")
    return(NULL)
  }
  
  # 메타데이터 준비 (문자형을 factor로 변환 후 숫자로 변환)
  meta_numeric <- meta.data[, meta_vars_numeric, drop = FALSE]
  for (i in 1:ncol(meta_numeric)) {
    if (is.character(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.factor(meta_numeric[[i]])
    }
    if (is.factor(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.numeric(meta_numeric[[i]])
    }
  }
  
  sv_names <- colnames(svs_final)
  n_sv <- length(sv_names)
  
  # ============================================================================
  # 1. (metadata + SV들)^2 matrix (전체 상관관계 행렬)
  # ============================================================================
  message("... 1/3: 전체 상관관계 행렬 생성 중...")
  
  # SV와 메타데이터를 결합
  combined_data <- cbind(meta_numeric, svs_final)
  combined_names <- c(meta_vars_numeric, sv_names)
  
  # 전체 상관관계 행렬 계산
  cor_matrix_full <- matrix(NA, nrow = length(combined_names), ncol = length(combined_names))
  rownames(cor_matrix_full) <- combined_names
  colnames(cor_matrix_full) <- combined_names
  
  for (i in 1:length(combined_names)) {
    for (j in 1:length(combined_names)) {
      vec1 <- combined_data[, i]
      vec2 <- combined_data[, j]
      
      valid_idx <- !is.na(vec1) & !is.na(vec2)
      if (sum(valid_idx) > 1) {
        cor_matrix_full[i, j] <- cor(vec1[valid_idx], vec2[valid_idx], method = method)
      }
    }
  }
  
  # 대각선은 1로 설정
  diag(cor_matrix_full) <- 1
  
  # Heatmap 1 저장
  output_file_1 <- file.path(output_dir, paste0(prefix, "_heatmap_full.png"))
  png(output_file_1, width = 2000, height = 2000, res = 150)
  tryCatch({
    pheatmap::pheatmap(
      cor_matrix_full,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 7,
      na_col = "gray90"
    )
  }, error = function(e) {
    # 클러스터링 실패 시 NA를 0으로 대체
    cor_matrix_full_clust <- cor_matrix_full
    cor_matrix_full_clust[is.na(cor_matrix_full_clust)] <- 0
    diag(cor_matrix_full_clust) <- 1
    pheatmap::pheatmap(
      cor_matrix_full_clust,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 7
    )
  })
  dev.off()
  message("  ✓ 저장: ", output_file_1)
  
  # ============================================================================
  # 2. metadata x SV (모든 메타데이터 변수)
  # ============================================================================
  message("... 2/3: 메타데이터 x SV 상관관계 행렬 생성 중...")
  
  # 메타데이터(행) x SV(열) 상관관계 계산
  cor_matrix_meta_sv <- matrix(NA, nrow = length(meta_vars_numeric), ncol = n_sv)
  rownames(cor_matrix_meta_sv) <- meta_vars_numeric
  colnames(cor_matrix_meta_sv) <- sv_names
  
  for (i in seq_along(meta_vars_numeric)) {
    for (j in seq_along(sv_names)) {
      meta_vec <- meta_numeric[[i]]
      sv_vec <- svs_final[, j]
      
      valid_idx <- !is.na(meta_vec) & !is.na(sv_vec)
      if (sum(valid_idx) > 1) {
        cor_matrix_meta_sv[i, j] <- cor(meta_vec[valid_idx], sv_vec[valid_idx], method = method)
      }
    }
  }
  
  # Heatmap 2 저장
  output_file_2 <- file.path(output_dir, paste0(prefix, "_heatmap_metadata_x_sv.png"))
  png(output_file_2, width = 1200, height = max(800, length(meta_vars_numeric) * 15), res = 150)
  tryCatch({
    pheatmap::pheatmap(
      cor_matrix_meta_sv,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = FALSE,  # SV는 순서 유지
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 8,
      fontsize_row = 8,
      fontsize_col = 10,
      na_col = "gray90"
    )
  }, error = function(e) {
    cor_matrix_meta_sv_clust <- cor_matrix_meta_sv
    cor_matrix_meta_sv_clust[is.na(cor_matrix_meta_sv_clust)] <- 0
    pheatmap::pheatmap(
      cor_matrix_meta_sv_clust,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 8,
      fontsize_row = 8,
      fontsize_col = 10
    )
  })
  dev.off()
  message("  ✓ 저장: ", output_file_2)
  
  # ============================================================================
  # 3. top n metadata x SV (p-value < 0.05 기준)
  # ============================================================================
  message("... 3/3: 상위 n개 메타데이터 x SV 상관관계 행렬 생성 중 (p<", p_threshold, ")...")
  
  # p-value 계산
  p_value_matrix <- matrix(NA, nrow = length(meta_vars_numeric), ncol = n_sv)
  rownames(p_value_matrix) <- meta_vars_numeric
  colnames(p_value_matrix) <- sv_names
  
  for (i in seq_along(meta_vars_numeric)) {
    for (j in seq_along(sv_names)) {
      meta_vec <- meta_numeric[[i]]
      sv_vec <- svs_final[, j]
      
      valid_idx <- !is.na(meta_vec) & !is.na(sv_vec)
      if (sum(valid_idx) > 2) {
        cor_test <- tryCatch({
          cor.test(meta_vec[valid_idx], sv_vec[valid_idx], method = method)
        }, error = function(e) NULL)
        
        if (!is.null(cor_test)) {
          p_value_matrix[i, j] <- cor_test$p.value
        }
      }
    }
  }
  
  # p-value < threshold인 변수 찾기
  significant_vars <- rownames(p_value_matrix)[rowSums(p_value_matrix < p_threshold, na.rm = TRUE) > 0]
  
  if (length(significant_vars) == 0) {
    message("  ⚠ p < ", p_threshold, "인 메타데이터 변수가 없습니다.")
    # 절댓값 상관관계가 높은 변수 선택
    abs_cor <- abs(cor_matrix_meta_sv)
    abs_cor[is.na(abs_cor)] <- 0
    max_cor_per_var <- apply(abs_cor, 1, max, na.rm = TRUE)
    significant_vars <- names(sort(max_cor_per_var, decreasing = TRUE))[1:min(top_n, length(max_cor_per_var))]
    message("  → 절댓값 상관관계가 높은 상위 ", length(significant_vars), "개 변수를 선택합니다.")
  } else {
    # p-value가 유의한 변수 중 상위 n개 선택 (절댓값 상관관계 기준)
    abs_cor_sig <- abs(cor_matrix_meta_sv[significant_vars, , drop = FALSE])
    abs_cor_sig[is.na(abs_cor_sig)] <- 0
    max_cor_per_var <- apply(abs_cor_sig, 1, max, na.rm = TRUE)
    significant_vars <- names(sort(max_cor_per_var, decreasing = TRUE))[1:min(top_n, length(significant_vars))]
  }
  
  if (length(significant_vars) > 0) {
    cor_matrix_top <- cor_matrix_meta_sv[significant_vars, , drop = FALSE]
    p_value_matrix_top <- p_value_matrix[significant_vars, , drop = FALSE]
    
    # Heatmap 3 저장 (p-value 표시)
    output_file_3 <- file.path(output_dir, paste0(prefix, "_heatmap_top", top_n, "_metadata_x_sv.png"))
    png(output_file_3, width = 1000, height = max(600, length(significant_vars) * 20), res = 150)
    
    # p-value를 별표로 표시
    display_matrix <- matrix("", nrow = nrow(cor_matrix_top), ncol = ncol(cor_matrix_top))
    rownames(display_matrix) <- rownames(cor_matrix_top)
    colnames(display_matrix) <- colnames(cor_matrix_top)
    
    for (i in 1:nrow(p_value_matrix_top)) {
      for (j in 1:ncol(p_value_matrix_top)) {
        if (!is.na(p_value_matrix_top[i, j])) {
          if (p_value_matrix_top[i, j] < 0.001) {
            display_matrix[i, j] <- "***"
          } else if (p_value_matrix_top[i, j] < 0.01) {
            display_matrix[i, j] <- "**"
          } else if (p_value_matrix_top[i, j] < 0.05) {
            display_matrix[i, j] <- "*"
          }
        }
      }
    }
    
    tryCatch({
      pheatmap::pheatmap(
        cor_matrix_top,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-1, 1, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        display_numbers = display_matrix,
        number_format = "%s",
        fontsize = 9,
        fontsize_row = 9,
        fontsize_col = 10,
        fontsize_number = 8,
        na_col = "gray90"
      )
    }, error = function(e) {
      cor_matrix_top_clust <- cor_matrix_top
      cor_matrix_top_clust[is.na(cor_matrix_top_clust)] <- 0
      pheatmap::pheatmap(
        cor_matrix_top_clust,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = seq(-1, 1, length.out = 101),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        display_numbers = display_matrix,
        number_format = "%s",
        fontsize = 9,
        fontsize_row = 9,
        fontsize_col = 10,
        fontsize_number = 8
      )
    })
    dev.off()
    message("  ✓ 저장: ", output_file_3, " (", length(significant_vars), "개 변수)")
  } else {
    message("  ⚠ 표시할 변수가 없습니다.")
    output_file_3 <- NULL
  }
  
  # 결과 반환
  result <- list(
    heatmap_files = c(
      full = output_file_1,
      metadata_x_sv = output_file_2,
      top_metadata_x_sv = output_file_3
    ),
    correlation_matrices = list(
      full = cor_matrix_full,
      metadata_x_sv = cor_matrix_meta_sv,
      top_metadata_x_sv = if (length(significant_vars) > 0) cor_matrix_top else NULL
    ),
    p_value_matrix = p_value_matrix,
    significant_vars = significant_vars
  )
  
  if (save_intermediate) {
    output_path <- file.path(output_dir, paste0(prefix, "_08_heatmaps.qs"))
    qs::qsave(result, output_path)
    message("... 중간 결과 저장: ", output_path)
  }
  
  message("8/8: Heatmap 생성 완료.")
  
  return(result)
}

