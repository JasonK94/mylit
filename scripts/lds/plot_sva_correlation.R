# SVA와 메타데이터 상관관계 플롯 생성
# 사용법: source("scripts/lds/plot_sva_correlation.R")

suppressPackageStartupMessages({
  library(qs)
  library(corrplot)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
})

# 1. 결과 로드
cat("결과 로드 중...\n")
result <- qs::qread("/data/user3/sobj/test_lds_final_result.qs")

# 2. 원본 데이터 로드 (메타데이터용)
cat("원본 데이터 로드 중...\n")
is5s <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
is5s_test <- is5s
is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]

# 3. SV 추출
sv_used <- result$svs_used
cat("사용된 SV:", ncol(sv_used), "개\n")
cat("SV 이름:", paste(colnames(sv_used), collapse=", "), "\n")

# 4. 메타데이터 준비
meta <- is5s_test@meta.data
cat("메타데이터 변수:", ncol(meta), "개\n")

# 5. SV와 메타데이터 결합
if (nrow(sv_used) == nrow(meta)) {
  # 숫자형/팩터형 변수만 선택
  numeric_vars <- sapply(colnames(meta), function(v) {
    is.numeric(meta[[v]]) || is.factor(meta[[v]])
  })
  meta_numeric <- meta[, numeric_vars, drop = FALSE]
  
  # 팩터를 숫자로 변환
  for (i in 1:ncol(meta_numeric)) {
    if (is.factor(meta_numeric[[i]])) {
      meta_numeric[[i]] <- as.numeric(meta_numeric[[i]])
    }
  }
  
  cat("숫자형/팩터형 변수:", ncol(meta_numeric), "개\n")
  
  # 6. 상관관계 계산
  cat("상관관계 계산 중...\n")
  sv_names <- colnames(sv_used)
  meta_var_names <- colnames(meta_numeric)
  
  cor_matrix <- matrix(NA, nrow = length(sv_names), ncol = length(meta_var_names))
  rownames(cor_matrix) <- sv_names
  colnames(cor_matrix) <- meta_var_names
  
  for (i in seq_along(sv_names)) {
    for (j in seq_along(meta_var_names)) {
      sv_vec <- sv_used[, i]
      meta_vec <- meta_numeric[[j]]
      
      # NA 제거 후 상관관계 계산
      valid_idx <- !is.na(sv_vec) & !is.na(meta_vec)
      if (sum(valid_idx) > 1) {
        cor_matrix[i, j] <- cor(sv_vec[valid_idx], meta_vec[valid_idx], method = "pearson")
      }
    }
  }
  
  cat("상관관계 행렬 생성 완료:", nrow(cor_matrix), "x", ncol(cor_matrix), "\n")
  
  # 7. 상관관계 플롯 생성
  
  # 7-1. corrplot 사용
  cat("\ncorrplot 생성 중...\n")
  png("/data/user3/sobj/sva_correlation_corrplot.png", width = 2000, height = 1200, res = 150)
  corrplot(cor_matrix, 
           method = "color",
           type = "full",
           order = "hclust",
           tl.cex = 0.6,
           tl.col = "black",
           col = colorRampPalette(c("blue", "white", "red"))(200),
           cl.lim = c(-1, 1),
           is.corr = TRUE)
  dev.off()
  cat("✓ 저장: /data/user3/sobj/sva_correlation_corrplot.png\n")
  
  # 7-2. pheatmap 사용 (NA 처리)
  cat("pheatmap 생성 중...\n")
  # NA를 0으로 대체 (클러스터링을 위해)
  cor_matrix_for_clust <- cor_matrix
  cor_matrix_for_clust[is.na(cor_matrix_for_clust)] <- 0
  
  png("/data/user3/sobj/sva_correlation_pheatmap.png", width = 2000, height = 1200, res = 150)
  tryCatch({
    pheatmap(cor_matrix,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-1, 1, length.out = 101),
             cluster_rows = FALSE,  # NA 때문에 클러스터링 비활성화
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize = 8,
             fontsize_row = 10,
             fontsize_col = 8,
             na_col = "gray90")  # NA를 회색으로 표시
  }, error = function(e) {
    # 클러스터링 실패 시 NA를 0으로 대체하여 재시도
    pheatmap(cor_matrix_for_clust,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-1, 1, length.out = 101),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize = 8,
             fontsize_row = 10,
             fontsize_col = 8)
  })
  dev.off()
  cat("✓ 저장: /data/user3/sobj/sva_correlation_pheatmap.png\n")
  
  # 7-3. 상위 상관관계만 선택하여 플롯 (변수가 너무 많을 경우)
  if (ncol(cor_matrix) > 50) {
    cat("\n변수가 많아 상위 상관관계만 선택합니다...\n")
    
    # 각 SV에 대해 가장 높은 상관관계를 가진 변수 선택
    top_vars <- c()
    for (i in 1:nrow(cor_matrix)) {
      abs_cor <- abs(cor_matrix[i, ])
      abs_cor[is.na(abs_cor)] <- 0
      top_idx <- which.max(abs_cor)
      if (abs_cor[top_idx] > 0.1) {  # 최소 상관관계 0.1
        top_vars <- c(top_vars, colnames(cor_matrix)[top_idx])
      }
    }
    
    # 중복 제거
    top_vars <- unique(top_vars)
    
    if (length(top_vars) > 0) {
      cor_matrix_top <- cor_matrix[, top_vars, drop = FALSE]
      
      png("/data/user3/sobj/sva_correlation_top.png", width = 1500, height = 800, res = 150)
      pheatmap(cor_matrix_top,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-1, 1, length.out = 101),
               cluster_rows = FALSE,
               cluster_cols = TRUE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize = 10,
               na_col = "gray90")
      dev.off()
      cat("✓ 저장: /data/user3/sobj/sva_correlation_top.png (", length(top_vars), "개 변수)\n", sep="")
    }
  }
  
  # 8. 결과 저장
  saveRDS(cor_matrix, "/data/user3/sobj/sva_metadata_correlation.rds")
  cat("\n✓ 상관관계 행렬 저장: /data/user3/sobj/sva_metadata_correlation.rds\n")
  
  # 9. 상위 상관관계 출력
  cat("\n=== 상위 상관관계 (절댓값 기준) ===\n")
  abs_cor <- abs(cor_matrix)
  abs_cor[is.na(abs_cor)] <- 0
  top_indices <- order(abs_cor, decreasing = TRUE)[1:min(20, length(abs_cor))]
  
  for (idx in top_indices) {
    row_idx <- ((idx - 1) %/% ncol(cor_matrix)) + 1
    col_idx <- ((idx - 1) %% ncol(cor_matrix)) + 1
    if (row_idx <= nrow(cor_matrix) && col_idx <= ncol(cor_matrix)) {
      sv_name <- rownames(cor_matrix)[row_idx]
      meta_name <- colnames(cor_matrix)[col_idx]
      cor_val <- cor_matrix[row_idx, col_idx]
      if (!is.na(cor_val) && abs(cor_val) > 0.1) {
        cat(sprintf("%s <-> %s: %.3f\n", sv_name, meta_name, cor_val))
      }
    }
  }
  
} else {
  cat("⚠ SV와 메타데이터의 행 수가 일치하지 않습니다.\n")
  cat("  SV:", nrow(sv_used), ", 메타데이터:", nrow(meta), "\n")
}

cat("\n=== 완료 ===\n")

