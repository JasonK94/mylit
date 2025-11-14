# LDS 함수 인터랙티브 테스트
# R 세션에서 직접 실행: source("scripts/lds/test_lds_interactive.R")

cat("=== LDS 함수 테스트 시작 ===\n\n")

# 1. 패키지 로드
cat("1. 패키지 로드 중...\n")
suppressPackageStartupMessages({
  library(qs)
  library(limma)
  library(edgeR)
  library(lme4)
  library(BiocParallel)
  library(sva)
  library(dplyr)
  library(Seurat)
})
cat("✓ 패키지 로드 완료\n\n")

# 2. myR 패키지 로드
cat("2. myR 패키지 로드 중...\n")
# 워크트리에서 직접 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/lds/myR", quiet = FALSE)
if (!exists("LDS")) {
  stop("LDS 함수를 찾을 수 없습니다.")
}
cat("✓ LDS 함수 확인 완료\n")
cat("  함수 파라미터:", paste(names(formals(LDS)), collapse = ", "), "\n\n")

# 3. 데이터 로드
cat("3. 데이터 로드 중...\n")
data_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
if (!file.exists(data_path)) {
  stop("데이터 파일을 찾을 수 없습니다: ", data_path)
}
is5s <- qs::qread(data_path)
cat("✓ 데이터 로드 완료\n")
cat("  - 셀 수:", ncol(is5s), "\n")
cat("  - 유전자 수:", nrow(is5s), "\n\n")

# 4. 데이터 전처리
cat("4. 데이터 전처리 중...\n")
is5s_test <- is5s
is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
n_before <- ncol(is5s_test)
is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]
n_after <- ncol(is5s_test)
cat("✓ 전처리 완료\n")
cat("  - NA 제거: ", n_before, " -> ", n_after, " 셀\n", sep = "")
cat("  - g3_clean 값:", paste(sort(unique(is5s_test$g3_clean)), collapse = ", "), "\n")
cat("  - hos_no 샘플 수:", length(unique(is5s_test@meta.data$hos_no)), "\n")
cat("  - GEM 배치 수:", length(unique(is5s_test@meta.data$GEM)), "\n\n")

# 5. LDS 실행
cat("5. LDS 함수 실행 중...\n")
cat("   포뮬러: ~ g3_clean + (1|hos_no) + (1|GEM)\n")
cat("   SV: 자동 결정 (sv_var_cutoff = 0.5)\n")
cat("   필터링 옵션: min.count=5, min.prop=0.05 (완화)\n\n")

result <- tryCatch({
  LDS(
    sobj = is5s_test,
    formula = ~ g3_clean + (1|hos_no) + (1|GEM),
    n_sv = NULL,
    sv_var_cutoff = 0.5,
    n_cores = 4,
    remove_na = TRUE,
    min.count = 5,        # 완화된 필터링
    min.total.count = 10,
    min.prop = 0.05,      # 5%로 완화
    large.n = 10,
    plot_sva_correlation = TRUE
  )
}, error = function(e) {
  cat("\n❌ 오류 발생:\n")
  cat("  메시지:", conditionMessage(e), "\n")
  traceback()
  return(NULL)
})

if (!is.null(result)) {
  cat("\n✓ LDS 실행 성공!\n\n")
  
  cat("6. 결과 확인:\n")
  cat("  - fit 클래스:", class(result$fit), "\n")
  cat("  - voom 클래스:", class(result$voom), "\n")
  cat("  - sva_obj 클래스:", class(result$sva_obj), "\n")
  
  if (!is.null(result$svs_used)) {
    cat("  - 사용된 SV 개수:", ncol(result$svs_used), "\n")
  } else {
    cat("  - 사용된 SV: 없음\n")
  }
  
  cat("  - 최종 포뮬러:", as.character(result$final_formula), "\n\n")
  
  # topTable 테스트
  cat("7. TopTable 테스트:\n")
  
  # dream의 경우 coef를 지정해야 p-value가 나옴
  # 고정 효과 변수 찾기
  fixed_vars <- all.vars(result$final_formula)
  fixed_vars <- fixed_vars[!grepl("^SV\\d+$", fixed_vars) & !grepl("^\\(1\\|", fixed_vars)]
  
  coef_name <- NULL
  if (length(fixed_vars) > 0) {
    # 첫 번째 고정 효과 변수를 coef로 사용
    coef_name <- fixed_vars[1]
    cat("  coef 사용:", coef_name, "\n")
    top_genes <- tryCatch({
      limma::topTable(result$fit, coef = coef_name, number = 10)
    }, error = function(e) {
      cat("  coef 지정 실패, coef 없이 시도\n")
      limma::topTable(result$fit, number = 10)
    })
  } else {
    # coef 없이 시도
    cat("  coef 없이 topTable 실행\n")
    top_genes <- limma::topTable(result$fit, number = 10)
  }
  
  cat("  상위 10개 유전자:\n")
  print(head(top_genes, 5))
  
  # p-value 확인
  if ("P.Value" %in% colnames(top_genes)) {
    n_valid_p <- sum(!is.na(top_genes$P.Value))
    cat("  ✓ P.Value 컬럼 존재\n")
    cat("  NA가 아닌 P.Value 개수:", n_valid_p, "/", nrow(top_genes), "\n")
    if (n_valid_p == 0) {
      cat("  ⚠ 모든 P.Value가 NA입니다. fit 객체 확인 필요\n")
      if (!is.null(result$fit$coefficients)) {
        cat("  fit 객체의 컬럼:", paste(colnames(result$fit$coefficients), collapse = ", "), "\n")
      }
    }
  } else {
    cat("  ⚠ P.Value 컬럼 없음. 컬럼:", paste(colnames(top_genes), collapse = ", "), "\n")
  }
  
  # 결과 저장
  cat("\n8. 결과 저장 중...\n")
  output_path <- "/data/user3/sobj/test_lds_result.qs"
  qs::qsave(result, output_path)
  cat("  ✓ 결과 저장:", output_path, "\n")
  
  # 전체 결과 저장 (coef 지정)
  if (!is.null(coef_name)) {
    top_all <- limma::topTable(result$fit, coef = coef_name, number = Inf)
  } else {
    top_all <- limma::topTable(result$fit, number = Inf)
  }
  output_csv <- "/data/user3/sobj/test_lds_topTable.csv"
  write.csv(top_all, output_csv, row.names = TRUE)
  cat("  ✓ TopTable 저장:", output_csv, "\n")
  
  # SVA 상관관계 확인
  if (!is.null(result$sv_correlation)) {
    cat("\n9. SVA 상관관계 확인:\n")
    cat("  상관관계 행렬 크기:", nrow(result$sv_correlation$correlation_matrix), "x", 
        ncol(result$sv_correlation$correlation_matrix), "\n")
    cat("  메타데이터 변수:", length(result$sv_correlation$metadata_vars), "개\n")
    cat("  상위 상관관계 (절댓값 기준):\n")
    cor_mat <- result$sv_correlation$correlation_matrix
    if (length(cor_mat) > 0 && !all(is.na(cor_mat))) {
      abs_cor <- abs(cor_mat)
      max_idx <- which(abs_cor == max(abs_cor, na.rm = TRUE), arr.ind = TRUE)
      if (nrow(max_idx) > 0) {
        for (i in 1:min(5, nrow(max_idx))) {
          sv_name <- rownames(cor_mat)[max_idx[i, 1]]
          meta_name <- colnames(cor_mat)[max_idx[i, 2]]
          cor_val <- cor_mat[max_idx[i, 1], max_idx[i, 2]]
          cat(sprintf("    %s <-> %s: %.3f\n", sv_name, meta_name, cor_val))
        }
      }
    }
  }
  
  cat("\n=== 테스트 완료 ===\n")
} else {
  cat("\n=== 테스트 실패 ===\n")
}

