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
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
if (!exists("LDS")) {
  stop("LDS 함수를 찾을 수 없습니다.")
}
cat("✓ LDS 함수 확인 완료\n\n")

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
cat("   SV: 자동 결정 (sv_var_cutoff = 0.5)\n\n")

result <- tryCatch({
  LDS(
    sobj = is5s_test,
    formula = ~ g3_clean + (1|hos_no) + (1|GEM),
    n_sv = NULL,
    sv_var_cutoff = 0.5,
    n_cores = 4
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
  top_genes <- limma::topTable(result$fit, number = 10)
  cat("  상위 10개 유전자:\n")
  print(head(top_genes, 5))
  
  # 결과 저장
  cat("\n8. 결과 저장 중...\n")
  output_path <- "/data/user3/sobj/test_lds_result.qs"
  qs::qsave(result, output_path)
  cat("  ✓ 결과 저장:", output_path, "\n")
  
  top_all <- limma::topTable(result$fit, number = Inf)
  output_csv <- "/data/user3/sobj/test_lds_topTable.csv"
  write.csv(top_all, output_csv, row.names = TRUE)
  cat("  ✓ TopTable 저장:", output_csv, "\n")
  
  cat("\n=== 테스트 완료 ===\n")
} else {
  cat("\n=== 테스트 실패 ===\n")
}

