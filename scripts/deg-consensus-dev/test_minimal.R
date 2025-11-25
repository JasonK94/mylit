# ============================================================================
# 최소 테스트: 하나의 방법론만 실행하여 기본 동작 확인
# ============================================================================

cat("=== 최소 테스트 시작 ===\n\n")

# 환경 확인
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs 패키지가 필요합니다.")
}

# 함수 로드
cat("1. 함수 로드 중...\n")
tryCatch({
  source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")
  cat("   ✓ test_analysis.R\n")
}, error = function(e) {
  cat("   ✗ test_analysis.R:", conditionMessage(e), "\n")
  stop("함수 로드 실패")
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
  cat("   ✓ deg_methods_limma.R\n")
}, error = function(e) {
  cat("   ✗ deg_methods_limma.R:", conditionMessage(e), "\n")
  stop("함수 로드 실패")
})

tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
  cat("   ✓ deg_standardize.R\n")
}, error = function(e) {
  cat("   ✗ deg_standardize.R:", conditionMessage(e), "\n")
  stop("함수 로드 실패")
})

# 데이터 로드
cat("\n2. 데이터 로드 중...\n")
tryCatch({
  is5s <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
  cat(sprintf("   ✓ 데이터 로드: %d cells, %d genes\n", ncol(is5s), nrow(is5s)))
}, error = function(e) {
  cat("   ✗ 데이터 로드 실패:", conditionMessage(e), "\n")
  stop("데이터 로드 실패")
})

# 메타데이터 확인
cat("\n3. 메타데이터 확인...\n")
if ("anno3.scvi" %in% colnames(is5s@meta.data)) {
  cluster_key <- "anno3.scvi"
} else if ("seurat_clusters" %in% colnames(is5s@meta.data)) {
  cluster_key <- "seurat_clusters"
} else {
  stop("cluster_id 컬럼을 찾을 수 없습니다.")
}

sample_key <- "hos_no"
group_key <- "g3"
batch_key <- if ("GEM" %in% colnames(is5s@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

cat(sprintf("   cluster_id: %s\n", cluster_key))
cat(sprintf("   sample_id: %s\n", sample_key))
cat(sprintf("   group_id: %s\n", group_key))
cat(sprintf("   contrast: %s\n", contrast_str))

# 테스트: limma-voom 하나만 실행
cat("\n4. limma-voom 테스트 실행...\n")
result_limma <- tryCatch({
  runLIMMA_voom_v1(
    sobj = is5s,
    cluster_id = cluster_key,
    sample_id = sample_key,
    group_id = group_key,
    batch_id = batch_key,
    contrast = contrast_str,
    remove_na_groups = TRUE
  )
}, error = function(e) {
  cat("   ✗ limma-voom 실패:", conditionMessage(e), "\n")
  traceback()
  NULL
})

if (!is.null(result_limma)) {
  cat(sprintf("   ✓ limma-voom 완료: %d rows\n", nrow(result_limma)))
  cat(sprintf("   컬럼: %s\n", paste(colnames(result_limma), collapse=", ")))
  
  # 결과 저장
  qs::qsave(result_limma, "/data/user3/sobj/test_minimal_limma_voom.qs")
  cat("   ✓ 결과 저장 완료\n")
  
  # 표준화 테스트
  cat("\n5. 결과 표준화 테스트...\n")
  standardized <- tryCatch({
    standardize_deg_results(result_limma, "limma-voom")
  }, error = function(e) {
    cat("   ✗ 표준화 실패:", conditionMessage(e), "\n")
    traceback()
    NULL
  })
  
  if (!is.null(standardized)) {
    cat(sprintf("   ✓ 표준화 완료: %d rows\n", nrow(standardized)))
    cat(sprintf("   컬럼: %s\n", paste(colnames(standardized), collapse=", ")))
    cat(sprintf("   상위 5개 유전자:\n"))
    print(head(standardized[, c("gene", "logFC", "pvalue", "pvalue_adj")], 5))
  }
} else {
  stop("limma-voom 실행 실패")
}

cat("\n=== 최소 테스트 완료 ===\n")

