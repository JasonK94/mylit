# ============================================================================
# Phase 2 테스트 스크립트: limma 계열 방법론 (Interactive)
# ============================================================================
# 사용법: R 세션에서 source()로 실행
# ============================================================================

# 환경 설정
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs 패키지가 필요합니다.")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat 패키지가 필요합니다. devtools::load_all() 또는 library(Seurat)로 로드하세요.")
}

# 함수 로드
cat("=== 함수 로드 ===\n")
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")  # runMUSCAT2_v1, runNEBULA2_v1
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")  # limma 함수들
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")  # run_deg_consensus

# 데이터 로드
cat("\n=== 데이터 로드 ===\n")
is5s <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
cat(sprintf("셀 수: %d\n", ncol(is5s)))
cat(sprintf("유전자 수: %d\n", nrow(is5s)))

# 메타데이터 확인
cat("\n=== 메타데이터 확인 ===\n")
# 먼저 사용 가능한 컬럼 확인
cat("사용 가능한 메타데이터 컬럼 (처음 20개):\n")
print(head(colnames(is5s@meta.data), 20))

# cluster_id 확인 (여러 가능성 체크)
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
contrast_str <- "2 - 1"  # g3==2 vs g3==1

cat(sprintf("\ncluster_id: %s\n", cluster_key))
cat(sprintf("cluster levels: %d\n", length(unique(is5s@meta.data[[cluster_key]]))))
cat(sprintf("sample_id (hos_no): %d samples\n", length(unique(is5s@meta.data[[sample_key]]))))
cat(sprintf("group_id (g3): %s\n", paste(unique(is5s@meta.data[[group_key]]), collapse=", ")))
if (!is.null(batch_key)) {
  cat(sprintf("batch_id (GEM): %d levels\n", length(unique(is5s@meta.data[[batch_key]]))))
} else {
  cat("batch_id (GEM): 없음 - NULL로 설정\n")
}

# ============================================================================
# 테스트 1: limma-voom 직접 테스트
# ============================================================================
cat("\n=== 테스트 1: limma-voom 직접 실행 ===\n")
result_limma_voom <- tryCatch({
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
  cat(sprintf("✗ limma-voom 실패: %s\n", conditionMessage(e)))
  traceback()
  NULL
})

if (!is.null(result_limma_voom)) {
  cat(sprintf("✓ limma-voom 완료: %d 행\n", nrow(result_limma_voom)))
  cat(sprintf("  컬럼: %s\n", paste(colnames(result_limma_voom), collapse=", ")))
  cat(sprintf("  클러스터 수: %d\n", length(unique(result_limma_voom$cluster_id))))
  
  # 결과 저장
  qs::qsave(result_limma_voom, "/data/user3/sobj/test_limma_voom_v1_result.qs")
  cat("✓ 결과 저장: /data/user3/sobj/test_limma_voom_v1_result.qs\n")
}

# ============================================================================
# 테스트 2: limma-trend 직접 테스트
# ============================================================================
cat("\n=== 테스트 2: limma-trend 직접 실행 ===\n")
result_limma_trend <- tryCatch({
  runLIMMA_trend_v1(
    sobj = is5s,
    cluster_id = cluster_key,
    sample_id = sample_key,
    group_id = group_key,
    batch_id = batch_key,
    contrast = contrast_str,
    remove_na_groups = TRUE
  )
}, error = function(e) {
  cat(sprintf("✗ limma-trend 실패: %s\n", conditionMessage(e)))
  traceback()
  NULL
})

if (!is.null(result_limma_trend)) {
  cat(sprintf("✓ limma-trend 완료: %d 행\n", nrow(result_limma_trend)))
  cat(sprintf("  컬럼: %s\n", paste(colnames(result_limma_trend), collapse=", ")))
  cat(sprintf("  클러스터 수: %d\n", length(unique(result_limma_trend$cluster_id))))
  
  # 결과 저장
  qs::qsave(result_limma_trend, "/data/user3/sobj/test_limma_trend_v1_result.qs")
  cat("✓ 결과 저장: /data/user3/sobj/test_limma_trend_v1_result.qs\n")
}

# ============================================================================
# 테스트 3: run_deg_consensus로 통합 테스트
# ============================================================================
cat("\n=== 테스트 3: run_deg_consensus 통합 테스트 ===\n")
result_consensus <- tryCatch({
  run_deg_consensus(
    sobj = is5s,
    contrast = contrast_str,
    methods = c("limma-voom", "limma-trend"),
    cluster_id = cluster_key,
    sample_id = sample_key,
    group_id = group_key,
    batch_id = batch_key,
    remove_na_groups = TRUE,
    verbose = TRUE
  )
}, error = function(e) {
  cat(sprintf("✗ run_deg_consensus 실패: %s\n", conditionMessage(e)))
  traceback()
  NULL
})

if (!is.null(result_consensus)) {
  cat(sprintf("✓ run_deg_consensus 완료\n"))
  cat(sprintf("  성공한 방법론: %s\n", paste(result_consensus$methods_run, collapse=", ")))
  if (length(result_consensus$methods_failed) > 0) {
    cat(sprintf("  실패한 방법론: %s\n", paste(result_consensus$methods_failed, collapse=", ")))
  }
  
  # 결과 저장
  qs::qsave(result_consensus, "/data/user3/sobj/test_deg_consensus_phase2_limma.qs")
  cat("✓ 결과 저장: /data/user3/sobj/test_deg_consensus_phase2_limma.qs\n")
}

cat("\n=== 모든 테스트 완료 ===\n")

