# ============================================================================
# Phase 2 테스트 스크립트: limma 계열 방법론
# ============================================================================
# is5s 데이터로 limma-voom, limma-trend 테스트
# ============================================================================

# 환경 설정
library(qs)
library(Seurat)

# 함수 로드
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")  # runMUSCAT2_v1, runNEBULA2_v1
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")  # limma 함수들
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")  # run_deg_consensus

# 데이터 로드
message("=== 데이터 로드 ===")
is5s <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
message(sprintf("셀 수: %d", ncol(is5s)))
message(sprintf("유전자 수: %d", nrow(is5s)))

# 메타데이터 확인
message("\n=== 메타데이터 확인 ===")
# 먼저 사용 가능한 컬럼 확인
message("사용 가능한 메타데이터 컬럼:")
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
batch_key <- "GEM"
contrast_str <- "2 - 1"  # g3==2 vs g3==1

message(sprintf("cluster_id: %s", cluster_key))
message(sprintf("cluster levels: %d", length(unique(is5s@meta.data[[cluster_key]]))))
message(sprintf("sample_id (hos_no): %d samples", length(unique(is5s@meta.data[[sample_key]]))))
message(sprintf("group_id (g3): %s", paste(unique(is5s@meta.data[[group_key]]), collapse=", ")))
if (batch_key %in% colnames(is5s@meta.data)) {
  message(sprintf("batch_id (GEM): %d levels", length(unique(is5s@meta.data[[batch_key]]))))
} else {
  message("batch_id (GEM): 없음 - NULL로 설정")
  batch_key <- NULL
}

# ============================================================================
# 테스트 1: limma-voom 직접 테스트
# ============================================================================
message("\n=== 테스트 1: limma-voom 직접 실행 ===")
tryCatch({
  result_limma_voom <- runLIMMA_voom_v1(
    sobj = is5s,
    cluster_id = cluster_key,
    sample_id = sample_key,
    group_id = group_key,
    batch_id = batch_key,
    contrast = contrast_str,
    remove_na_groups = TRUE
  )
  
  message(sprintf("✓ limma-voom 완료: %d 행", nrow(result_limma_voom)))
  message(sprintf("  컬럼: %s", paste(colnames(result_limma_voom), collapse=", ")))
  message(sprintf("  클러스터 수: %d", length(unique(result_limma_voom$cluster_id))))
  
  # 결과 저장
  qs::qsave(result_limma_voom, "/data/user3/sobj/test_limma_voom_v1_result.qs")
  message("✓ 결과 저장: /data/user3/sobj/test_limma_voom_v1_result.qs")
  
}, error = function(e) {
  message(sprintf("✗ limma-voom 실패: %s", conditionMessage(e)))
})

# ============================================================================
# 테스트 2: limma-trend 직접 테스트
# ============================================================================
message("\n=== 테스트 2: limma-trend 직접 실행 ===")
tryCatch({
  result_limma_trend <- runLIMMA_trend_v1(
    sobj = is5s,
    cluster_id = cluster_key,
    sample_id = sample_key,
    group_id = group_key,
    batch_id = batch_key,
    contrast = contrast_str,
    remove_na_groups = TRUE
  )
  
  message(sprintf("✓ limma-trend 완료: %d 행", nrow(result_limma_trend)))
  message(sprintf("  컬럼: %s", paste(colnames(result_limma_trend), collapse=", ")))
  message(sprintf("  클러스터 수: %d", length(unique(result_limma_trend$cluster_id))))
  
  # 결과 저장
  qs::qsave(result_limma_trend, "/data/user3/sobj/test_limma_trend_v1_result.qs")
  message("✓ 결과 저장: /data/user3/sobj/test_limma_trend_v1_result.qs")
  
}, error = function(e) {
  message(sprintf("✗ limma-trend 실패: %s", conditionMessage(e)))
})

# ============================================================================
# 테스트 3: run_deg_consensus로 통합 테스트
# ============================================================================
message("\n=== 테스트 3: run_deg_consensus 통합 테스트 ===")
tryCatch({
  result_consensus <- run_deg_consensus(
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
  
  message(sprintf("✓ run_deg_consensus 완료"))
  message(sprintf("  성공한 방법론: %s", paste(result_consensus$methods_run, collapse=", ")))
  if (length(result_consensus$methods_failed) > 0) {
    message(sprintf("  실패한 방법론: %s", paste(result_consensus$methods_failed, collapse=", ")))
  }
  
  # 결과 저장
  qs::qsave(result_consensus, "/data/user3/sobj/test_deg_consensus_phase2_limma.qs")
  message("✓ 결과 저장: /data/user3/sobj/test_deg_consensus_phase2_limma.qs")
  
}, error = function(e) {
  message(sprintf("✗ run_deg_consensus 실패: %s", conditionMessage(e)))
})

message("\n=== 모든 테스트 완료 ===")

