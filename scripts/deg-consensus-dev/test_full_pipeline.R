# ============================================================================
# 전체 파이프라인 테스트 스크립트
# ============================================================================
# Phase 1-6 전체 테스트
# ============================================================================

# 환경 설정
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs 패키지가 필요합니다.")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat 패키지가 필요합니다.")
}

# 함수 로드
cat("=== 함수 로드 ===\n")
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")  # runMUSCAT2_v1, runNEBULA2_v1
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 데이터 로드
cat("\n=== 데이터 로드 ===\n")
is5s <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
cat(sprintf("셀 수: %d\n", ncol(is5s)))
cat(sprintf("유전자 수: %d\n", nrow(is5s)))

# 메타데이터 확인
cat("\n=== 메타데이터 확인 ===\n")
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

cat(sprintf("cluster_id: %s\n", cluster_key))
cat(sprintf("sample_id: %s\n", sample_key))
cat(sprintf("group_id: %s\n", group_key))
cat(sprintf("batch_id: %s\n", ifelse(is.null(batch_key), "NULL", batch_key)))
cat(sprintf("contrast: %s\n", contrast_str))

# ============================================================================
# Phase 1-2: 여러 방법론 실행
# ============================================================================
cat("\n=== Phase 1-2: 여러 방법론 실행 ===\n")
# 처음에는 적은 수의 방법론으로 테스트
methods_to_test <- c(
  "muscat-edgeR",
  "limma-voom",
  "edgeR-LRT"
)
# 성공하면 추가
# "muscat-DESeq2",
# "limma-trend",
# "edgeR-QLF",
# "DESeq2-Wald"

result_consensus <- tryCatch({
  run_deg_consensus(
    sobj = is5s,
    contrast = contrast_str,
    methods = methods_to_test,
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

if (is.null(result_consensus)) {
  stop("run_deg_consensus 실패")
}

cat(sprintf("✓ 성공한 방법론: %s\n", paste(result_consensus$methods_run, collapse=", ")))
if (length(result_consensus$methods_failed) > 0) {
  cat(sprintf("✗ 실패한 방법론: %s\n", paste(result_consensus$methods_failed, collapse=", ")))
}

# 결과 저장
qs::qsave(result_consensus, "/data/user3/sobj/test_deg_consensus_full_pipeline.qs")
cat("✓ Phase 1-2 결과 저장 완료\n")

# ============================================================================
# Phase 3: 결과 표준화 및 행렬 구성
# ============================================================================
cat("\n=== Phase 3: 결과 표준화 및 행렬 구성 ===\n")

standardized_results <- list()
for (method in result_consensus$methods_run) {
  cat(sprintf("표준화 중: %s\n", method))
  tryCatch({
    standardized_results[[method]] <- standardize_deg_results(
      deg_result = result_consensus$results[[method]],
      method_name = method
    )
    cat(sprintf("  ✓ %d 행\n", nrow(standardized_results[[method]])))
  }, error = function(e) {
    cat(sprintf("  ✗ 실패: %s\n", conditionMessage(e)))
  })
}

if (length(standardized_results) == 0) {
  stop("표준화된 결과가 없습니다.")
}

cat(sprintf("✓ %d 개 방법론 표준화 완료\n", length(standardized_results)))

# 행렬 구성
deg_matrices <- tryCatch({
  build_deg_matrices(
    standardized_results_list = standardized_results,
    genes = NULL,
    fdr_threshold = 0.05
  )
}, error = function(e) {
  cat(sprintf("✗ 행렬 구성 실패: %s\n", conditionMessage(e)))
  traceback()
  NULL
})

if (is.null(deg_matrices)) {
  stop("행렬 구성 실패")
}

cat(sprintf("✓ 행렬 구성 완료: %d 유전자 × %d 방법론\n", 
            nrow(deg_matrices$beta), ncol(deg_matrices$beta)))

# 결과 저장
qs::qsave(list(
  standardized_results = standardized_results,
  deg_matrices = deg_matrices
), "/data/user3/sobj/test_deg_consensus_phase3.qs")
cat("✓ Phase 3 결과 저장 완료\n")

# ============================================================================
# Phase 4-5: Consensus 분석
# ============================================================================
cat("\n=== Phase 4-5: Consensus 분석 ===\n")

# Agreement scores
agreement_scores <- compute_agreement_scores(deg_matrices$significance)
cat(sprintf("✓ Agreement scores 계산 완료: %d 유전자\n", length(agreement_scores)))

# PCA
pca_result <- tryCatch({
  perform_deg_pca(deg_matrices)
}, error = function(e) {
  cat(sprintf("✗ PCA 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(pca_result)) {
  cat(sprintf("✓ PCA 완료: %d components\n", ncol(pca_result$coordinates)))
}

# Clustering
clustering_result <- tryCatch({
  cluster_deg_methods(deg_matrices, method = "hierarchical")
}, error = function(e) {
  cat(sprintf("✗ Clustering 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(clustering_result)) {
  cat(sprintf("✓ Clustering 완료: %d clusters\n", clustering_result$k))
  cat(sprintf("  클러스터 할당: %s\n", paste(clustering_result$clusters, collapse=", ")))
}

# Consensus scores
consensus_scores <- tryCatch({
  compute_consensus_scores(
    deg_matrices = deg_matrices,
    agreement_scores = agreement_scores,
    clustering_results = clustering_result
  )
}, error = function(e) {
  cat(sprintf("✗ Consensus scores 계산 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(consensus_scores)) {
  cat(sprintf("✓ Consensus scores 계산 완료: %d 유전자\n", nrow(consensus_scores)))
}

# Consensus DEG list
consensus_deg_list <- tryCatch({
  generate_consensus_deg_list(
    consensus_scores = consensus_scores,
    fdr_threshold = 0.05,
    agreement_threshold = 0.5,
    min_methods = ceiling(length(result_consensus$methods_run) * 0.5)
  )
}, error = function(e) {
  cat(sprintf("✗ Consensus DEG list 생성 실패: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(consensus_deg_list)) {
  cat(sprintf("✓ Consensus DEG list 생성 완료: %d 유전자\n", nrow(consensus_deg_list)))
  cat(sprintf("  상위 10개 유전자:\n"))
  print(head(consensus_deg_list[, c("gene", "agreement", "consensus_score", "n_significant")], 10))
}

# 최종 결과 저장
final_result <- list(
  raw_results = result_consensus,
  standardized_results = standardized_results,
  deg_matrices = deg_matrices,
  agreement_scores = agreement_scores,
  pca_result = pca_result,
  clustering_result = clustering_result,
  consensus_scores = consensus_scores,
  consensus_deg_list = consensus_deg_list
)

qs::qsave(final_result, "/data/user3/sobj/test_deg_consensus_final_result.qs")
cat("\n✓ 최종 결과 저장 완료: /data/user3/sobj/test_deg_consensus_final_result.qs\n")

cat("\n=== 전체 파이프라인 테스트 완료 ===\n")

