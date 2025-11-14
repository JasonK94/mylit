# ============================================================================
# Multi-Model DEG Consensus Engine - 실행 스크립트
# ============================================================================
# R console에서 직접 실행 가능한 스크립트
# is5 (또는 다른 Seurat 객체)가 이미 로드되어 있다고 가정
# ============================================================================

# 함수 로드 (한 번만 실행)
cat("=== 함수 로드 ===\n")
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")
cat("✓ 함수 로드 완료\n\n")

# ============================================================================
# 메타데이터 확인 및 설정
# ============================================================================
# is5 (또는 다른 Seurat 객체)가 이미 로드되어 있다고 가정
# 필요시 아래 주석을 해제하여 데이터 로드:
# is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

cat("=== 메타데이터 확인 ===\n")
# cluster_id 자동 감지
if ("anno3.scvi" %in% colnames(is5@meta.data)) {
  cluster_key <- "anno3.scvi"
} else if ("seurat_clusters" %in% colnames(is5@meta.data)) {
  cluster_key <- "seurat_clusters"
} else {
  stop("cluster_id 컬럼을 찾을 수 없습니다.")
}

sample_key <- "hos_no"
group_key <- "g3"
batch_key <- if ("GEM" %in% colnames(is5@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

cat(sprintf("cluster_id: %s\n", cluster_key))
cat(sprintf("sample_id: %s\n", sample_key))
cat(sprintf("group_id: %s\n", group_key))
cat(sprintf("batch_id: %s\n", ifelse(is.null(batch_key), "NULL", batch_key)))
cat(sprintf("contrast: %s\n", contrast_str))
cat(sprintf("셀 수: %d\n", ncol(is5)))
cat(sprintf("유전자 수: %d\n", nrow(is5)))
cat("\n")

# ============================================================================
# 여러 방법론 실행
# ============================================================================
cat("=== 여러 방법론 실행 ===\n")

methods_to_run <- c(
  "muscat-edgeR",
  "muscat-DESeq2",
  "muscat-limma-voom",
  "muscat-limma-trend",
  "limma-voom",
  "limma-trend",
  "edgeR-LRT",
  "edgeR-QLF",
  "DESeq2-Wald",
  "DESeq2-LRT"
)

result_consensus <- run_deg_consensus(
  sobj = is5,
  contrast = contrast_str,
  methods = methods_to_run,
  cluster_id = cluster_key,
  sample_id = sample_key,
  group_id = group_key,
  batch_id = batch_key,
  remove_na_groups = TRUE,
  verbose = TRUE
)

cat(sprintf("\n✓ 성공: %d 개 방법론\n", length(result_consensus$methods_run)))
cat(sprintf("✗ 실패: %d 개 방법론\n", length(result_consensus$methods_failed)))
if (length(result_consensus$methods_failed) > 0) {
  cat(sprintf("실패한 방법론: %s\n", paste(result_consensus$methods_failed, collapse=", ")))
}
cat("\n")

# ============================================================================
# 결과 표준화 및 행렬 구성
# ============================================================================
cat("=== 결과 표준화 및 행렬 구성 ===\n")

standardized_results <- list()
for (method in result_consensus$methods_run) {
  cat(sprintf("표준화 중: %s\n", method))
  standardized_results[[method]] <- standardize_deg_results(
    deg_result = result_consensus$results[[method]],
    method_name = method
  )
  cat(sprintf("  ✓ %d rows\n", nrow(standardized_results[[method]])))
}

deg_matrices <- build_deg_matrices(
  standardized_results_list = standardized_results,
  genes = NULL,
  fdr_threshold = 0.05
)

cat(sprintf("\n✓ 행렬 구성 완료: %d genes × %d methods\n", 
            nrow(deg_matrices$beta), ncol(deg_matrices$beta)))
cat("\n")

# ============================================================================
# Consensus 분석
# ============================================================================
cat("=== Consensus 분석 ===\n")

# Agreement scores
agreement_scores <- compute_agreement_scores(deg_matrices$significance)
cat(sprintf("✓ Agreement scores: %d genes\n", length(agreement_scores)))
cat(sprintf("  범위: %.3f ~ %.3f\n", min(agreement_scores, na.rm=TRUE), max(agreement_scores, na.rm=TRUE)))
cat(sprintf("  평균: %.3f\n", mean(agreement_scores, na.rm=TRUE)))

# PCA (선택적)
pca_result <- tryCatch({
  perform_deg_pca(deg_matrices)
}, error = function(e) {
  cat(sprintf("PCA 실패: %s\n", conditionMessage(e)))
  NULL
})

# Clustering (선택적)
clustering_result <- tryCatch({
  cluster_deg_methods(deg_matrices, method = "hierarchical")
}, error = function(e) {
  cat(sprintf("Clustering 실패: %s\n", conditionMessage(e)))
  NULL
})

# Consensus scores
consensus_scores <- compute_consensus_scores(
  deg_matrices = deg_matrices,
  agreement_scores = agreement_scores,
  clustering_results = clustering_result
)

cat(sprintf("✓ Consensus scores: %d genes\n", nrow(consensus_scores)))

# Consensus DEG list
consensus_deg_list <- generate_consensus_deg_list(
  consensus_scores = consensus_scores,
  fdr_threshold = 0.05,
  agreement_threshold = 0.5,
  min_methods = ceiling(length(result_consensus$methods_run) * 0.5)
)

cat(sprintf("✓ Consensus DEG list: %d genes\n", nrow(consensus_deg_list)))

if (nrow(consensus_deg_list) > 0) {
  cat("\n상위 20개 유전자:\n")
  print(head(consensus_deg_list[, c("gene", "agreement", "consensus_score", "n_significant", "mean_beta")], 20))
}
cat("\n")

# ============================================================================
# 결과 저장
# ============================================================================
cat("=== 결과 저장 ===\n")

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

output_file <- "/data/user3/sobj/deg_consensus_final_result.qs"
qs::qsave(final_result, output_file)
cat(sprintf("✓ 결과 저장: %s\n", output_file))

cat("\n=== 분석 완료 ===\n")

