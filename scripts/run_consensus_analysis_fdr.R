# ============================================================================
# Multi-Model DEG Consensus Analysis - Full Pipeline with FDR Control
# ============================================================================
# FDR threshold를 조정 가능하게 한 전체 파이프라인 실행 스크립트
# ============================================================================

# --- 1. Setup ---
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")

# DEG 방법론 함수들 로드
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")

# Consensus 분석 함수들 로드
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

library(qs)

# --- 2. 데이터 로드 ---
cat("데이터 로드 중...\n")
is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# --- 3. 파라미터 설정 ---
cluster_key <- if ("anno3.scvi" %in% colnames(is5@meta.data)) "anno3.scvi" else "seurat_clusters"
sample_key <- "hos_no"
group_key <- "g3"
batch_key <- if ("GEM" %in% colnames(is5@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

# 실행할 방법론
methods_to_run <- c(
  "muscat-edgeR", "muscat-DESeq2", "muscat-limma-voom", "muscat-limma-trend",
  "limma-voom", "edgeR-LRT"
)

# FDR threshold (완화된 값 사용)
fdr_threshold <- 0.1  # 기본값: 0.1 (0.05보다 완화)

# Consensus 필터링 파라미터
agreement_threshold <- 0.3  # Agreement score 최소값 (0~1)
min_methods <- 2  # 최소 몇 개의 방법론에서 유의해야 하는지

# --- 4. DEG 방법론 실행 ---
cat("\n=== DEG 방법론 실행 ===\n")
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

cat(sprintf("\n성공한 방법론: %d개\n", length(result_consensus$methods_run)))
cat(sprintf("실패한 방법론: %d개\n", length(result_consensus$methods_failed)))

if (length(result_consensus$methods_failed) > 0) {
  cat("\n실패한 방법론:\n")
  for (m in result_consensus$methods_failed) {
    cat(sprintf("  - %s: %s\n", m, result_consensus$errors[[m]]))
  }
}

# --- 5. 결과 표준화 ---
cat("\n=== 결과 표준화 ===\n")
standardized_results <- lapply(result_consensus$methods_run, function(m) {
  standardize_deg_results(result_consensus$results[[m]], m)
})
names(standardized_results) <- result_consensus$methods_run

# 각 방법론별 유의한 유전자 수 확인
cat("\n각 방법론별 유의한 유전자 수 (FDR < ", fdr_threshold, "):\n", sep="")
for (method in names(standardized_results)) {
  n_sig <- sum(standardized_results[[method]]$pvalue_adj < fdr_threshold, na.rm=TRUE)
  n_total <- nrow(standardized_results[[method]])
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", method, n_sig, n_total, 100*n_sig/n_total))
}

# --- 6. 행렬 구성 ---
cat("\n=== DEG 행렬 구성 ===\n")
deg_matrices <- build_deg_matrices(
  standardized_results, 
  genes = NULL, 
  fdr_threshold = fdr_threshold
)

cat(sprintf("행렬 크기: %d genes × %d methods\n", 
            nrow(deg_matrices$significance), ncol(deg_matrices$significance)))

# Significance 행렬 요약
cat("\nSignificance 행렬 요약:\n")
for (method in colnames(deg_matrices$significance)) {
  n_sig <- sum(deg_matrices$significance[, method], na.rm=TRUE)
  cat(sprintf("  %s: %d genes\n", method, n_sig))
}

# --- 7. Consensus 분석 ---
cat("\n=== Consensus 분석 ===\n")

# Agreement scores
agreement_scores <- compute_agreement_scores(deg_matrices$significance)
cat(sprintf("Agreement score 범위: %.3f ~ %.3f\n", 
            min(agreement_scores, na.rm=TRUE), max(agreement_scores, na.rm=TRUE)))
cat(sprintf("Agreement score 평균: %.3f\n", mean(agreement_scores, na.rm=TRUE)))

# Consensus scores
consensus_scores <- compute_consensus_scores(deg_matrices, agreement_scores)

cat("\nn_significant 분포:\n")
print(summary(consensus_scores$n_significant))
cat(sprintf("\nn_significant >= 1: %d genes\n", 
            sum(consensus_scores$n_significant >= 1, na.rm=TRUE)))
cat(sprintf("n_significant >= 2: %d genes\n", 
            sum(consensus_scores$n_significant >= 2, na.rm=TRUE)))
cat(sprintf("n_significant >= 3: %d genes\n", 
            sum(consensus_scores$n_significant >= 3, na.rm=TRUE)))

# --- 8. Consensus DEG 리스트 생성 ---
cat("\n=== Consensus DEG 리스트 생성 ===\n")
cat(sprintf("필터링 조건:\n"))
cat(sprintf("  - FDR threshold: %.2f\n", fdr_threshold))
cat(sprintf("  - Agreement threshold: %.2f\n", agreement_threshold))
cat(sprintf("  - min_methods: %d\n", min_methods))

consensus_deg_list <- generate_consensus_deg_list(
  consensus_scores,
  fdr_threshold = fdr_threshold,
  agreement_threshold = agreement_threshold,
  min_methods = min_methods
)

cat(sprintf("\n최종 Consensus DEG 수: %d genes\n", nrow(consensus_deg_list)))

if (nrow(consensus_deg_list) > 0) {
  cat("\n상위 20개 Consensus DEG:\n")
  print(head(consensus_deg_list[, c("gene", "n_significant", "agreement", 
                                     "consensus_score", "weighted_avg_beta")], 20))
}

# --- 9. 방법론 클러스터링 (선택적) ---
cat("\n=== 방법론 클러스터링 ===\n")
method_clustering <- cluster_deg_methods(deg_matrices, method = "hierarchical")
cat("클러스터링 완료\n")

# --- 10. 결과 저장 ---
cat("\n=== 결과 저장 ===\n")
final_result <- list(
  raw_results = result_consensus,
  standardized_results = standardized_results,
  deg_matrices = deg_matrices,
  agreement_scores = agreement_scores,
  consensus_scores = consensus_scores,
  consensus_deg_list = consensus_deg_list,
  method_clustering = method_clustering,
  parameters = list(
    fdr_threshold = fdr_threshold,
    agreement_threshold = agreement_threshold,
    min_methods = min_methods,
    methods_run = result_consensus$methods_run
  )
)

output_file <- sprintf("/data/user3/sobj/deg_consensus_fdr%.0f_min%d.qs", 
                       fdr_threshold * 100, min_methods)
qs::qsave(final_result, output_file)
cat(sprintf("결과 저장: %s\n", output_file))

cat("\n=== 완료 ===\n")

