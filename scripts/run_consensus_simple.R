# ============================================================================
# R Console에서 직접 실행 가능한 간단한 스크립트
# ============================================================================
# is5 (또는 다른 Seurat 객체)가 이미 로드되어 있다고 가정
# ============================================================================

# 1. 함수 로드 (한 번만 실행)
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 2. 메타데이터 설정 (is5가 이미 로드되어 있다고 가정)
cluster_key <- if ("anno3.scvi" %in% colnames(is5@meta.data)) "anno3.scvi" else "seurat_clusters"
sample_key <- "hos_no"
group_key <- "g3"
batch_key <- if ("GEM" %in% colnames(is5@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

# 3. 여러 방법론 실행
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

# 4. 결과 표준화
standardized_results <- lapply(result_consensus$methods_run, function(m) {
  standardize_deg_results(result_consensus$results[[m]], m)
})
names(standardized_results) <- result_consensus$methods_run

# 5. 행렬 구성
# per-method significance 설정:
# - significance_mode = "fdr": 각 방법론 내부의 adj.p (pvalue_adj)와 fdr_threshold 사용
# - significance_mode = "pvalue": raw p.value와 pvalue_threshold 사용 (조금 더 관대함)
#
# 여기서는 adj.p 대신 raw p.value 기반으로 "유의 여부"를 정의하고,
# 최종 FDR 제어는 meta-analysis p-value (meta_p_adj)에 대해 수행한다.
fdr_threshold <- 0.1        # meta-analysis FDR threshold (generate_consensus_deg_list에서 사용)
pvalue_threshold <- 0.01    # per-method raw p-value threshold (build_deg_matrices에서 사용)

deg_matrices <- build_deg_matrices(
  standardized_results,
  genes = NULL,
  fdr_threshold = fdr_threshold,          # significance_mode = "fdr"인 경우에만 사용
  significance_mode = "pvalue",          # raw p.value 기반 유의성
  pvalue_threshold = pvalue_threshold    # per-method p.value threshold
)

# 6. Consensus 분석
agreement_scores <- compute_agreement_scores(deg_matrices$significance)
consensus_scores <- compute_consensus_scores(deg_matrices, agreement_scores)

# Consensus DEG 리스트 생성
# min_methods: 최소 몇 개의 방법론에서 유의해야 하는지
#   - 예: min_methods = 2 → 최소 2개 방법론에서 유의한 유전자만 선택
#   - 예: min_methods = NULL → 방법론 수 제한 없음
# agreement_threshold: Agreement score 최소값 (0~1, 높을수록 엄격)
min_methods <- 2  # 최소 2개 방법론에서 유의
agreement_threshold <- 0.3  # Agreement score >= 0.3

consensus_deg_list <- generate_consensus_deg_list(
  consensus_scores,
  fdr_threshold = fdr_threshold,
  agreement_threshold = agreement_threshold,
  min_methods = min_methods
)

# 7. 결과 확인
cat(sprintf("\n성공한 방법론: %d 개\n", length(result_consensus$methods_run)))
cat(sprintf("Consensus DEG: %d genes\n", nrow(consensus_deg_list)))
if (nrow(consensus_deg_list) > 0) {
  print(head(consensus_deg_list[, c("gene", "agreement", "consensus_score", "n_significant")], 20))
}

# 8. 결과 저장
final_result <- list(
  raw_results = result_consensus,
  standardized_results = standardized_results,
  deg_matrices = deg_matrices,
  agreement_scores = agreement_scores,
  consensus_scores = consensus_scores,
  consensus_deg_list = consensus_deg_list
)
qs::qsave(final_result, "/data/user3/sobj/deg_consensus_final_result.qs")

