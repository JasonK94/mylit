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
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_base.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_dream.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 1-1. 출력 파일 prefix 설정 (데이터 크기 + 타임스탬프 기반)
out_dir <- "/data/user3/sobj"
cons_dir <- file.path(out_dir, "consensus")
if (!dir.exists(cons_dir)) {
  dir.create(cons_dir, recursive = TRUE, showWarnings = FALSE)
}
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
n_cells_is5 <- tryCatch(ncol(is5), error = function(e) NA_integer_)
if (!exists("dataset_tag")) {
  dataset_tag <- if (!is.na(n_cells_is5) && n_cells_is5 <= 3000) "ds" else "full"
}
file_prefix <- file.path(cons_dir, paste0("deg_consensus_", dataset_tag, "_", timestamp))

# 2. 메타데이터 설정 (is5가 이미 로드되어 있다고 가정)
cluster_key <- if ("anno3.scvi" %in% colnames(is5@meta.data)) "anno3.scvi" else "seurat_clusters"
sample_key <- "hos_no"
group_key <- "g3"
batch_key <- if ("GEM" %in% colnames(is5@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

# 3. 여러 방법론 실행
base_methods <- c(
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

# 고급(선택) 방법론 토글
include_dream <- TRUE
include_nebula <- FALSE
include_nebula_pb <- FALSE

optional_methods <- c(
  if (isTRUE(include_dream)) "dream",
  if (isTRUE(include_nebula)) "nebula",
  if (isTRUE(include_nebula_pb)) "nebula-pb"
)

methods_to_run <- unique(c(base_methods, optional_methods))

result_consensus <- run_deg_consensus(
  sobj = is5,
  contrast = contrast_str,
  methods = methods_to_run,
  cluster_id = cluster_key,
  sample_id = sample_key,
  group_id = group_key,
  batch_id = batch_key,
  covar_effects = "sex",  # sex를 fixed effect로 추가
  remove_na_groups = TRUE,
  verbose = TRUE
)

# 중간 결과 저장: raw method-level 결과
qs::qsave(result_consensus, paste0(file_prefix, "_raw_results.qs"))

# 4. 결과 표준화
standardized_results <- lapply(result_consensus$methods_run, function(m) {
  standardize_deg_results(result_consensus$results[[m]], m)
})
names(standardized_results) <- result_consensus$methods_run

# 중간 결과 저장: 표준화된 결과
qs::qsave(standardized_results, paste0(file_prefix, "_standardized_results.qs"))

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

# 중간 결과 저장: gene × method 행렬
qs::qsave(deg_matrices, paste0(file_prefix, "_deg_matrices.qs"))

# 6. Consensus 분석
agreement_scores <- compute_agreement_scores(deg_matrices$significance)
consensus_scores <- compute_consensus_scores(deg_matrices, agreement_scores)

# 중간 결과 저장: agreement 및 consensus score
qs::qsave(agreement_scores, paste0(file_prefix, "_agreement_scores.qs"))
qs::qsave(consensus_scores, paste0(file_prefix, "_consensus_scores.qs"))

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

# 중간 결과 저장: 최종 consensus DEG 리스트만 별도 저장
qs::qsave(consensus_deg_list, paste0(file_prefix, "_consensus_deg_list.qs"))

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


qs::qsave(final_result, paste0(file_prefix, "_final_result.qs"))

