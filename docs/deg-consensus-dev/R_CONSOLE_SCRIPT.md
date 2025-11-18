# R Console에서 직접 실행 가능한 스크립트

## 현재 브랜치
- **브랜치**: `deg-consensus-dev`

## 실행 스크립트

### 전체 파이프라인 실행 (is5가 이미 로드되어 있다고 가정)

```r
# 1. 함수 로드 (한 번만 실행)
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 2. 메타데이터 설정
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
deg_matrices <- build_deg_matrices(standardized_results, genes = NULL, fdr_threshold = 0.05)

# 6. Consensus 분석
agreement_scores <- compute_agreement_scores(deg_matrices$significance)
consensus_scores <- compute_consensus_scores(deg_matrices, agreement_scores)
consensus_deg_list <- generate_consensus_deg_list(
  consensus_scores,
  agreement_threshold = 0.5,
  min_methods = ceiling(length(result_consensus$methods_run) * 0.5)
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
```

### 또는 스크립트 파일 실행

```r
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/deg-consensus/run_consensus_simple.R")
```

