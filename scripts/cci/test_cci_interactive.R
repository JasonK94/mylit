# ============================================================================
# CCI 분석 도구 대화형 테스트 스크립트
# R 세션에서 직접 실행: source("scripts/cci/test_cci_interactive.R")
# ============================================================================

# 환경 설정
# ============================================================================
message("=== CCI Tool Interactive Test ===")
message("This script should be run in an R session, not via Rscript")
message("")

# 필요한 패키지
library(Seurat)
library(dplyr)
library(qs)

# 함수 로드
message("Loading CCI functions...")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")

# run_nichenet_analysis 로드
cci_core_worktree <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
cci_core_mainrepo <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
  message("✓ run_nichenet_analysis loaded (worktree)")
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
  message("✓ run_nichenet_analysis loaded (main repo)")
} else {
  warning("CCI.R not found. Please source it manually.")
}

# ============================================================================
# 데이터 로드
# ============================================================================
message("\nLoading test data...")
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 메타데이터 확인
message("Checking metadata...")
colnames(sobj_test@meta.data)
table(sobj_test@meta.data$anno3.scvi, useNA = "ifany")
table(sobj_test@meta.data$g3, useNA = "ifany")

# 클러스터 목록
clusters <- unique(sobj_test@meta.data$anno3.scvi)
clusters <- clusters[!is.na(clusters)]
message("Available clusters: ", paste(clusters, collapse = ", "))

# ============================================================================
# DEG 리스트 준비
# ============================================================================
message("\nPreparing DEG list...")
receiver_cluster_test <- clusters[1]

# 예시 DEG 데이터프레임
deg_df_example <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
  cluster = rep(receiver_cluster_test, 5),
  avg_log2FC = c(1.5, 2.0, -1.2, 0.8, -0.5),
  p_val_adj = c(0.001, 0.0001, 0.01, 0.05, 0.1),
  stringsAsFactors = FALSE
)

message("Example DEG list prepared with ", nrow(deg_df_example), " genes")

# ============================================================================
# 테스트 1: 기본 기능 테스트
# ============================================================================
message("\n=== Test 1: Basic CCI Analysis ===")
message("Running run_cci_analysis...")
message("(This may take several minutes)")

results_test1 <- run_cci_analysis(
  sobj = sobj_test,
  cluster_col = "anno3.scvi",
  deg_df = deg_df_example,
  receiver_cluster = receiver_cluster_test,
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  verbose = TRUE,
  auto_save = TRUE
)

message("\nTest 1 completed!")
message("Results structure:")
str(results_test1, max.level = 2)

# 결과 확인
message("\nKey results:")
message("- Top ligands: ", length(results_test1$nichenet_results$best_upstream_ligands))
message("- Sender clusters: ", paste(results_test1$sender_clusters, collapse = ", "))
message("- Receiver cluster: ", results_test1$receiver_cluster)

if (!is.null(results_test1$saved_path)) {
  message("- Results saved to: ", results_test1$saved_path)
}

message("\n=== Test Complete ===")

