# Shared options
nichenet_data_dir_default <- "/data/user3/git_repo/human"

# ============================================================================
# CCI 분석 도구 테스트 스크립트
# ============================================================================

# 환경 설정
# ============================================================================
# 패키지 로드
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 함수 소스 직접 로드 (devtools::load_all은 패키지 구조가 필요하므로 소스 로드 사용)
# CCI 모듈 함수들
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")

# 기존 CCI.R 함수도 필요 (run_nichenet_analysis)
cci_core_worktree <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
cci_core_mainrepo <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
} else {
  warning("CCI.R not found. run_nichenet_analysis may not be available.")
}

# 필요한 패키지
library(Seurat)
library(dplyr)
library(qs)
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect)
  conflicted::conflicts_prefer(base::setdiff)
}

# ============================================================================
# 데이터 로드
# ============================================================================
message("Loading test data...")

# 테스트 데이터 (다운샘플)
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 메타데이터 확인
message("Checking metadata...")
colnames(sobj_test@meta.data)
table(sobj_test@meta.data$anno3.scvi, useNA = "ifany")
table(sobj_test@meta.data$g3, useNA = "ifany")

# 클러스터 목록 확인
clusters <- unique(sobj_test@meta.data$anno3.scvi)
clusters <- clusters[!is.na(clusters)]
message("Available clusters: ", paste(clusters, collapse = ", "))

# ============================================================================
# DEG 리스트 준비 (예시)
# ============================================================================
message("Preparing DEG list...")

# 예시 1: 간단한 테스트용 DEG 리스트
# 실제 사용 시에는 muscat, edgeR, DESeq2 등으로 분석한 결과를 사용
receiver_cluster_test <- clusters[1]  # 첫 번째 클러스터를 receiver로 사용

# receiver 클러스터 내에서 발현된 유전자 + NicheNet 타깃 매트릭스 교집합 확보
receiver_cells <- sobj_test@meta.data[["anno3.scvi"]] == receiver_cluster_test
receiver_expr <- Seurat::GetAssayData(sobj_test, assay = "RNA", slot = "counts")[, receiver_cells, drop = FALSE]
expressed_genes <- rownames(receiver_expr)[Matrix::rowSums(receiver_expr > 0) > 0]

nichenet_gene_pool <- character(0)
ligand_target_path <- file.path(nichenet_data_dir_default, "ligand_target_matrix_nsga2r_final.rds")
if (file.exists(ligand_target_path)) {
  nichenet_gene_pool <- tryCatch({
    rownames(readRDS(ligand_target_path))
  }, error = function(e) {
    message("  ! Warning: Failed to read ligand_target_matrix from ", ligand_target_path, ": ", e$message)
    character(0)
  })
}

candidate_genes <- expressed_genes
if (length(nichenet_gene_pool) > 0) {
  candidate_genes <- base::intersect(expressed_genes, nichenet_gene_pool)
}
if (length(candidate_genes) == 0) {
  candidate_genes <- expressed_genes
}

deg_genes <- head(candidate_genes, n = 5)
if (length(deg_genes) < 3) {
  warning("Fewer than 3 overlap genes found between receiver expression and NicheNet target matrix; tests may fail.")
}

# 예시 DEG 데이터프레임 생성 (실제로는 분석 결과를 사용)
deg_df_example <- data.frame(
  gene = deg_genes,
  cluster = rep(receiver_cluster_test, length(deg_genes)),
  avg_log2FC = rep(1.5, length(deg_genes)),
  p_val_adj = rep(0.01, length(deg_genes)),
  stringsAsFactors = FALSE
)

message("Example DEG list prepared with ", nrow(deg_df_example), " genes")

# ============================================================================
# 테스트 1: 기본 기능 테스트
# ============================================================================
message("\n=== Test 1: Basic CCI Analysis ===")

tryCatch({
  results_test1 <- run_cci_analysis(
    sobj = sobj_test,
    cluster_col = "anno3.scvi",
    deg_df = deg_df_example,
    receiver_cluster = receiver_cluster_test,
    condition_col = "g3",
    condition_oi = "2",
    condition_ref = "1",
    species = "human",
    nichenet_data_dir = nichenet_data_dir_default,
    verbose = TRUE,
    auto_save = TRUE
  )
  
  message("Test 1 completed successfully!")
  message("Results structure:")
  str(results_test1, max.level = 2)
  
}, error = function(e) {
  message("Test 1 failed: ", e$message)
})

# ============================================================================
# 테스트 2: Sender 클러스터 지정
# ============================================================================
message("\n=== Test 2: Specified Sender Clusters ===")

if (length(clusters) >= 2) {
  sender_clusters_test <- clusters[2:min(3, length(clusters))]
  
  tryCatch({
    results_test2 <- run_cci_analysis(
      sobj = sobj_test,
      cluster_col = "anno3.scvi",
      deg_df = deg_df_example,
      receiver_cluster = receiver_cluster_test,
      sender_clusters = sender_clusters_test,
      condition_col = "g3",
      condition_oi = "2",
      condition_ref = "1",
      species = "human",
    nichenet_data_dir = nichenet_data_dir_default,
      verbose = TRUE
    )
    
    message("Test 2 completed successfully!")
    
  }, error = function(e) {
    message("Test 2 failed: ", e$message)
  })
}

# ============================================================================
# 테스트 3: 자동 Sender 식별
# ============================================================================
message("\n=== Test 3: Auto-identify Sender Clusters ===")

tryCatch({
  results_test3 <- run_cci_analysis(
    sobj = sobj_test,
    cluster_col = "anno3.scvi",
    deg_df = deg_df_example,
    receiver_cluster = receiver_cluster_test,
    sender_clusters = NULL,  # 자동 식별
    condition_col = "g3",
    condition_oi = "2",
    condition_ref = "1",
    species = "human",
    nichenet_data_dir = nichenet_data_dir_default,
    verbose = TRUE
  )
  
  message("Test 3 completed successfully!")
  
}, error = function(e) {
  message("Test 3 failed: ", e$message)
})

# ============================================================================
# 결과 확인
# ============================================================================
message("\n=== Checking Results ===")

# 저장된 파일 확인
saved_files <- list.files("/data/user3/sobj", pattern = "cci.*\\.qs$", full.names = TRUE)
if (length(saved_files) > 0) {
  message("Found ", length(saved_files), " saved result file(s):")
  for (f in saved_files) {
    message("  - ", f)
  }
} else {
  message("No saved result files found")
}

# ============================================================================
# 요약
# ============================================================================
message("\n=== Test Summary ===")
message("CCI analysis tool testing completed.")
message("Check the results above for any errors.")
message("Saved results can be loaded using:")
message('  results <- qs::qread("/data/user3/sobj/cci_analysis_results_YYYYMMDD_HHMMSS.qs")')

