# ============================================================================
# CCI 분석 타이밍 테스트 - 병목 지점 확인
# ============================================================================

# 환경 설정
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R")

library(Seurat)
library(dplyr)
library(qs)

# 다운샘플 데이터 로드
message("Loading downsampled data...")
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
message("  Cells: ", ncol(sobj_test), ", Genes: ", nrow(sobj_test))

# 실제 발현 유전자로 DEG 리스트 준비
clusters <- unique(sobj_test@meta.data$anno3.scvi)
clusters <- clusters[!is.na(clusters)]
receiver_cluster_test <- clusters[1]

message("\n=== Preparing DEG list from actual expressed genes ===")
# Receiver cell에서 발현되는 유전자 가져오기
receiver_cells <- sobj_test@meta.data[["anno3.scvi"]] == receiver_cluster_test
receiver_expr <- Seurat::GetAssayData(sobj_test, assay = "RNA", layer = "counts")[, receiver_cells]
expressed_genes <- rownames(receiver_expr)[Matrix::rowSums(receiver_expr > 0) > 0]
# NicheNet 데이터 로드 (이미 로드되어 있을 수 있음)
if (!exists("NicheNetData")) {
  # 간단한 방법: 실제 유전자 중 일부 사용
  actual_genes <- head(expressed_genes, 20)
} else {
  # NicheNet ligand_target_matrix에 있는 유전자만 사용
  if ("ligand_target_matrix" %in% names(NicheNetData)) {
    nichenet_genes <- rownames(NicheNetData$ligand_target_matrix)
    actual_genes <- base::intersect(expressed_genes, nichenet_genes)[1:min(20, length(base::intersect(expressed_genes, nichenet_genes)))]
  } else {
    actual_genes <- head(expressed_genes, 20)
  }
}

deg_df_example <- data.frame(
  gene = actual_genes,
  cluster = rep(receiver_cluster_test, length(actual_genes)),
  avg_log2FC = c(rep(1.5, length(actual_genes) %/% 2), rep(-1.2, length(actual_genes) - length(actual_genes) %/% 2)),
  p_val_adj = rep(0.01, length(actual_genes)),
  stringsAsFactors = FALSE
)
message("  Prepared ", nrow(deg_df_example), " DEGs from expressed genes")

# 타이밍 측정
message("\n=== Starting CCI Analysis with Timing ===")
start_total <- Sys.time()

tryCatch({
  results_test <- run_cci_analysis(
    sobj = sobj_test,
    cluster_col = "anno3.scvi",
    deg_df = deg_df_example,
    receiver_cluster = receiver_cluster_test,
    condition_col = "g3",
    condition_oi = "2",
    condition_ref = "1",
    species = "human",
    nichenet_data_dir = "/data/user3/git_repo/human",
    output_dir = NULL,  # 파일 저장 안 함
    verbose = TRUE,
    top_n_targets_per_ligand = 50  # 작게 설정
  )
  
  elapsed_total <- difftime(Sys.time(), start_total, units = "secs")
  message("\n=== Total Time: ", round(elapsed_total, 1), " seconds ===")
  
}, error = function(e) {
  elapsed_total <- difftime(Sys.time(), start_total, units = "secs")
  message("\n=== Failed after ", round(elapsed_total, 1), " seconds ===")
  message("Error: ", e$message)
  traceback()
})

