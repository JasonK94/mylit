# ============================================================================
# CCI 분석 도구 간단 테스트 스크립트 (기본 검증만)
# ============================================================================

# 환경 설정
# ============================================================================
message("=== CCI Tool Simple Test ===")

# 필요한 패키지 확인
required_pkgs <- c("Seurat", "dplyr", "qs")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

library(Seurat)
library(dplyr)
library(qs)

# 함수 로드
message("Loading CCI functions...")
tryCatch({
  source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
  source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
  source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
  source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
  message("  ✓ CCI functions loaded")
}, error = function(e) {
  stop("Failed to load CCI functions: ", e$message)
})

# run_nichenet_analysis 로드
message("Loading run_nichenet_analysis...")
tryCatch({
  if (file.exists("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")) {
    source("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")
    message("  ✓ run_nichenet_analysis loaded from mylit")
  } else {
    warning("CCI.R not found. run_nichenet_analysis may not be available.")
  }
}, error = function(e) {
  warning("Failed to load run_nichenet_analysis: ", e$message)
})

# 데이터 로드 테스트
message("\nTesting data loading...")
tryCatch({
  sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
  message("  ✓ Data loaded successfully")
  message("  - Cells: ", ncol(sobj_test))
  message("  - Genes: ", nrow(sobj_test))
}, error = function(e) {
  stop("Failed to load data: ", e$message)
})

# 메타데이터 확인
message("\nChecking metadata...")
colnames_meta <- colnames(sobj_test@meta.data)
message("  - Metadata columns: ", length(colnames_meta))

if ("anno3.scvi" %in% colnames_meta) {
  clusters <- unique(sobj_test@meta.data$anno3.scvi)
  clusters <- clusters[!is.na(clusters)]
  message("  - Clusters found: ", length(clusters))
  message("  - Cluster IDs: ", paste(head(clusters, 5), collapse = ", "), if(length(clusters) > 5) "..." else "")
} else {
  stop("anno3.scvi column not found in metadata")
}

if ("g3" %in% colnames_meta) {
  g3_values <- unique(sobj_test@meta.data$g3)
  g3_values <- g3_values[!is.na(g3_values)]
  message("  - g3 values: ", paste(g3_values, collapse = ", "))
} else {
  warning("g3 column not found in metadata")
}

# 함수 존재 확인
message("\nChecking function availability...")
functions_to_check <- c(
  "validate_cci_inputs",
  "extract_receiver_degs",
  "identify_sender_clusters",
  "filter_expressed_genes",
  "format_deg_summary",
  "create_sender_receiver_map",
  "save_cci_intermediate",
  "save_cci_final",
  "run_cci_analysis"
)

for (func_name in functions_to_check) {
  if (exists(func_name)) {
    message("  ✓ ", func_name)
  } else {
    message("  ✗ ", func_name, " (NOT FOUND)")
  }
}

# run_nichenet_analysis 확인
if (exists("run_nichenet_analysis")) {
  message("  ✓ run_nichenet_analysis")
} else {
  message("  ✗ run_nichenet_analysis (NOT FOUND - will cause error in full test)")
}

# 간단한 입력 검증 테스트
message("\nTesting input validation...")
receiver_cluster_test <- clusters[1]
deg_df_example <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  cluster = rep(receiver_cluster_test, 3),
  avg_log2FC = c(1.5, 2.0, -1.2),
  p_val_adj = c(0.001, 0.0001, 0.01),
  stringsAsFactors = FALSE
)

tryCatch({
  validation <- validate_cci_inputs(
    sobj = sobj_test,
    cluster_col = "anno3.scvi",
    deg_df = deg_df_example,
    receiver_cluster = receiver_cluster_test,
    sender_clusters = NULL
  )
  message("  ✓ Input validation passed")
}, error = function(e) {
  message("  ✗ Input validation failed: ", e$message)
})

# DEG 추출 테스트
message("\nTesting DEG extraction...")
tryCatch({
  receiver_degs <- extract_receiver_degs(
    deg_df_example,
    receiver_cluster_test,
    p_val_adj_cutoff = 0.05,
    logfc_cutoff = 0.25
  )
  message("  ✓ DEG extraction passed")
  message("  - Extracted DEGs: ", nrow(receiver_degs))
}, error = function(e) {
  message("  ✗ DEG extraction failed: ", e$message)
})

message("\n=== Simple Test Complete ===")
message("If all checks passed, you can proceed with full test using test_cci.R")

