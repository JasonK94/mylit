# ============================================================================
# Step 1: 데이터 로드 및 메타데이터 확인
# ============================================================================

# 함수 소스 로드 (패키지 대신 직접 소스)
source("myR/R/test_analysis.R")

# 필요한 패키지 로드
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs 패키지가 필요합니다.")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat 패키지가 필요합니다.")
}

message("========================================")
message("Step 1: 데이터 로드 및 메타데이터 확인")
message("========================================")

# 데이터 로드
test_data_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
if (!file.exists(test_data_path)) {
  stop(sprintf("테스트 데이터 파일이 없습니다: %s", test_data_path))
}

message("데이터 로드 중...")
sobj <- qs::qread(test_data_path)

message(sprintf("데이터 로드 완료: %d 세포, %d 유전자", 
                ncol(sobj), nrow(sobj)))

# 메타데이터 확인
message("\n메타데이터 컬럼:")
meta_cols <- colnames(sobj@meta.data)
message(paste(meta_cols, collapse = ", "))

# 필수 컬럼 확인
required_cols <- c("seurat_clusters", "hos_no", "type")
missing_cols <- setdiff(required_cols, meta_cols)
if (length(missing_cols) > 0) {
  warning(sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse=", ")))
} else {
  message("\n필수 컬럼 확인 완료")
}

# g3 결측치 확인
if ("g3" %in% meta_cols) {
  n_na_g3 <- sum(is.na(sobj@meta.data$g3))
  n_total <- nrow(sobj@meta.data)
  pct_na <- 100 * n_na_g3 / n_total
  message(sprintf("\ng3 결측치: %d / %d (%.2f%%)", 
                  n_na_g3, n_total, pct_na))
  
  # g3 레벨 확인
  g3_levels <- unique(sobj@meta.data$g3)
  g3_levels <- g3_levels[!is.na(g3_levels)]
  message(sprintf("g3 레벨: %s", paste(g3_levels, collapse=", ")))
} else {
  warning("g3 컬럼이 없습니다.")
}

# type 레벨 확인
if ("type" %in% meta_cols) {
  type_levels <- unique(sobj@meta.data$type)
  type_levels <- type_levels[!is.na(type_levels)]
  message(sprintf("\ntype 레벨: %s", paste(type_levels, collapse=", ")))
  
  # 각 type별 세포 수
  type_counts <- table(sobj@meta.data$type, useNA = "ifany")
  message("\ntype별 세포 수:")
  print(type_counts)
} else {
  warning("type 컬럼이 없습니다.")
}

# seurat_clusters 확인
if ("seurat_clusters" %in% meta_cols) {
  cluster_levels <- unique(sobj@meta.data$seurat_clusters)
  cluster_levels <- cluster_levels[!is.na(cluster_levels)]
  message(sprintf("\nseurat_clusters 레벨: %s", paste(sort(cluster_levels), collapse=", ")))
  message(sprintf("클러스터 수: %d", length(cluster_levels)))
} else {
  warning("seurat_clusters 컬럼이 없습니다.")
}

# hos_no 확인
if ("hos_no" %in% meta_cols) {
  n_samples <- length(unique(sobj@meta.data$hos_no))
  message(sprintf("\nhos_no 샘플 수: %d", n_samples))
} else {
  warning("hos_no 컬럼이 없습니다.")
}

# counts layer 확인
if ("counts" %in% Seurat::Layers(sobj)) {
  message("\ncounts layer 존재함")
} else {
  warning("counts layer가 없을 수 있습니다.")
}

message("\n========================================")
message("Step 1 완료")
message("========================================")

# 데이터 저장 (다음 단계에서 사용)
save(sobj, file = "/tmp/sobj_test.RData")
message("\n데이터를 /tmp/sobj_test.RData에 저장했습니다.")

