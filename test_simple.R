# ============================================================================
# 간단한 데이터 확인 및 함수 테스트
# ============================================================================

# 작업 디렉터리 설정
if (dir.exists("/home/user3/GJC_KDW_250721")) {
  setwd("/home/user3/GJC_KDW_250721")
}

# start.R 실행 (패키지 로드)
if (file.exists("st/start.R")) {
  source("st/start.R")
} else if (file.exists("../st/start.R")) {
  source("../st/start.R")
} else if (file.exists("/home/user3/data_user3/git_repo/mylit-main2/st/start.R")) {
  source("/home/user3/data_user3/git_repo/mylit-main2/st/start.R")
}

# 함수 소스 로드
repo_root <- "/home/user3/data_user3/git_repo/mylit-main2"
if (dir.exists(repo_root)) {
  source(file.path(repo_root, "myR/R/test_analysis.R"))
}

# qs 패키지 확인
if (!requireNamespace("qs", quietly = TRUE)) {
  stop("qs 패키지가 필요합니다.")
}

message("========================================")
message("데이터 확인 및 함수 테스트")
message("========================================")

# 데이터 로드
test_data_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
if (!file.exists(test_data_path)) {
  stop("테스트 데이터 파일이 없습니다: ", test_data_path)
}

message("데이터 로드 중...")
sobj <- qs::qread(test_data_path)
message(sprintf("데이터 로드 완료: %d 세포, %d 유전자", ncol(sobj), nrow(sobj)))

# 메타데이터 확인
message("\n메타데이터 컬럼:")
meta_cols <- colnames(sobj@meta.data)
message(paste(head(meta_cols, 20), collapse=", "))

# 필수 컬럼 확인
required_cols <- c("seurat_clusters", "hos_no", "type")
missing_cols <- setdiff(required_cols, meta_cols)
if (length(missing_cols) > 0) {
  stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse=", "))
} else {
  message("\n필수 컬럼 확인 완료")
}

# type 레벨 확인
type_levels <- unique(sobj@meta.data$type)
type_levels <- type_levels[!is.na(type_levels)]
message(sprintf("\ntype 레벨: %s", paste(type_levels, collapse=", ")))

# type별 세포 수
type_counts <- table(sobj@meta.data$type, useNA = "ifany")
message("\ntype별 세포 수:")
print(type_counts)

# g3 확인
if ("g3" %in% meta_cols) {
  n_na_g3 <- sum(is.na(sobj@meta.data$g3))
  n_total <- nrow(sobj@meta.data)
  pct_na <- 100 * n_na_g3 / n_total
  message(sprintf("\ng3 결측치: %d / %d (%.2f%%)", n_na_g3, n_total, pct_na))
  
  g3_levels <- unique(sobj@meta.data$g3)
  g3_levels <- g3_levels[!is.na(g3_levels)]
  message(sprintf("g3 레벨: %s", paste(g3_levels, collapse=", ")))
}

# 클러스터 확인
clusters <- unique(sobj@meta.data$seurat_clusters)
clusters <- clusters[!is.na(clusters)]
message(sprintf("\n클러스터 수: %d", length(clusters)))
message(sprintf("클러스터: %s", paste(sort(clusters)[1:min(10, length(clusters))], collapse=", ")))

# 샘플 수 확인
if ("hos_no" %in% meta_cols) {
  n_samples <- length(unique(sobj@meta.data$hos_no))
  message(sprintf("\n샘플 수: %d", n_samples))
}

# 함수 존재 확인
message("\n========================================")
message("함수 존재 확인")
message("========================================")

functions_to_check <- c("runMUSCAT2_v1", "runNEBULA2_v1", "runNEBULA2_v1_with_pseudobulk")

for (func_name in functions_to_check) {
  if (exists(func_name)) {
    message(sprintf("✓ %s: 존재함", func_name))
  } else {
    message(sprintf("✗ %s: 존재하지 않음", func_name))
  }
}

# 패키지 확인
message("\n========================================")
message("패키지 확인")
message("========================================")

packages_to_check <- c("Seurat", "muscat", "SingleCellExperiment", "SummarizedExperiment", "nebula", "dplyr", "Matrix")

for (pkg in packages_to_check) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("✓ %s: 설치됨", pkg))
  } else {
    message(sprintf("✗ %s: 설치되지 않음", pkg))
  }
}

message("\n========================================")
message("데이터 확인 완료")
message("========================================")
message("\n다음 단계:")
message("1. 필요한 패키지가 모두 설치되어 있는지 확인")
message("2. runMUSCAT2_v1 함수 테스트 실행")
message("3. 결과 확인")

