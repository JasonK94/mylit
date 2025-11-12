# ============================================================================
# 대화형 테스트 스크립트
# ============================================================================
# R 세션에서 직접 실행하여 테스트할 수 있는 스크립트
# 
# 사용법:
#   1. R 세션 시작 (alias st 디렉터리에서)
#   2. source("test_interactive.R") 실행
#   3. 각 테스트를 순차적으로 실행
# ============================================================================

# 작업 디렉터리 설정
if (dir.exists("/home/user3/GJC_KDW_250721")) {
  setwd("/home/user3/GJC_KDW_250721")
  message("Working directory: ", getwd())
}

# start.R 실행 (패키지 로드)
if (file.exists("st/start.R")) {
  source("st/start.R")
} else {
  warning("st/start.R을 찾을 수 없습니다. 패키지를 수동으로 로드하세요.")
}

# 함수 소스 로드
repo_root <- "/home/user3/data_user3/git_repo/mylit-main2"
if (dir.exists(repo_root)) {
  source(file.path(repo_root, "myR/R/test_analysis.R"))
  message("Functions loaded from: ", file.path(repo_root, "myR/R/test_analysis.R"))
} else {
  stop("Repository root not found: ", repo_root)
}

# 필요한 패키지 확인
check_packages <- function(packages) {
  missing <- character(0)
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  if (length(missing) > 0) {
    warning("다음 패키지가 설치되지 않았습니다: ", paste(missing, collapse=", "))
    return(FALSE)
  }
  return(TRUE)
}

# 필수 패키지 확인
required_packages <- c("qs", "Seurat", "muscat", "SingleCellExperiment", "SummarizedExperiment", "dplyr", "Matrix")
package_check <- check_packages(required_packages)

if (!package_check) {
  message("\n일부 패키지가 누락되었습니다. 테스트를 계속할 수 없습니다.")
  message("필요한 패키지를 설치하세요.")
  stop("Required packages missing")
}

# ============================================================================
# 1. 데이터 로드
# ============================================================================
message("\n========================================")
message("1. 데이터 로드")
message("========================================")

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
}

# type 레벨 확인
type_levels <- unique(sobj@meta.data$type)
type_levels <- type_levels[!is.na(type_levels)]
message(sprintf("\ntype 레벨: %s", paste(type_levels, collapse=", ")))

# IS와 SAH가 모두 있는지 확인
if (all(c("IS", "SAH") %in% type_levels)) {
  contrast_str <- "IS - SAH"
  message(sprintf("\n사용할 contrast: %s", contrast_str))
} else {
  if (length(type_levels) >= 2) {
    contrast_groups <- type_levels[1:2]
    contrast_str <- paste(contrast_groups, collapse=" - ")
    message(sprintf("\nIS 또는 SAH가 없습니다. 대신 사용할 contrast: %s", contrast_str))
  } else {
    stop("비교할 그룹이 2개 이상 필요합니다.")
  }
}

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

# 테스트 클러스터 선택
test_clusters <- as.character(sort(clusters)[1:min(3, length(clusters))])
message(sprintf("\n테스트 클러스터: %s", paste(test_clusters, collapse=", ")))

# ============================================================================
# 2. runMUSCAT2_v1 테스트
# ============================================================================
message("\n========================================")
message("2. runMUSCAT2_v1 테스트")
message("========================================")

message("runMUSCAT2_v1 실행 중...")
message("이 작업은 몇 분에서 수십 분이 걸릴 수 있습니다.")

start_time <- Sys.time()

tryCatch({
  res_muscat2 <- runMUSCAT2_v1(
    sobj = sobj,
    cluster_id = "seurat_clusters",
    sample_id = "hos_no",
    group_id = "type",
    contrast = contrast_str,
    method = "edgeR",
    remove_na_groups = TRUE,
    keep_clusters = test_clusters
  )
  
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "secs")
  
  message(sprintf("\nrunMUSCAT2_v1 완료 (소요 시간: %.2f초)", as.numeric(elapsed_time)))
  message(sprintf("결과 행 수: %d", nrow(res_muscat2)))
  
  if (nrow(res_muscat2) > 0) {
    message("\n결과 컬럼:")
    message(paste(colnames(res_muscat2), collapse=", "))
    
    message("\n결과 미리보기:")
    print(head(res_muscat2, 5))
    
    # 결과 저장
    output_path <- "/data/user3/sobj/test_muscat2_v1_result.qs"
    qs::qsave(res_muscat2, output_path)
    message(sprintf("\n결과 저장: %s", output_path))
    message(sprintf("파일 크기: %.2f MB", file.size(output_path) / 1024 / 1024))
  } else {
    warning("결과가 비어있습니다.")
  }
  
}, error = function(e) {
  message(sprintf("\nrunMUSCAT2_v1 오류: %s", conditionMessage(e)))
  traceback()
})

# ============================================================================
# 3. runNEBULA2_v1_with_pseudobulk 테스트 (선택사항)
# ============================================================================
message("\n========================================")
message("3. runNEBULA2_v1_with_pseudobulk 테스트 (선택사항)")
message("========================================")

message("이 테스트를 실행하시겠습니까? (y/n)")
message("실행하려면 다음 코드를 실행하세요:")

cat("
res_nebula2_pb <- runNEBULA2_v1_with_pseudobulk(
  sobj = sobj,
  layer = 'counts',
  cluster_id = 'seurat_clusters',
  sample_id = 'hos_no',
  group_id = 'type',
  fixed_effects = 'g3',
  covar_effects = NULL,
  patient_col = 'hos_no',
  offset_method = 'sum',
  min_count = 10,
  min_cells_per_pb = 3,
  remove_na_cells = TRUE,
  keep_clusters = test_clusters
)

# 결과 확인
str(res_nebula2_pb, max.level = 2)

# 결과 저장
qs::qsave(res_nebula2_pb, '/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs')
")

message("\n========================================")
message("테스트 완료")
message("========================================")

