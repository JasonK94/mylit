# ============================================================================
# Formula 1 경량 테스트 스크립트
# Formula: ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/patient)
# ============================================================================
# 이 스크립트는 /home/user3/GJC_KDW_250721 디렉터리에서 실행해야 합니다.
# 
# 사용법:
#   cd /home/user3/GJC_KDW_250721
#   Rscript "/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_formula1_light.R"
# ============================================================================

# 작업 디렉터리를 /home/user3/GJC_KDW_250721로 설정 (중요!)
original_wd <- getwd()
target_wd <- "/home/user3/GJC_KDW_250721"
if (dir.exists(target_wd)) {
  setwd(target_wd)
  message("Working directory set to: ", getwd())
} else {
  stop("Required working directory not found: ", target_wd)
}

# start.R 실행 (패키지 로드)
start_r_path <- NULL
if (file.exists("start.R")) {
  start_r_path <- "start.R"
} else if (file.exists("mylit/st/start.R")) {
  start_r_path <- "mylit/st/start.R"
} else if (file.exists("/data/user3/git_repo/mylit/st/start.R")) {
  start_r_path <- "/data/user3/git_repo/mylit/st/start.R"
}

if (!is.null(start_r_path)) {
  message("Loading start.R from: ", start_r_path)
  source(start_r_path)
} else {
  stop("start.R not found. Current directory: ", getwd())
}

# 함수 소스 로드 (main2 워크트리에서)
repo_root <- "/home/user3/data_user3/git_repo/_wt/analysis"
if (dir.exists(repo_root)) {
  pkg_root <- file.path(repo_root, "myR")
  if ("package:myR" %in% search()) {
    detach("package:myR", unload = TRUE, character.only = TRUE)
  }
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    source(file.path(pkg_root, "R", "analysis.R"))
  }
  message("Functions loaded from: ", pkg_root)
} else {
  stop("Repository root not found: ", repo_root)
}

# 필요한 패키지 확인
required_packages <- c("Seurat", "nebula", "qs", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("필요한 패키지가 없습니다: ", paste(missing_packages, collapse = ", "))
}

if (getOption("expressions") < 50000) {
  options(expressions = 50000)
}

options(warn = 1)

# 함수 존재 확인
if (!exists("runNEBULA")) {
  stop("runNEBULA 함수가 로드되지 않았습니다.")
}
message("✓ runNEBULA 함수 확인됨")

# ============================================================================
# 데이터 로드 및 서브샘플링 (경량 테스트용)
# ============================================================================
message("\n========================================")
message("데이터 로드 및 서브샘플링 (경량 테스트)")
message("========================================\n")

data_path <- "/data/user3/sobj/IS6_sex_added_251110.qs"
if (!file.exists(data_path)) {
  stop("데이터 파일이 없습니다: ", data_path)
}

message("데이터 파일: ", data_path)
is6 <- qs::qread(data_path)
message("✓ 원본 데이터 로드 완료")
message(sprintf("  - 원본 세포 수: %d", ncol(is6)))
message(sprintf("  - 원본 유전자 수: %d", nrow(is6)))

# 경량 테스트를 위한 서브샘플링
message("\n경량 테스트를 위한 서브샘플링 중...")

# 1. 세포 서브샘플링: 각 그룹에서 최대 1000개씩
meta <- is6@meta.data
if ("g3" %in% colnames(meta)) {
  g3_levels <- unique(meta$g3)
  g3_levels <- g3_levels[!is.na(g3_levels)]
  
  cell_indices <- c()
  max_cells_per_group <- 300
  
  for (g in g3_levels) {
    g_idx <- which(meta$g3 == g & !is.na(meta$g3))
    if (length(g_idx) > max_cells_per_group) {
      g_idx <- sample(g_idx, max_cells_per_group)
    }
    cell_indices <- c(cell_indices, g_idx)
  }
  
  is6_sub <- is6[, cell_indices]
  message(sprintf("  - 세포 서브샘플링: %d -> %d", ncol(is6), ncol(is6_sub)))
} else {
  # g3가 없으면 랜덤 샘플링
  max_cells <- 600
  if (ncol(is6) > max_cells) {
    cell_indices <- sample(1:ncol(is6), max_cells)
    is6_sub <- is6[, cell_indices]
    message(sprintf("  - 세포 랜덤 샘플링: %d -> %d", ncol(is6), ncol(is6_sub)))
  } else {
    is6_sub <- is6
  }
}

# 2. 유전자 필터링: 많이 발현된 유전자만 선택 (처음 1000개)
message("  - 유전자 필터링 중...")
gene_counts <- Matrix::rowSums(GetAssayData(is6_sub, layer = "counts") > 0)
keep_genes <- names(sort(gene_counts, decreasing = TRUE))[1:min(500, length(gene_counts))]
is6_sub <- is6_sub[keep_genes, ]
message(sprintf("  - 유전자 필터링: %d -> %d", nrow(is6), nrow(is6_sub)))

message("✓ 서브샘플링 완료")
message(sprintf("  - 최종 테스트 데이터: %d 세포, %d 유전자", ncol(is6_sub), nrow(is6_sub)))

# ============================================================================
# 메타데이터 확인
# ============================================================================
message("\n========================================")
message("메타데이터 확인")
message("========================================\n")

meta_sub <- is6_sub@meta.data
required_vars <- c("g3", "sex", "anno3.scvi", "GEM", "hos_no", "nCount_RNA")
missing_vars <- required_vars[!required_vars %in% colnames(meta_sub)]
if (length(missing_vars) > 0) {
  stop("필수 변수가 없습니다: ", paste(missing_vars, collapse = ", "))
}
message("✓ 모든 필수 변수 확인됨")

# ============================================================================
# Formula 1 경량 테스트 실행
# ============================================================================
message("\n========================================")
message("Formula 1 경량 테스트 실행")
message("========================================\n")

# Formula 1 정의
formula1 <- ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)

message("Formula:")
message(deparse(formula1))

# 분석 시작 시간 기록
start_time <- Sys.time()
message(sprintf("\n분석 시작 시간: %s", format(start_time)))

# runNEBULA 실행 (formula 모드)
message("\nrunNEBULA 실행 중...")
message("(경량 테스트: 시간이 적게 걸립니다.)")
message(" - complete separation guard: 레벨당 최소 15개 세포")
message(" - 그룹 균형 필터: g3 기준, 레벨당 최소 3 셀")
message(" - interaction simplification: 활성화\n")

result <- tryCatch({
  runNEBULA(
    sobj = is6_sub,
    formula = formula1,
    patient_col = "hos_no",
    offset = "nCount_RNA",
    min_count = 5,  # 경량 테스트: 더 낮은 threshold
    remove_na_cells = TRUE,
    layer = "counts",
    separation_min_cells = 15,
    simplify_interactions = TRUE,
    separation_group_var = "g3",
    separation_min_cells_per_group = 3
  )
}, error = function(e) {
  message("\n✗ 분석 실패:")
  message(conditionMessage(e))
  traceback()
  stop("분석 실패: ", conditionMessage(e))
})

# 분석 완료 시간 기록
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")
message(sprintf("\n분석 완료 시간: %s", format(end_time)))
message(sprintf("소요 시간: %.2f 분", as.numeric(elapsed_time)))

# ============================================================================
# 결과 확인
# ============================================================================
message("\n========================================")
message("결과 확인")
message("========================================\n")

message("결과 구조:")
str(result, max.level = 2)

if (!is.null(result$summary)) {
  message("\n결과 요약 (처음 10개):")
  print(head(result$summary, 10))
  message(sprintf("\n총 유전자 수: %d", nrow(result$summary)))
}

if (!is.null(result$formula)) {
  message("\n사용된 Formula:")
  message(result$formula)
}

# ============================================================================
# 결과 저장
# ============================================================================
message("\n========================================")
message("결과 저장")
message("========================================\n")

output_dir <- "/data/user3/sobj"
if (!dir.exists(output_dir)) {
  stop("출력 디렉터리가 없습니다: ", output_dir)
}

# 결과 파일명 (경량 테스트 표시)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- file.path(output_dir, paste0("IS6_formula1_nebula_test_light_", timestamp, ".qs"))

message("결과 저장 중...")
message(sprintf("  파일: %s", output_file))

qs::qsave(result, file = output_file)
message("✓ 결과 저장 완료")

# ============================================================================
# 요약 출력
# ============================================================================
message("\n========================================")
message("경량 테스트 완료 요약")
message("========================================\n")
message(sprintf("입력 데이터: %s", data_path))
message(sprintf("테스트 데이터: %d 세포, %d 유전자", ncol(is6_sub), nrow(is6_sub)))
message(sprintf("Formula: %s", deparse(formula1)))
message(sprintf("분석 시간: %.2f 분", as.numeric(elapsed_time)))
message(sprintf("결과 파일: %s", output_file))
message("\n✓ 경량 테스트 완료!")

