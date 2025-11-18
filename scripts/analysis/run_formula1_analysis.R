# ============================================================================
# Formula 1 분석 스크립트
# Formula: ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/patient)
# ============================================================================
# 이 스크립트는 /home/user3/GJC_KDW_250721 디렉터리에서 실행해야 합니다.
# 
# 사용법:
#   cd /home/user3/GJC_KDW_250721
#   Rscript "/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_analysis.R"
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

# 함수 존재 확인
if (!exists("runNEBULA")) {
  stop("runNEBULA 함수가 로드되지 않았습니다.")
}
message("✓ runNEBULA 함수 확인됨")

# ============================================================================
# 데이터 로드
# ============================================================================
message("\n========================================")
message("데이터 로드 중...")
message("========================================\n")

data_path <- "/data/user3/sobj/IS6_sex_added_251110.qs"
if (!file.exists(data_path)) {
  stop("데이터 파일이 없습니다: ", data_path)
}

message("데이터 파일: ", data_path)
is6 <- qs::qread(data_path)
message("✓ 데이터 로드 완료")
message(sprintf("  - 세포 수: %d", ncol(is6)))
message(sprintf("  - 유전자 수: %d", nrow(is6)))

# ============================================================================
# 메타데이터 확인
# ============================================================================
message("\n========================================")
message("메타데이터 확인 중...")
message("========================================\n")

meta <- is6@meta.data
message("메타데이터 컬럼:")
print(colnames(meta))

# 필수 변수 확인
required_vars <- c("g3", "sex", "anno3.scvi", "GEM", "hos_no", "nCount_RNA")
missing_vars <- required_vars[!required_vars %in% colnames(meta)]
if (length(missing_vars) > 0) {
  stop("필수 변수가 없습니다: ", paste(missing_vars, collapse = ", "))
}
message("✓ 모든 필수 변수 확인됨")

# 각 변수 확인
message("\n변수별 레벨 확인:")
message("\n1. g3 (그룹):")
print(table(meta$g3, useNA = "ifany"))

message("\n2. sex (성별):")
print(table(meta$sex, useNA = "ifany"))

message("\n3. anno3.scvi (셀 타입):")
print(table(meta$anno3.scvi, useNA = "ifany"))

message("\n4. GEM:")
print(table(meta$GEM, useNA = "ifany"))

message("\n5. hos_no (샘플 ID):")
n_samples <- length(unique(meta$hos_no))
message(sprintf("  - 고유 샘플 수: %d", n_samples))

# GEM과 hos_no 관계 확인 (nested structure 확인)
message("\n6. GEM x hos_no 관계 (nested structure):")
gem_hos_table <- table(meta$GEM, meta$hos_no)
message(sprintf("  - GEM 레벨 수: %d", nrow(gem_hos_table)))
message(sprintf("  - hos_no 레벨 수: %d", ncol(gem_hos_table)))

# ============================================================================
# Formula 1 분석 실행
# ============================================================================
message("\n========================================")
message("Formula 1 분석 실행")
message("========================================\n")

# Formula 1 정의
# ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/patient)
# Note: (1|GEM/patient)는 함수 내에서 (1|patient)로 변환되고 GEM이 fixed effect로 추가됨

formula1 <- ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)

message("Formula:")
message(deparse(formula1))

# 분석 시작 시간 기록
start_time <- Sys.time()
message(sprintf("\n분석 시작 시간: %s", format(start_time)))

# runNEBULA 실행
message("\nrunNEBULA 실행 중...")
message("(이 작업은 시간이 오래 걸릴 수 있습니다.)")
message(" - complete separation guard: 레벨당 최소 25개 세포")
message(" - 그룹 균형 필터: g3 기준, 레벨당 최소 5 셀")
message(" - interaction simplification: 활성화\n")

result <- tryCatch({
  runNEBULA(
    sobj = is6,
    formula = formula1,
    patient_col = "hos_no",  # patient 컬럼명 명시
    offset = "nCount_RNA",
    min_count = 10,  # 최소 10개 세포에서 발현된 유전자만
    remove_na_cells = TRUE,  # NA 값이 있는 세포 제거
    layer = "counts",
    separation_min_cells = 25,
    simplify_interactions = TRUE,
    separation_group_var = "g3",
    separation_min_cells_per_group = 5
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
  message("\n결과 요약:")
  print(head(result$summary, 10))
  message(sprintf("\n총 유전자 수: %d", nrow(result$summary)))
}

if (!is.null(result$formula)) {
  message("\n사용된 Formula:")
  message(result$formula)
}

if (!is.null(result$design_formula)) {
  message("\nDesign Formula:")
  message(result$design_formula)
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

# 결과 파일명
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- file.path(output_dir, paste0("IS6_formula1_nebula_result_", timestamp, ".qs"))

message("결과 저장 중...")
message(sprintf("  파일: %s", output_file))

qs::qsave(result, file = output_file)
message("✓ 결과 저장 완료")

# 메타데이터도 함께 저장 (선택사항)
meta_output_file <- file.path(output_dir, paste0("IS6_formula1_metadata_", timestamp, ".qs"))
meta_info <- list(
  data_file = data_path,
  formula = deparse(formula1),
  design_formula = if (!is.null(result$design_formula)) result$design_formula else NULL,
  patient_col = if (!is.null(result$patient_col)) result$patient_col else "hos_no",
  fixed_effects = if (!is.null(result$fixed_effects)) result$fixed_effects else NULL,
  analysis_time = format(start_time),
  elapsed_time_minutes = as.numeric(elapsed_time),
  n_cells = ncol(is6),
  n_genes = nrow(is6),
  n_samples = n_samples
)

qs::qsave(meta_info, file = meta_output_file)
message(sprintf("✓ 메타데이터 저장 완료: %s", meta_output_file))

# ============================================================================
# 요약 출력
# ============================================================================
message("\n========================================")
message("분석 완료 요약")
message("========================================\n")
message(sprintf("입력 데이터: %s", data_path))
message(sprintf("Formula: %s", deparse(formula1)))
message(sprintf("분석 시간: %.2f 분", as.numeric(elapsed_time)))
message(sprintf("결과 파일: %s", output_file))
message(sprintf("메타데이터 파일: %s", meta_output_file))
message("\n✓ 모든 작업 완료!")

