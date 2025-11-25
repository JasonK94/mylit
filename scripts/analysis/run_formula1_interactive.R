# ============================================================================
# Formula 1 인터랙티브 분석 스크립트
# ============================================================================
# 이 스크립트는 R 인터랙티브 세션에서 직접 실행할 수 있습니다.
# 
# 사용법:
#   1. R 세션 시작
#      cd /home/user3/GJC_KDW_250721
#      R
#   
#   2. 이 스크립트 소스
#      source("/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_interactive.R")
#
#   3. 또는 직접 함수 호출
#      result <- run_formula1_analysis_interactive()
# ============================================================================

run_formula1_analysis_interactive <- function(
  data_path = "/data/user3/sobj/IS6_sex_added_251110.qs",
  output_dir = "/data/user3/sobj",
  min_count = 10,
  light_test = FALSE,
  max_genes = NULL,
  max_cells_per_group = NULL
) {
  
  cat("\n")
  cat("========================================\n")
  cat("Formula 1 인터랙티브 분석\n")
  cat("========================================\n\n")
  
  # ============================================================================
  # 환경 확인
  # ============================================================================
  cat("1. 환경 확인 중...\n")
  
  # 작업 디렉터리 확인
  current_wd <- getwd()
  target_wd <- "/home/user3/GJC_KDW_250721"
  if (current_wd != target_wd && dir.exists(target_wd)) {
    cat(sprintf("   작업 디렉터리 변경: %s -> %s\n", current_wd, target_wd))
    setwd(target_wd)
  }
  
  # start.R 로드 확인
  if (!exists(".myR_loaded", envir = .GlobalEnv)) {
    cat("   start.R 로드 확인 중...\n")
    start_r_path <- NULL
    if (file.exists("start.R")) {
      start_r_path <- "start.R"
    } else if (file.exists("mylit/st/start.R")) {
      start_r_path <- "mylit/st/start.R"
    } else if (file.exists("/data/user3/git_repo/mylit/st/start.R")) {
      start_r_path <- "/data/user3/git_repo/mylit/st/start.R"
    }
    
    if (!is.null(start_r_path)) {
      source(start_r_path, local = FALSE)
      assign(".myR_loaded", TRUE, envir = .GlobalEnv)
      cat(sprintf("   ✓ start.R 로드 완료: %s\n", start_r_path))
    } else {
      warning("start.R을 찾을 수 없습니다. 패키지가 로드되지 않을 수 있습니다.")
    }
  } else {
    cat("   ✓ start.R 이미 로드됨\n")
  }
  
  # 함수 소스 로드 (main2 워크트리에서)
  repo_root <- "/home/user3/data_user3/git_repo/_wt/analysis"
  if (dir.exists(repo_root)) {
    source(file.path(repo_root, "myR/R/test_analysis.R"))
    cat(sprintf("   ✓ 함수 로드 완료: %s\n", file.path(repo_root, "myR/R/test_analysis.R")))
  } else {
    stop("Repository root not found: ", repo_root)
  }
  
  # 함수 존재 확인
  if (!exists("runNEBULA")) {
    stop("runNEBULA 함수가 로드되지 않았습니다.")
  }
  cat("   ✓ runNEBULA 함수 확인됨\n\n")
  
  # ============================================================================
  # 데이터 로드
  # ============================================================================
  cat("2. 데이터 로드 중...\n")
  
  if (!file.exists(data_path)) {
    stop("데이터 파일이 없습니다: ", data_path)
  }
  
  cat(sprintf("   파일: %s\n", data_path))
  sobj <- qs::qread(data_path)
  cat(sprintf("   ✓ 데이터 로드 완료\n"))
  cat(sprintf("     - 세포 수: %d\n", ncol(sobj)))
  cat(sprintf("     - 유전자 수: %d\n", nrow(sobj)))
  
  # 경량 테스트 옵션
  if (light_test || !is.null(max_genes) || !is.null(max_cells_per_group)) {
    cat("\n   경량 테스트 모드: 데이터 서브샘플링 중...\n")
    
    original_n_cells <- ncol(sobj)
    original_n_genes <- nrow(sobj)
    
    # 세포 서브샘플링
    if (!is.null(max_cells_per_group)) {
      meta <- sobj@meta.data
      if ("g3" %in% colnames(meta)) {
        g3_levels <- unique(meta$g3)
        g3_levels <- g3_levels[!is.na(g3_levels)]
        
        cell_indices <- c()
        for (g in g3_levels) {
          g_idx <- which(meta$g3 == g & !is.na(meta$g3))
          if (length(g_idx) > max_cells_per_group) {
            g_idx <- sample(g_idx, max_cells_per_group)
          }
          cell_indices <- c(cell_indices, g_idx)
        }
        sobj <- sobj[, cell_indices]
        cat(sprintf("     - 세포 서브샘플링: %d -> %d\n", original_n_cells, ncol(sobj)))
      }
    } else if (light_test && ncol(sobj) > 2000) {
      cell_indices <- sample(1:ncol(sobj), 2000)
      sobj <- sobj[, cell_indices]
      cat(sprintf("     - 세포 랜덤 샘플링: %d -> %d\n", original_n_cells, ncol(sobj)))
    }
    
    # 유전자 필터링
    if (!is.null(max_genes)) {
      gene_counts <- Matrix::rowSums(GetAssayData(sobj, layer = "counts") > 0)
      keep_genes <- names(sort(gene_counts, decreasing = TRUE))[1:min(max_genes, length(gene_counts))]
      sobj <- sobj[keep_genes, ]
      cat(sprintf("     - 유전자 필터링: %d -> %d\n", original_n_genes, nrow(sobj)))
    } else if (light_test && nrow(sobj) > 1000) {
      gene_counts <- Matrix::rowSums(GetAssayData(sobj, layer = "counts") > 0)
      keep_genes <- names(sort(gene_counts, decreasing = TRUE))[1:1000]
      sobj <- sobj[keep_genes, ]
      cat(sprintf("     - 유전자 필터링: %d -> %d\n", original_n_genes, nrow(sobj)))
    }
  }
  
  cat("\n")
  
  # ============================================================================
  # 메타데이터 확인
  # ============================================================================
  cat("3. 메타데이터 확인 중...\n")
  
  meta <- sobj@meta.data
  required_vars <- c("g3", "sex", "anno3.scvi", "GEM", "hos_no", "nCount_RNA")
  missing_vars <- required_vars[!required_vars %in% colnames(meta)]
  if (length(missing_vars) > 0) {
    stop("필수 변수가 없습니다: ", paste(missing_vars, collapse = ", "))
  }
  cat("   ✓ 모든 필수 변수 확인됨\n")
  
  cat(sprintf("     - g3 레벨: %s\n", paste(unique(meta$g3[!is.na(meta$g3)]), collapse = ", ")))
  cat(sprintf("     - sex 레벨: %s\n", paste(unique(meta$sex[!is.na(meta$sex)]), collapse = ", ")))
  cat(sprintf("     - anno3.scvi 레벨 수: %d\n", length(unique(meta$anno3.scvi[!is.na(meta$anno3.scvi)]))))
  cat(sprintf("     - GEM 레벨 수: %d\n", length(unique(meta$GEM[!is.na(meta$GEM)]))))
  cat(sprintf("     - hos_no 샘플 수: %d\n", length(unique(meta$hos_no[!is.na(meta$hos_no)]))))
  cat("\n")
  
  # ============================================================================
  # Formula 1 분석 실행
  # ============================================================================
  cat("4. Formula 1 분석 실행\n")
  cat("========================================\n\n")
  
  # Formula 1 정의
  formula1 <- ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)
  
  cat("Formula:\n")
  cat(deparse(formula1))
  cat("\n\n")
  
  # 분석 시작 시간 기록
  start_time <- Sys.time()
  cat(sprintf("분석 시작 시간: %s\n\n", format(start_time)))
  
  # runNEBULA 실행
  cat("runNEBULA 실행 중...\n")
  if (light_test) {
    cat("(경량 테스트 모드)\n")
  }
  cat("\n")
  
  result <- tryCatch({
    runNEBULA(
      sobj = sobj,
      formula = formula1,
      patient_col = "hos_no",
      offset = "nCount_RNA",
      min_count = min_count,
      remove_na_cells = TRUE,
      layer = "counts"
    )
  }, error = function(e) {
    cat("\n✗ 분석 실패:\n")
    cat(conditionMessage(e), "\n")
    traceback()
    stop("분석 실패: ", conditionMessage(e))
  })
  
  # 분석 완료 시간 기록
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("\n분석 완료 시간: %s\n", format(end_time)))
  cat(sprintf("소요 시간: %.2f 분\n\n", as.numeric(elapsed_time)))
  
  # ============================================================================
  # 결과 확인
  # ============================================================================
  cat("========================================\n")
  cat("결과 확인\n")
  cat("========================================\n\n")
  
  cat("결과 구조:\n")
  str(result, max.level = 2)
  cat("\n")
  
  if (!is.null(result$summary)) {
    cat("결과 요약 (처음 10개):\n")
    print(head(result$summary, 10))
    cat(sprintf("\n총 유전자 수: %d\n\n", nrow(result$summary)))
  }
  
  if (!is.null(result$formula)) {
    cat("사용된 Formula:\n")
    cat(result$formula)
    cat("\n\n")
  }
  
  # ============================================================================
  # 결과 저장
  # ============================================================================
  cat("========================================\n")
  cat("결과 저장\n")
  cat("========================================\n\n")
  
  if (!dir.exists(output_dir)) {
    warning("출력 디렉터리가 없습니다: ", output_dir, "\n결과를 저장하지 않습니다.")
  } else {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    suffix <- if (light_test) "test_light" else "full"
    output_file <- file.path(output_dir, paste0("IS6_formula1_nebula_", suffix, "_", timestamp, ".qs"))
    
    cat("결과 저장 중...\n")
    cat(sprintf("  파일: %s\n", output_file))
    
    qs::qsave(result, file = output_file)
    cat("✓ 결과 저장 완료\n\n")
    
    # 메타데이터도 저장
    meta_output_file <- file.path(output_dir, paste0("IS6_formula1_metadata_", suffix, "_", timestamp, ".qs"))
    meta_info <- list(
      data_file = data_path,
      formula = deparse(formula1),
      design_formula = if (!is.null(result$design_formula)) result$design_formula else NULL,
      patient_col = if (!is.null(result$patient_col)) result$patient_col else "hos_no",
      fixed_effects = if (!is.null(result$fixed_effects)) result$fixed_effects else NULL,
      analysis_time = format(start_time),
      elapsed_time_minutes = as.numeric(elapsed_time),
      n_cells = ncol(sobj),
      n_genes = nrow(sobj),
      light_test = light_test,
      max_genes = max_genes,
      max_cells_per_group = max_cells_per_group
    )
    
    qs::qsave(meta_info, file = meta_output_file)
    cat(sprintf("✓ 메타데이터 저장 완료: %s\n\n", meta_output_file))
  }
  
  # ============================================================================
  # 요약 출력
  # ============================================================================
  cat("========================================\n")
  cat("분석 완료 요약\n")
  cat("========================================\n\n")
  cat(sprintf("입력 데이터: %s\n", data_path))
  cat(sprintf("분석 데이터: %d 세포, %d 유전자\n", ncol(sobj), nrow(sobj)))
  cat(sprintf("Formula: %s\n", deparse(formula1)))
  cat(sprintf("분석 시간: %.2f 분\n", as.numeric(elapsed_time)))
  if (dir.exists(output_dir)) {
    cat(sprintf("결과 파일: %s\n", output_file))
  }
  cat("\n✓ 모든 작업 완료!\n\n")
  
  # 결과 반환
  return(invisible(result))
}

# 함수를 바로 사용할 수 있도록 메시지 출력
cat("\n")
cat("========================================\n")
cat("Formula 1 인터랙티브 분석 함수 로드 완료\n")
cat("========================================\n\n")
cat("사용법:\n")
cat("  # 전체 데이터로 분석:\n")
cat("  result <- run_formula1_analysis_interactive()\n\n")
cat("  # 경량 테스트 (2000 세포, 1000 유전자):\n")
cat("  result <- run_formula1_analysis_interactive(light_test = TRUE)\n\n")
cat("  # 커스텀 옵션:\n")
cat("  result <- run_formula1_analysis_interactive(\n")
cat("    data_path = \"/data/user3/sobj/IS6_sex_added_251110.qs\",\n")
cat("    min_count = 10,\n")
cat("    max_genes = 500,      # 최대 500개 유전자만\n")
cat("    max_cells_per_group = 500  # 그룹당 최대 500개 세포\n")
cat("  )\n\n")
cat("결과는 result 변수에 저장되고 /data/user3/sobj/ 에 자동으로 저장됩니다.\n\n")

