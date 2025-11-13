# ============================================================================
# 테스트 스크립트: MUSCAT2_v1, NEBULA2_v1, runNEBULA2_v1_with_pseudobulk
# ============================================================================
# 이 스크립트는 /home/user3/GJC_KDW_250721 디렉터리에서 실행해야 합니다.
# 
# 사용법:
#   cd /home/user3/GJC_KDW_250721
#   Rscript "/home/user3/data_user3/git_repo/mylit-main2/test_interactive.R"
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

# start.R 실행 (패키지 로드) - 수정하지 않음
# start.R은 현재 디렉터리 또는 mylit/st/start.R에 있을 수 있음
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
# start.R이 myR을 로드했을 수도 있지만, main2 워크트리의 함수를 사용해야 함
repo_root <- "/home/user3/data_user3/git_repo/mylit-main2"
if (dir.exists(repo_root)) {
  # 함수 소스 직접 로드 (start.R이 로드한 것과 별개로)
  source(file.path(repo_root, "myR/R/test_analysis.R"))
  message("Functions loaded from: ", file.path(repo_root, "myR/R/test_analysis.R"))
} else {
  stop("Repository root not found: ", repo_root)
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
    stop(sprintf("함수 %s를 찾을 수 없습니다.", func_name))
  }
}

# ============================================================================
# 1. 데이터 로드 및 메타데이터 확인
# ============================================================================
message("\n========================================")
message("1. 데이터 로드 및 메타데이터 확인")
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
message(paste(head(meta_cols, 30), collapse=", "))

# 필수 컬럼 확인 (type은 선택사항)
required_cols <- c("seurat_clusters", "hos_no")
missing_cols <- setdiff(required_cols, meta_cols)
if (length(missing_cols) > 0) {
  stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse=", "))
}
message("✓ 필수 컬럼 확인 완료")

# 그룹 변수 확인 (type이 없으면 g3 또는 set 사용)
group_id_var <- NULL
if ("type" %in% meta_cols) {
  group_id_var <- "type"
  message("\n그룹 변수: type 사용")
} else if ("g3" %in% meta_cols) {
  group_id_var <- "g3"
  message("\n그룹 변수: g3 사용 (type 컬럼이 없음)")
} else if ("set" %in% meta_cols) {
  group_id_var <- "set"
  message("\n그룹 변수: set 사용 (type, g3 컬럼이 없음)")
} else {
  stop("그룹 변수를 찾을 수 없습니다. type, g3, 또는 set 컬럼이 필요합니다.")
}

# 그룹 레벨 확인
group_levels <- unique(sobj@meta.data[[group_id_var]])
group_levels <- group_levels[!is.na(group_levels)]
message(sprintf("%s 레벨: %s", group_id_var, paste(group_levels, collapse=", ")))

# 그룹별 세포 수
group_counts <- table(sobj@meta.data[[group_id_var]], useNA = "ifany")
message(sprintf("\n%s별 세포 수:", group_id_var))
print(group_counts)

# contrast 생성
if (length(group_levels) >= 2) {
  # 숫자형이면 문자로 변환
  if (is.numeric(group_levels)) {
    contrast_groups <- sort(group_levels)[1:2]
    contrast_str <- paste(contrast_groups, collapse=" - ")
  } else {
    contrast_groups <- group_levels[1:2]
    contrast_str <- paste(contrast_groups, collapse=" - ")
  }
  message(sprintf("\n사용할 contrast: %s", contrast_str))
} else {
  stop("비교할 그룹이 2개 이상 필요합니다.")
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
} else {
  warning("g3 컬럼이 없습니다.")
}

# 클러스터 확인
clusters <- unique(sobj@meta.data$seurat_clusters)
clusters <- clusters[!is.na(clusters)]
clusters <- sort(clusters)
message(sprintf("\n클러스터 수: %d", length(clusters)))
message(sprintf("클러스터: %s", paste(clusters[1:min(10, length(clusters))], collapse=", ")))

# 테스트 클러스터 선택 (처음 3개)
test_clusters <- as.character(clusters[1:min(3, length(clusters))])
message(sprintf("\n테스트 클러스터: %s", paste(test_clusters, collapse=", ")))

# 샘플 수 확인
if ("hos_no" %in% meta_cols) {
  n_samples <- length(unique(sobj@meta.data$hos_no))
  message(sprintf("샘플 수: %d", n_samples))
}

# ============================================================================
# 2. runMUSCAT2_v1 테스트
# ============================================================================
message("\n========================================")
message("2. runMUSCAT2_v1 테스트")
message("========================================")

message("runMUSCAT2_v1 실행 중...")
message("이 작업은 몇 분에서 수십 분이 걸릴 수 있습니다.")

start_time <- Sys.time()

test_result_muscat2 <- tryCatch({
  res_muscat2 <- runMUSCAT2_v1(
    sobj = sobj,
    cluster_id = "seurat_clusters",
    sample_id = "hos_no",
    group_id = group_id_var,  # 동적으로 결정된 그룹 변수 사용
    contrast = contrast_str,
    method = "edgeR",
    remove_na_groups = TRUE,
    keep_clusters = test_clusters
  )
  
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "secs")
  
  message(sprintf("\n✓ runMUSCAT2_v1 완료 (소요 시간: %.2f초)", as.numeric(elapsed_time)))
  message(sprintf("결과 행 수: %d", nrow(res_muscat2)))
  
  if (nrow(res_muscat2) > 0) {
    message("\n결과 컬럼:")
    message(paste(colnames(res_muscat2), collapse=", "))
    
    message("\n결과 미리보기 (상위 5개):")
    print(head(res_muscat2, 5))
    
    # 결과 저장
    output_path <- "/data/user3/sobj/test_muscat2_v1_result.qs"
    qs::qsave(res_muscat2, output_path)
    message(sprintf("\n결과 저장: %s", output_path))
    if (file.exists(output_path)) {
      message(sprintf("파일 크기: %.2f MB", file.size(output_path) / 1024 / 1024))
    }
    
    list(success = TRUE, result = res_muscat2, elapsed_time = elapsed_time)
  } else {
    warning("결과가 비어있습니다.")
    list(success = FALSE, error = "Empty result", elapsed_time = elapsed_time)
  }
  
}, error = function(e) {
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "secs")
  message(sprintf("\n✗ runMUSCAT2_v1 오류 (소요 시간: %.2f초)", as.numeric(elapsed_time)))
  message(sprintf("오류 메시지: %s", conditionMessage(e)))
  traceback()
  list(success = FALSE, error = conditionMessage(e), elapsed_time = elapsed_time)
})

# ============================================================================
# 3. runNEBULA2_v1_with_pseudobulk 테스트
# ============================================================================
message("\n========================================")
message("3. runNEBULA2_v1_with_pseudobulk 테스트")
message("========================================")

message("runNEBULA2_v1_with_pseudobulk 실행 중...")
message("이 작업은 몇 분에서 수십 분이 걸릴 수 있습니다.")

start_time_pb <- Sys.time()

test_result_nebula2_pb <- tryCatch({
  # g3가 그룹 변수로 사용된 경우, fixed_effects에는 다른 변수 사용
  # 그렇지 않으면 g3를 fixed_effects로 사용
  if (group_id_var == "g3") {
    # g3가 그룹 변수로 사용되면, fixed_effects에는 다른 변수 사용 (예: set)
    if ("set" %in% meta_cols) {
      fixed_effects_var <- "set"
      message("g3가 그룹 변수로 사용되므로, fixed_effects에는 'set' 사용")
    } else {
      # set도 없으면 g3를 그대로 사용 (group_id와 동일)
      fixed_effects_var <- "g3"
      message("g3가 그룹 변수로 사용되며, 다른 변수가 없으므로 fixed_effects에도 'g3' 사용")
    }
  } else {
    # g3가 그룹 변수가 아니면 fixed_effects로 사용
    fixed_effects_var <- "g3"
    message(sprintf("그룹 변수는 '%s'이므로, fixed_effects에는 'g3' 사용", group_id_var))
  }
  
  # g3에 결측치가 있는지 확인
  if (fixed_effects_var == "g3" && "g3" %in% meta_cols) {
    n_na_g3 <- sum(is.na(sobj@meta.data$g3))
    if (n_na_g3 > 0) {
      message(sprintf("경고: g3에 %d개의 결측치가 있습니다. remove_na_cells=TRUE로 제거됩니다.", n_na_g3))
    }
  }
  
  # 유전자 수가 많으면 NEBULA가 실패할 수 있으므로, 더 강한 필터링 적용
  # 또는 처음 1000개 유전자만 테스트
  # min_count를 높이거나, 유전자를 제한할 수 있음
  res_nebula2_pb <- runNEBULA2_v1_with_pseudobulk(
    sobj = sobj,
    layer = "counts",
    cluster_id = "seurat_clusters",
    sample_id = "hos_no",
    group_id = group_id_var,  # 동적으로 결정된 그룹 변수 사용
    fixed_effects = fixed_effects_var,
    covar_effects = NULL,
    patient_col = "hos_no",
    offset_method = "sum",
    min_count = 20,  # 더 높은 min_count로 필터링 강화
    min_cells_per_pb = 3,
    remove_na_cells = TRUE,
    keep_clusters = test_clusters
  )
  
  end_time_pb <- Sys.time()
  elapsed_time_pb <- difftime(end_time_pb, start_time_pb, units = "secs")
  
  message(sprintf("\n✓ runNEBULA2_v1_with_pseudobulk 완료 (소요 시간: %.2f초)", as.numeric(elapsed_time_pb)))
  
  # 결과 확인
  message("\n결과 구조:")
  message(sprintf("  - nebula_result: %s", class(res_nebula2_pb$nebula_result)[1]))
  message(sprintf("  - pseudobulk_meta: %d 행", nrow(res_nebula2_pb$pseudobulk_meta)))
  message(sprintf("  - pseudobulk_counts: %d 유전자 x %d 샘플", 
                  nrow(res_nebula2_pb$pseudobulk_counts), 
                  ncol(res_nebula2_pb$pseudobulk_counts)))
  
  # 결과 저장
  output_path_pb <- "/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs"
  qs::qsave(res_nebula2_pb, output_path_pb)
  message(sprintf("\n결과 저장: %s", output_path_pb))
  if (file.exists(output_path_pb)) {
    message(sprintf("파일 크기: %.2f MB", file.size(output_path_pb) / 1024 / 1024))
  }
  
  list(success = TRUE, result = res_nebula2_pb, elapsed_time = elapsed_time_pb)
  
}, error = function(e) {
  end_time_pb <- Sys.time()
  elapsed_time_pb <- difftime(end_time_pb, start_time_pb, units = "secs")
  message(sprintf("\n✗ runNEBULA2_v1_with_pseudobulk 오류 (소요 시간: %.2f초)", as.numeric(elapsed_time_pb)))
  message(sprintf("오류 메시지: %s", conditionMessage(e)))
  traceback()
  list(success = FALSE, error = conditionMessage(e), elapsed_time = elapsed_time_pb)
})

# ============================================================================
# 4. 테스트 결과 요약
# ============================================================================
message("\n========================================")
message("4. 테스트 결과 요약")
message("========================================")

# runMUSCAT2_v1 결과
if (test_result_muscat2$success) {
  message("✓ runMUSCAT2_v1: 성공")
  message(sprintf("  - 소요 시간: %.2f초", as.numeric(test_result_muscat2$elapsed_time)))
  message(sprintf("  - 결과 행 수: %d", nrow(test_result_muscat2$result)))
} else {
  message("✗ runMUSCAT2_v1: 실패")
  message(sprintf("  - 오류: %s", test_result_muscat2$error))
}

# runNEBULA2_v1_with_pseudobulk 결과
if (test_result_nebula2_pb$success) {
  message("✓ runNEBULA2_v1_with_pseudobulk: 성공")
  message(sprintf("  - 소요 시간: %.2f초", as.numeric(test_result_nebula2_pb$elapsed_time)))
  message(sprintf("  - Pseudobulk 샘플 수: %d", nrow(test_result_nebula2_pb$result$pseudobulk_meta)))
  message(sprintf("  - Pseudobulk 유전자 수: %d", nrow(test_result_nebula2_pb$result$pseudobulk_counts)))
} else {
  message("✗ runNEBULA2_v1_with_pseudobulk: 실패")
  message(sprintf("  - 오류: %s", test_result_nebula2_pb$error))
}

message("\n========================================")
message("테스트 완료")
message("========================================")

# 모든 테스트가 성공했는지 확인
all_success <- test_result_muscat2$success && test_result_nebula2_pb$success
if (all_success) {
  message("\n✓ 모든 테스트가 성공적으로 완료되었습니다!")
  quit(status = 0)
} else {
  message("\n⚠ 일부 테스트가 실패했습니다. 오류 메시지를 확인하세요.")
  quit(status = 1)
}
