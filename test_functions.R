# ============================================================================
# 테스트 스크립트: MUSCAT2_v1, NEBULA2_v1, runNEBULA2_v1_with_pseudobulk
# ============================================================================
# 이 스크립트는 새로 개발된 함수들을 테스트합니다.
# 
# 사용법:
#   1. R 세션에서 이 스크립트를 실행
#   2. 테스트 데이터 경로 확인: /data/user3/sobj/IS_scvi_251107_ds2500.qs
#   3. 각 함수를 순차적으로 테스트
# ============================================================================

library(myR)
library(qs)

# ============================================================================
# 1. 데이터 로드
# ============================================================================
message("========================================")
message("테스트 데이터 로드 중...")
message("========================================")

test_data_path <- "/data/user3/sobj/IS_scvi_251107_ds2500.qs"
if (!file.exists(test_data_path)) {
  stop(sprintf("테스트 데이터 파일이 없습니다: %s", test_data_path))
}

sobj <- qs::qread(test_data_path)
message(sprintf("데이터 로드 완료: %d 세포, %d 유전자", 
                ncol(sobj), nrow(sobj)))

# 메타데이터 확인
message("\n메타데이터 컬럼:")
message(paste(colnames(sobj@meta.data), collapse = ", "))

# g3 결측치 확인
if ("g3" %in% colnames(sobj@meta.data)) {
  n_na_g3 <- sum(is.na(sobj@meta.data$g3))
  message(sprintf("\ng3 결측치: %d / %d (%.2f%%)", 
                  n_na_g3, nrow(sobj@meta.data), 
                  100 * n_na_g3 / nrow(sobj@meta.data)))
}

# ============================================================================
# 2. runMUSCAT2_v1 테스트
# ============================================================================
message("\n========================================")
message("runMUSCAT2_v1 테스트 시작...")
message("========================================")

# 테스트 실행 (작은 클러스터만)
test_muscat2 <- FALSE  # TRUE로 변경하여 테스트 실행
if (test_muscat2) {
  tryCatch({
    res_muscat2 <- runMUSCAT2_v1(
      sobj = sobj,
      cluster_id = "seurat_clusters",
      sample_id = "hos_no",
      group_id = "type",
      contrast = "IS - SAH",
      method = "edgeR",
      remove_na_groups = TRUE,
      keep_clusters = c("0", "1", "2")  # 처음 3개 클러스터만 테스트
    )
    
    message(sprintf("\nrunMUSCAT2_v1 완료: %d 행", nrow(res_muscat2)))
    message("결과 컬럼:")
    message(paste(colnames(res_muscat2), collapse = ", "))
    
    # 결과 저장
    output_path <- "/data/user3/sobj/test_muscat2_v1_result.qs"
    qs::qsave(res_muscat2, output_path)
    message(sprintf("\n결과 저장: %s", output_path))
    
  }, error = function(e) {
    message(sprintf("\nrunMUSCAT2_v1 오류: %s", conditionMessage(e)))
  })
} else {
  message("runMUSCAT2_v1 테스트 건너뜀 (test_muscat2 = FALSE)")
}

# ============================================================================
# 3. runNEBULA2_v1 테스트
# ============================================================================
message("\n========================================")
message("runNEBULA2_v1 테스트 시작...")
message("========================================")

# 테스트 실행
test_nebula2 <- FALSE  # TRUE로 변경하여 테스트 실행
if (test_nebula2) {
  tryCatch({
    # g3가 있는지 확인
    if (!"g3" %in% colnames(sobj@meta.data)) {
      message("g3 컬럼이 없습니다. 다른 변수를 사용하세요.")
    } else {
      # 작은 서브셋으로 테스트 (처음 1000개 유전자만)
      sobj_sub <- sobj[1:min(1000, nrow(sobj)), ]
      
      res_nebula2 <- runNEBULA2_v1(
        sobj = sobj_sub,
        layer = "counts",
        fixed_effects = "g3",
        covar_effects = NULL,
        patient_col = "hos_no",
        offset = "nCount_RNA",
        min_count = 10,
        remove_na_cells = TRUE
      )
      
      message(sprintf("\nrunNEBULA2_v1 완료"))
      message("결과 구조:")
      str(res_nebula2, max.level = 2)
      
      # 결과 저장
      output_path <- "/data/user3/sobj/test_nebula2_v1_result.qs"
      qs::qsave(res_nebula2, output_path)
      message(sprintf("\n결과 저장: %s", output_path))
    }
    
  }, error = function(e) {
    message(sprintf("\nrunNEBULA2_v1 오류: %s", conditionMessage(e)))
  })
} else {
  message("runNEBULA2_v1 테스트 건너뜀 (test_nebula2 = FALSE)")
}

# ============================================================================
# 4. runNEBULA2_v1_with_pseudobulk 테스트
# ============================================================================
message("\n========================================")
message("runNEBULA2_v1_with_pseudobulk 테스트 시작...")
message("========================================")

# 테스트 실행
test_nebula2_pb <- FALSE  # TRUE로 변경하여 테스트 실행
if (test_nebula2_pb) {
  tryCatch({
    # g3가 있는지 확인
    if (!"g3" %in% colnames(sobj@meta.data)) {
      message("g3 컬럼이 없습니다. 다른 변수를 사용하세요.")
    } else {
      # 작은 서브셋으로 테스트 (처음 몇 개 클러스터만)
      res_nebula2_pb <- runNEBULA2_v1_with_pseudobulk(
        sobj = sobj,
        layer = "counts",
        cluster_id = "seurat_clusters",
        sample_id = "hos_no",
        group_id = "type",
        fixed_effects = "g3",
        covar_effects = NULL,
        patient_col = "hos_no",
        offset_method = "sum",
        min_count = 10,
        min_cells_per_pb = 3,
        remove_na_cells = TRUE,
        keep_clusters = c("0", "1", "2")  # 처음 3개 클러스터만 테스트
      )
      
      message(sprintf("\nrunNEBULA2_v1_with_pseudobulk 완료"))
      message("결과 구조:")
      str(res_nebula2_pb, max.level = 2)
      
      message(sprintf("\nPseudobulk 샘플 수: %d", nrow(res_nebula2_pb$pseudobulk_meta)))
      message(sprintf("Pseudobulk 유전자 수: %d", nrow(res_nebula2_pb$pseudobulk_counts)))
      
      # 결과 저장
      output_path <- "/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs"
      qs::qsave(res_nebula2_pb, output_path)
      message(sprintf("\n결과 저장: %s", output_path))
    }
    
  }, error = function(e) {
    message(sprintf("\nrunNEBULA2_v1_with_pseudobulk 오류: %s", conditionMessage(e)))
  })
} else {
  message("runNEBULA2_v1_with_pseudobulk 테스트 건너뜀 (test_nebula2_pb = FALSE)")
}

# ============================================================================
# 완료
# ============================================================================
message("\n========================================")
message("테스트 완료")
message("========================================")
