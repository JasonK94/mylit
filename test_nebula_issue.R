# ============================================================================
# NEBULA 오류 문제 진단 스크립트
# ============================================================================
# 사용자가 보고한 오류를 재현하고 진단합니다.
# ============================================================================

# 작업 디렉터리 설정
setwd("/home/user3/GJC_KDW_250721")
source("start.R")

# 함수 로드
source("/home/user3/data_user3/git_repo/mylit-main2/myR/R/test_analysis.R")

# 데이터 로드
is5s <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

message("========================================")
message("문제 진단: GEM과 g3의 완전 분리")
message("========================================")

# NA 제거 후 데이터 확인
sobj_clean <- is5s[, !is.na(is5s@meta.data$g3) & 
                    !is.na(is5s@meta.data$GEM) & 
                    !is.na(is5s@meta.data$hos_no)]

message("Clean cells:", ncol(sobj_clean))
message("\nGEM x g3 contingency table:")
print(table(sobj_clean@meta.data$GEM, sobj_clean@meta.data$g3))

message("\n완전 분리된 조합:")
contingency <- table(sobj_clean@meta.data$GEM, sobj_clean@meta.data$g3)
zero_cells <- which(contingency == 0, arr.ind = TRUE)
if (nrow(zero_cells) > 0) {
  for (i in 1:nrow(zero_cells)) {
    gem_level <- rownames(contingency)[zero_cells[i, 1]]
    g3_level <- colnames(contingency)[zero_cells[i, 2]]
    message(sprintf("  - GEM=%s, g3=%s: 0개 세포", gem_level, g3_level))
  }
}

message("\n========================================")
message("해결 방법 테스트")
message("========================================")

# 방법 1: GEM을 covar_effects에서 제거 (g3만 사용)
message("\n방법 1: GEM을 제거하고 g3만 사용")
tryCatch({
  neb1_test1 <- runNEBULA2_v1(
    is5s, 
    fixed_effects = c("g3"), 
    covar_effects = NULL,  # GEM 제거
    patient_col = "hos_no", 
    offset = "nCount_RNA",
    min_count = 20  # 더 높은 min_count
  )
  message("✓ 방법 1 성공: GEM을 제거하고 g3만 사용")
}, error = function(e) {
  message("✗ 방법 1 실패:", conditionMessage(e))
})

# 방법 2: GEM을 fixed_effects로 사용 (g3 제거)
message("\n방법 2: GEM을 fixed_effects로 사용 (g3 제거)")
tryCatch({
  neb1_test2 <- runNEBULA2_v1(
    is5s, 
    fixed_effects = c("GEM"),  # GEM을 fixed_effects로
    covar_effects = NULL,
    patient_col = "hos_no", 
    offset = "nCount_RNA",
    min_count = 20
  )
  message("✓ 방법 2 성공: GEM을 fixed_effects로 사용")
}, error = function(e) {
  message("✗ 방법 2 실패:", conditionMessage(e))
})

# 방법 3: 완전 분리된 GEM 레벨 제거
message("\n방법 3: 완전 분리된 GEM 레벨 제거 (GEM2, GEM5, GEM6, GEM8 제거)")
gem_levels_to_keep <- c("GEM1", "GEM3", "GEM4", "GEM7")
sobj_filtered <- is5s[, is5s@meta.data$GEM %in% gem_levels_to_keep & 
                       !is.na(is5s@meta.data$g3) & 
                       !is.na(is5s@meta.data$hos_no)]
message(sprintf("Filtered cells: %d", ncol(sobj_filtered)))
message("GEM x g3 after filtering:")
print(table(sobj_filtered@meta.data$GEM, sobj_filtered@meta.data$g3))

tryCatch({
  neb1_test3 <- runNEBULA2_v1(
    sobj_filtered, 
    fixed_effects = c("g3"), 
    covar_effects = "GEM",
    patient_col = "hos_no", 
    offset = "nCount_RNA",
    min_count = 20
  )
  message("✓ 방법 3 성공: 완전 분리된 GEM 레벨 제거")
}, error = function(e) {
  message("✗ 방법 3 실패:", conditionMessage(e))
})

message("\n========================================")
message("진단 완료")
message("========================================")

