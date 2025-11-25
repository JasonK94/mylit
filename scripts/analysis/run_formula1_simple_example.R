# ============================================================================
# Formula 1 간단 실행 예시 (인터랙티브 세션용)
# ============================================================================
# R 인터랙티브 세션에서 하나씩 실행할 수 있는 간단한 예시
# ============================================================================

# 1. 데이터 로드
library(qs)
is6 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# 2. 함수 로드 (main2 워크트리에서)
source("/home/user3/data_user3/git_repo/_wt/analysis/myR/R/test_analysis.R")

# 3. Formula 정의
formula1 <- ~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)

# 4. 분석 실행
result <- runNEBULA(
  sobj = is6,
  formula = formula1,
  patient_col = "hos_no",
  offset = "nCount_RNA",
  min_count = 5,  # 또는 10, 20 등
  remove_na_cells = TRUE,
  layer = "counts"
)

# 5. 결과 확인
str(result, max.level = 2)
head(result$summary, 20)

# 6. 결과 저장
qs::qsave(result, file = "/data/user3/sobj/IS6_formula1_result_interactive.qs")

