# LDS 함수 테스트 및 실행 방법

## 모듈화된 LDS 함수

LDS 함수는 이제 `lds.R` 파일에 모듈화되어 있습니다. 각 단계별 함수를 개별적으로 사용하거나, 전체 파이프라인을 실행할 수 있습니다.

### 모듈화된 함수 목록
- `lds_01_extract_data()`: 데이터 추출
- `lds_01b_filter_na()`: NA 필터링
- `lds_02_parse_formula()`: 포뮬러 파싱
- `lds_03_preprocess_dge()`: DGEList 전처리
- `lds_04_run_sva()`: SVA 실행
- `lds_05_build_final_formula()`: 최종 포뮬러 생성
- `lds_06_run_dream()`: Dream 실행
- `lds_07_analyze_sva_correlation()`: SVA 상관관계 분석
- `LDS()`: 전체 파이프라인 (위 함수들을 순차 호출)

### 중간 결과 저장

각 단계별 함수는 `save_intermediate=TRUE` 옵션으로 중간 결과를 저장할 수 있습니다. 이를 통해 특정 단계에서 오류가 발생했을 때, 이전 단계의 결과를 다시 로드하여 디버깅할 수 있습니다.

## 환경 설정

### 1. R 세션 시작
```bash
# alias st를 치면 들어가는 디렉터리에서 R 시작
cd /home/user3/GJC_KDW_250721
R
```

### 2. R 세션에서 실행
```r
# 워크트리에서 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/lds/myR")

# LDS 함수 확인
exists("LDS")
exists("lds_01_extract_data")  # 모듈화된 함수들도 확인
```

## 데이터 준비

### 1. 데이터 로드
```r
library(qs)

# 다운샘플링된 데이터 (테스트용)
is5s <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# 전체 데이터
is5 <- qs::qread("/data/user3/sobj/IS_scvi_251107.qs")
```

### 2. 메타데이터 확인
```r
# 메타데이터 컬럼 확인
colnames(is5s@meta.data)

# 샘플 ID 확인
table(is5s@meta.data$hos_no, useNA = "ifany")

# 그룹 변수 확인
table(is5s@meta.data$g3, useNA = "ifany")  # "NA", 1, 2

# 배치 변수 확인
table(is5s@meta.data$GEM, useNA = "ifany")

# 클러스터 확인
table(is5s@meta.data$anno3.scvi, useNA = "ifany")

# 임상 메타데이터 (환자별)
meta_clinical <- is5@meta.data %>%
  distinct(hos_no, .keep_all = TRUE)

# 환자-그룹 매핑 확인
table(meta_clinical$hos_no, meta_clinical$g3)
table(meta_clinical$GEM, meta_clinical$g3)
```

## LDS 테스트

### 모듈화된 함수로 단계별 테스트 (권장)

각 단계를 개별적으로 실행하고 중간 결과를 저장하여 디버깅이 용이합니다.

```r
# 테스트 스크립트 실행
source("/home/user3/data_user3/git_repo/_wt/lds/scripts/lds/test_lds_modular.R")
```

또는 수동으로 단계별 실행:

```r
# 데이터 로드
library(qs)
is5s <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
is5s_test <- is5s
is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]

# 중간 결과 저장 디렉터리
output_dir <- "/data/user3/sobj/lds_intermediate"
dir.create(output_dir, recursive = TRUE)

# 단계 1: 데이터 추출
step1 <- lds_01_extract_data(
  sobj = is5s_test,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 1b: NA 필터링
step1b <- lds_01b_filter_na(
  counts_matrix = step1$counts_matrix,
  meta.data = step1$meta.data,
  formula = ~ g3_clean + (1|hos_no) + (1|GEM),
  remove_na = TRUE,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 2: 포뮬러 파싱
step2 <- lds_02_parse_formula(
  formula = ~ g3_clean + (1|hos_no) + (1|GEM),
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 3: DGEList 전처리
step3 <- lds_03_preprocess_dge(
  counts_matrix = step1b$counts_matrix,
  meta.data = step1b$meta.data,
  fixed_effects_formula = step2$fixed_effects_formula,
  min.count = 5,
  min.prop = 0.05,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 4: SVA 실행
step4 <- lds_04_run_sva(
  dge = step3$dge,
  meta.data = step1b$meta.data,
  fixed_effects_formula = step2$fixed_effects_formula,
  n_sv = NULL,
  sv_var_cutoff = 0.5,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 5: 최종 포뮬러 생성
step5 <- lds_05_build_final_formula(
  original_formula = step2$original_formula,
  svs_final = step4$svs_final,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 6: Dream 실행
step6 <- lds_06_run_dream(
  dge = step3$dge,
  final_formula = step5$final_formula,
  meta.data = step4$meta.data_with_sv,
  n_cores = 4,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 단계 7: SVA 상관관계 분석
step7 <- lds_07_analyze_sva_correlation(
  svs_final = step4$svs_final,
  meta.data = step4$meta.data_with_sv,
  plot_sva_correlation = TRUE,
  save_intermediate = TRUE,
  output_dir = output_dir,
  prefix = "lds_test"
)

# 결과 확인
result <- list(
  fit = step6$fit_ebayes,
  voom = step6$v_dream,
  sva_obj = step4$sva_obj,
  svs_used = step4$svs_final,
  final_formula = step5$final_formula,
  dge = step3$dge
)
```

### 중간 결과 로드 및 재시작

특정 단계에서 오류가 발생했을 때, 이전 단계의 결과를 로드하여 해당 단계부터 다시 시작할 수 있습니다:

```r
# 예: 단계 4에서 오류 발생 → 단계 3까지의 결과 로드
step3 <- qs::qread("/data/user3/sobj/lds_intermediate/lds_test_03_preprocess_dge.qs")
step1b <- qs::qread("/data/user3/sobj/lds_intermediate/lds_test_01b_filter_na.qs")
step2 <- qs::qread("/data/user3/sobj/lds_intermediate/lds_test_02_parse_formula.qs")

# 단계 4부터 다시 시작
step4 <- lds_04_run_sva(
  dge = step3$dge,
  meta.data = step1b$meta.data,
  fixed_effects_formula = step2$fixed_effects_formula,
  ...
)
```

### 기본 테스트 (전체 파이프라인, 다운샘플링 데이터)

```r
# 변수 설정
sample_key <- "hos_no"
group_key <- "g3"
batch_key <- "GEM"

# g3에서 "NA" 제거 (필요시)
is5s_test <- is5s
is5s_test$g3_clean <- as.numeric(as.character(is5s_test$g3))
is5s_test <- is5s_test[, !is.na(is5s_test$g3_clean)]

# 또는 remove_na_groups 옵션이 있다면 사용
# (LDS 함수에는 없지만, 사전 필터링 필요)

# 기본 포뮬러: g3 효과, 환자와 배치를 임의 효과로
result_lds <- LDS(
  sobj = is5s_test,
  formula = ~ g3_clean + (1|hos_no) + (1|GEM),
  n_sv = NULL,  # 자동 결정
  sv_var_cutoff = 0.5,
  n_cores = 4,
  remove_na = TRUE,
  min.count = 5,
  min.prop = 0.05,
  plot_sva_correlation = TRUE,
  save_intermediate = TRUE,  # 중간 결과 저장
  output_dir = "/data/user3/sobj/lds_intermediate",
  prefix = "lds_test"
)

# 결과 확인
str(result_lds, max.level = 2)

# 사용된 SV 확인
if (!is.null(result_lds$svs_used)) {
  cat("사용된 SV 개수:", ncol(result_lds$svs_used), "\n")
  head(result_lds$svs_used)
}

# 최종 포뮬러 확인
cat("최종 포뮬러:", result_lds$final_formula, "\n")

# topTable로 결과 추출
library(limma)
top_genes <- topTable(
  result_lds$fit,
  number = 100,
  sort.by = "P"
)

head(top_genes)
```

### 전체 데이터 테스트

```r
# 전체 데이터로 실행 (시간이 오래 걸릴 수 있음)
is5_test <- is5
is5_test$g3_clean <- as.numeric(as.character(is5_test$g3))
is5_test <- is5_test[, !is.na(is5_test$g3_clean)]

result_lds_full <- LDS(
  sobj = is5_test,
  formula = ~ g3_clean + (1|hos_no) + (1|GEM),
  n_sv = NULL,
  sv_var_cutoff = 0.5,
  n_cores = 8  # 더 많은 코어 사용
)

# 결과 저장
qs::qsave(result_lds_full, "/data/user3/sobj/test_lds_result_full.qs")
```

### SV 개수 지정 테스트

```r
# SV 개수를 명시적으로 지정
result_lds_manual <- LDS(
  sobj = is5s_test,
  formula = ~ g3_clean + (1|hos_no) + (1|GEM),
  n_sv = 3,  # 정확히 3개 SV 사용
  n_cores = 4
)
```

### 다른 포뮬러 테스트

```r
# 성별을 공변량으로 추가
result_lds_covar <- LDS(
  sobj = is5s_test,
  formula = ~ g3_clean + sex + (1|hos_no) + (1|GEM),
  n_sv = NULL,
  sv_var_cutoff = 0.5
)

# 임의 효과만 (고정 효과 없음 - 주의 필요)
result_lds_re_only <- LDS(
  sobj = is5s_test,
  formula = ~ (1|hos_no) + (1|GEM),
  n_sv = NULL
)
```

## 결과 분석

### 1. 결과 요약
```r
# 전체 결과 추출
all_results <- topTable(
  result_lds$fit,
  number = Inf,
  sort.by = "P"
)

# 유의한 유전자 필터링
significant <- all_results %>%
  filter(adj.P.Val < 0.05)

cat("유의한 유전자 수:", nrow(significant), "\n")

# logFC 분포 확인
hist(all_results$logFC, breaks = 50)
```

### 2. 결과 저장
```r
# 전체 결과 저장
qs::qsave(result_lds, "/data/user3/sobj/test_lds_result.qs")

# topTable 결과 저장
write.csv(all_results, "/data/user3/sobj/test_lds_topTable.csv")
```

### 3. 결과 로드
```r
# 저장된 결과 로드
result_lds <- qs::qread("/data/user3/sobj/test_lds_result.qs")

# topTable 재계산
top_genes <- topTable(result_lds$fit, number = 100)
```

## 문제 해결

### 패키지 오류
```r
# 필수 패키지 확인
required_packages <- c("limma", "edgeR", "lme4", "BiocParallel", "sva")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("패키지 설치 필요:", pkg, "\n")
    if (pkg == "sva") {
      BiocManager::install("sva")
    } else {
      BiocManager::install(pkg)
    }
  }
}
```

### 함수가 로드되지 않음
```r
# 워크트리에서 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/lds/myR")

# 또는 함수 소스 직접 로드
source("/home/user3/data_user3/git_repo/_wt/lds/myR/R/lds.R")
```

### 데이터 로드 오류
```r
# 데이터 파일 확인
file.exists("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
file.exists("/data/user3/sobj/IS_scvi_251107.qs")

# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
```

### 메모리 부족
```r
# 작은 데이터셋으로 먼저 테스트
# 또는 유전자 서브셋 사용
# (LDS 함수 내부에서 filterByExpr로 이미 필터링됨)
```

## 다음 단계

1. 테스트 성공 후 전체 데이터로 확장
2. 다양한 포뮬러 테스트
3. 결과 해석 및 시각화
4. 다른 분석 방법과 비교 (예: MUSCAT, NEBULA)

