# LDS 함수 모듈화 가이드

## 개요

LDS 함수를 7개의 독립적인 단계별 함수로 모듈화하여, 각 단계를 개별적으로 실행하고 중간 결과를 저장할 수 있도록 개선했습니다.

## 모듈화의 장점

1. **디버깅 용이성**: 특정 단계에서 오류가 발생했을 때, 해당 단계만 수정하고 재실행 가능
2. **중간 결과 저장**: 각 단계의 결과를 저장하여 시간 절약 (특히 SVA, Dream 등 시간이 오래 걸리는 단계)
3. **유연성**: 필요한 단계만 선택적으로 실행 가능
4. **재사용성**: 각 단계 함수를 다른 분석에서도 활용 가능

## 모듈화된 함수 구조

### 단계별 함수

1. **`lds_01_extract_data()`**: 데이터 추출
   - Seurat/DGEList/Matrix에서 count 데이터와 메타데이터 추출
   - 입력 검증

2. **`lds_01b_filter_na()`**: NA 필터링
   - 포뮬러 변수의 NA 값을 가진 셀 제거
   - `remove_na` 파라미터로 제어

3. **`lds_02_parse_formula()`**: 포뮬러 파싱
   - LMM 포뮬러에서 고정 효과와 임의 효과 분리
   - `lme4::nobars()` 사용

4. **`lds_03_preprocess_dge()`**: DGEList 전처리
   - DGEList 생성
   - `filterByExpr`로 유전자 필터링 (설정 가능)
   - 정규화 인자 계산

5. **`lds_04_run_sva()`**: SVA 실행
   - Surrogate Variable Analysis 실행
   - SV 개수 자동 결정 (잔차 분산 기반)
   - SV를 메타데이터에 추가

6. **`lds_05_build_final_formula()`**: 최종 포뮬러 생성
   - 원본 포뮬러에 SV 추가

7. **`lds_06_run_dream()`**: Dream 실행
   - `voomWithDreamWeights` 실행
   - `dream`으로 LMM 피팅
   - `eBayes`로 Empirical Bayes 조정
   - Contrast 생성 및 적용 (p-value 계산을 위해)

8. **`lds_07_analyze_sva_correlation()`**: SVA 상관관계 분석
   - SVA와 메타데이터 변수 간 상관관계 계산
   - 선택적 실행 (`plot_sva_correlation` 파라미터)

### 통합 함수

- **`LDS()`**: 전체 파이프라인
  - 위의 모든 단계를 순차적으로 호출
  - 기존 인터페이스와 호환

## 사용 방법

### 방법 1: 전체 파이프라인 실행 (기존 방식)

```r
result <- LDS(
  sobj = is5s_test,
  formula = ~ g3_clean + (1|hos_no) + (1|GEM),
  save_intermediate = TRUE,  # 중간 결과 저장
  output_dir = "/data/user3/sobj/lds_intermediate",
  prefix = "lds_test"
)
```

### 방법 2: 단계별 실행 (권장)

```r
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

# ... 나머지 단계들
```

### 방법 3: 중간 결과 로드 및 재시작

특정 단계에서 오류가 발생했을 때:

```r
# 이전 단계의 결과 로드
step3 <- qs::qread("/data/user3/sobj/lds_intermediate/lds_test_03_preprocess_dge.qs")
step1b <- qs::qread("/data/user3/sobj/lds_intermediate/lds_test_01b_filter_na.qs")
step2 <- qs::qread("/data/user3/sobj/lds_intermediate/lds_test_02_parse_formula.qs")

# 문제가 발생한 단계부터 다시 시작
step4 <- lds_04_run_sva(
  dge = step3$dge,
  meta.data = step1b$meta.data,
  fixed_effects_formula = step2$fixed_effects_formula,
  ...
)
```

## 중간 결과 파일

`save_intermediate=TRUE`로 설정하면 다음 파일들이 생성됩니다:

- `{prefix}_01_extract_data.qs`: 데이터 추출 결과
- `{prefix}_01b_filter_na.qs`: NA 필터링 결과
- `{prefix}_02_parse_formula.qs`: 포뮬러 파싱 결과
- `{prefix}_03_preprocess_dge.qs`: DGEList 전처리 결과
- `{prefix}_04_run_sva.qs`: SVA 실행 결과
- `{prefix}_05_build_final_formula.qs`: 최종 포뮬러 생성 결과
- `{prefix}_06_run_dream.qs`: Dream 실행 결과
- `{prefix}_07_sva_correlation.qs`: SVA 상관관계 분석 결과

## 테스트 스크립트

모듈화된 함수를 테스트하는 스크립트:

```r
source("/home/user3/data_user3/git_repo/_wt/lds/scripts/lds/test_lds_modular.R")
```

이 스크립트는:
- 각 단계를 순차적으로 실행
- 중간 결과를 저장
- 각 단계의 성공/실패를 확인
- 최종 결과를 저장

## 디버깅 팁

1. **특정 단계에서 오류 발생 시**:
   - 해당 단계의 입력 데이터 확인
   - 이전 단계의 중간 결과 로드하여 검증
   - 해당 단계 함수만 수정하여 재실행

2. **시간이 오래 걸리는 단계**:
   - SVA (`lds_04_run_sva`)와 Dream (`lds_06_run_dream`) 단계는 시간이 오래 걸림
   - 중간 결과를 저장하여 재실행 시 시간 절약

3. **메모리 부족 시**:
   - 각 단계 실행 후 불필요한 객체 제거
   - 중간 결과 저장 후 R 세션 재시작

## 파일 위치

- 모듈화된 함수: `/home/user3/data_user3/git_repo/_wt/lds/myR/R/lds.R`
- 테스트 스크립트: `/home/user3/data_user3/git_repo/_wt/lds/scripts/lds/test_lds_modular.R`
- 기존 통합 테스트: `/home/user3/data_user3/git_repo/_wt/lds/scripts/lds/test_lds_interactive.R`

## 참고

- 기존 `LDS()` 함수는 여전히 사용 가능하며, 내부적으로 모듈화된 함수들을 호출합니다.
- 각 단계 함수는 독립적으로 사용 가능하지만, 입력/출력 형식에 주의해야 합니다.
- 자세한 사용법은 `TEST_INSTRUCTIONS.md`를 참조하세요.

