# Formula 1 분석 사용 가이드

## 개요

Formula 1: `~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)`

이 문서는 Formula 1을 사용한 NEBULA 분석을 실행하는 방법을 설명합니다.

## 파일 위치

- **전체 데이터 분석**: `/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_analysis.R`
- **경량 테스트**: `/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_formula1_light.R`
- **인터랙티브 스크립트**: `/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_interactive.R`

## 사용 방법

### 1. 경량 테스트 (권장: 먼저 실행)

정상 수행 여부를 빠르게 확인하기 위한 경량 테스트입니다.

```bash
cd /home/user3/GJC_KDW_250721
Rscript /home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_formula1_light.R
```

**특징:**
- 세포 수: 그룹당 최대 1000개 (총 약 2000개)
- 유전자 수: 1000개 (많이 발현된 유전자만)
- 예상 소요 시간: 수 분 ~ 10분 정도
- 결과 파일: `/data/user3/sobj/IS6_formula1_nebula_test_light_[timestamp].qs`

### 2. 전체 데이터 분석

전체 데이터셋으로 완전한 분석을 실행합니다.

```bash
cd /home/user3/GJC_KDW_250721
Rscript /home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_analysis.R
```

**특징:**
- 세포 수: 전체 (약 28,706개)
- 유전자 수: 전체 (약 51,795개, min_count=10 필터링 후)
- 예상 소요 시간: 수 시간 (데이터 크기에 따라)
- 결과 파일: `/data/user3/sobj/IS6_formula1_nebula_result_[timestamp].qs`

### 3. 인터랙티브 R 세션에서 사용

R 인터랙티브 세션에서 직접 실행할 수 있습니다.

#### 3.1 R 세션 시작

```bash
cd /home/user3/GJC_KDW_250721
R
```

#### 3.2 함수 로드

```r
# 인터랙티브 스크립트 소스
source("/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_interactive.R")
```

#### 3.3 분석 실행

**전체 데이터 분석:**
```r
result <- run_formula1_analysis_interactive()
```

**경량 테스트:**
```r
result <- run_formula1_analysis_interactive(light_test = TRUE)
```

**커스텀 옵션:**
```r
result <- run_formula1_analysis_interactive(
  data_path = "/data/user3/sobj/IS6_sex_added_251110.qs",
  output_dir = "/data/user3/sobj",
  min_count = 10,                    # 최소 발현 세포 수
  light_test = FALSE,                # 경량 테스트 여부
  max_genes = NULL,                  # 최대 유전자 수 (NULL = 전체)
  max_cells_per_group = NULL         # 그룹당 최대 세포 수 (NULL = 전체)
)
```

#### 3.4 결과 확인

```r
# 결과 구조 확인
str(result, max.level = 2)

# 결과 요약 확인
if (!is.null(result$summary)) {
  head(result$summary, 20)
  nrow(result$summary)  # 총 유전자 수
}

# Formula 확인
result$formula
result$design_formula
result$fixed_effects
```

## 분석 파라미터

### Formula 설명

```
~ g3 + sex + anno3.scvi + GEM + g3:anno3.scvi + sex:anno3.scvi + (1|GEM/hos_no)
```

- **g3**: 주요 그룹 변수 (IS/SAH 등)
- **sex**: 성별 (공변량)
- **anno3.scvi**: 셀 타입/클러스터
- **GEM**: 배치/실험 배치
- **g3:anno3.scvi**: g3와 셀 타입의 교호작용 (셀 타입별 g3 효과 차이)
- **sex:anno3.scvi**: 성별과 셀 타입의 교호작용 (셀 타입별 성별 효과 차이)
- **(1|GEM/hos_no)**: Nested random effects
  - `GEM/hos_no`는 `hos_no`가 `GEM` 내에 nesting됨을 의미
  - NEBULA에서는 `hos_no`만 random effect로, `GEM`은 fixed effect로 변환됨

### 주요 파라미터

- **min_count**: 최소 발현 세포 수 (기본값: 10)
  - 이 값 이상의 세포에서 발현된 유전자만 분석에 포함
  - 높을수록 더 적은 유전자가 포함되어 빠르지만, 정보 손실 가능

- **remove_na_cells**: NA 값이 있는 세포 제거 여부 (기본값: TRUE)
  - 모델에 사용되는 모든 변수에 NA가 없는 세포만 분석에 포함

- **layer**: 사용할 assay layer (기본값: "counts")
  - 원시 count 데이터 사용

## 결과 파일

### 결과 파일 구조

결과는 `.qs` 형식으로 저장되며, 다음 정보를 포함합니다:

- **result$summary**: 각 유전자별 분석 결과 (p-value, coefficient 등)
- **result$formula**: 사용된 원본 formula
- **result$design_formula**: 설계 행렬을 위한 formula
- **result$patient_col**: patient 컬럼명
- **result$fixed_effects**: 사용된 fixed effects 목록

### 결과 파일 로드

```r
library(qs)

# 결과 로드
result <- qs::qread("/data/user3/sobj/IS6_formula1_nebula_result_20241114_160000.qs")

# 결과 확인
str(result, max.level = 2)
head(result$summary, 20)
```

## 트러블슈팅

### 1. "Formula에 random effects가 없습니다" 에러

- **원인**: Formula 파싱 실패
- **해결**: Formula에 `(1|hos_no)` 또는 `(1|GEM/hos_no)` 형식의 random effects가 포함되어 있는지 확인

### 2. "설계 행렬이 특이(singular)합니다" 경고

- **원인**: 완전 분리(complete separation) 문제
- **해결**: 
  - 교호작용 항 제거 또는 단순화
  - 완전 분리된 조합 제거
  - `min_count`를 높여서 더 적은 유전자로 분석

### 3. 분석 시간이 너무 오래 걸림

- **해결**:
  - 경량 테스트로 먼저 확인 (`light_test = TRUE`)
  - `min_count`를 높여서 더 적은 유전자 분석
  - `max_genes` 또는 `max_cells_per_group`로 데이터 제한

### 4. 메모리 부족 에러

- **해결**:
  - 데이터 서브샘플링 사용
  - `min_count`를 높여서 유전자 수 감소
  - 클러스터별로 분리하여 분석

## 예시

### 경량 테스트 후 전체 분석

```bash
# 1. 경량 테스트 (빠른 확인)
cd /home/user3/GJC_KDW_250721
Rscript /home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_formula1_light.R

# 2. 결과 확인 후 문제없으면 전체 분석
Rscript /home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_analysis.R
```

### 인터랙티브 세션에서 단계별 실행

```r
# R 세션에서
source("/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_interactive.R")

# 1. 경량 테스트
result_test <- run_formula1_analysis_interactive(light_test = TRUE)

# 2. 결과 확인
str(result_test)
head(result_test$summary)

# 3. 문제없으면 전체 분석
result_full <- run_formula1_analysis_interactive()

# 4. 결과 비교 및 분석
```

## 참고사항

1. **분석 시간**: 전체 데이터 분석은 수 시간이 걸릴 수 있습니다.
2. **결과 저장**: 모든 분석 결과는 자동으로 `/data/user3/sobj/`에 저장됩니다.
3. **로그 확인**: 스크립트 실행 시 로그는 터미널에 출력됩니다.
4. **백그라운드 실행**: 전체 분석은 `nohup` 또는 `screen`을 사용하여 백그라운드로 실행하는 것을 권장합니다.

```bash
# 백그라운드 실행 예시
nohup Rscript /home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/run_formula1_analysis.R > formula1_analysis.log 2>&1 &
```

