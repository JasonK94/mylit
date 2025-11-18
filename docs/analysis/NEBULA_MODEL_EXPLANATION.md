# NEBULA 모델 상세 설명: 가정, 복잡성, 완전 분리 문제

## 1. NEBULA 모델 개요

### 1.1 모델 구조

NEBULA (Negative Binomial mixed-effects model)는 단일 세포 RNA 시퀀싱 데이터의 차등 발현 분석을 위한 통계 모델입니다.

**모델 수식:**

```
Y_ij ~ Negative Binomial(μ_ij, φ)

log(μ_ij) = log(offset_ij) + X_ij^T β + b_i

where:
  Y_ij: 유전자 j의 세포 i에서의 count
  μ_ij: 유전자 j의 세포 i에서의 평균 count
  φ: Negative Binomial 분산 파라미터
  offset_ij: 세포 i에서의 offset (예: nCount_RNA)
  X_ij: 설계 행렬 (fixed effects)
  β: fixed effects 계수 (추정 대상)
  b_i: random effect (patient/sample i의 랜덤 절편)
      b_i ~ N(0, σ²_b)  (랜덤 효과 분산)
```

### 1.2 모델 구성 요소

1. **Fixed Effects (고정 효과)**
   - `fixed_effects`: 주요 관심 변수 (예: g3)
   - `covar_effects`: 보정하고 싶은 공변량 (예: GEM, batch)
   - 설계 행렬 X에 포함됨
   - 각 유전자별로 β 계수를 추정

2. **Random Effects (랜덤 효과)**
   - `patient_col` (예: hos_no): Patient/sample-level random intercept
   - 각 patient마다 다른 랜덤 효과 b_i를 가짐
   - b_i ~ N(0, σ²_b): 평균 0, 분산 σ²_b인 정규분포
   - σ²_b는 데이터에서 추정됨

3. **Offset**
   - `offset` (예: nCount_RNA): 세포별 sequencing depth
   - log(μ_ij) = log(offset_ij) + ... 형태로 사용
   - 각 세포의 total count로 정규화

## 2. 완전 분리(Complete Separation) 문제

### 2.1 완전 분리의 정의

완전 분리는 범주형 변수 간에 일부 조합이 완전히 없을 때 발생합니다.

**예시:**
```
GEM x g3 contingency table:
        g3
GEM     1   2
  GEM1 105 193  # ✓ 양쪽 모두 있음
  GEM2   0 327  # ✗ g3=1 없음 (완전 분리)
  GEM6 342   0  # ✗ g3=2 없음 (완전 분리)
```

### 2.2 완전 분리가 문제가 되는 경우

#### Case 1: Fixed Effects 간 완전 분리 (문제!)

**예시: GEM과 g3 (둘 다 fixed effects)**

```
model.matrix(~ g3 + GEM)

설계 행렬이 특이(singular)해짐:
- rank < columns
- 역행렬 계산 불가
- NEBULA 최적화 실패 (NA/NaN 발생)
```

**왜 문제인가?**
- Fixed effects는 설계 행렬 X에 포함됨
- X가 특이하면 β 계수를 추정할 수 없음
- NEBULA의 최적화 과정에서 NA/NaN 발생

#### Case 2: Fixed Effects와 Random Effects 간 완전 분리 (문제 없음!)

**예시: hos_no와 g3 (hos_no는 random effect, g3는 fixed effect)**

```
각 hos_no는 g3=1 또는 g3=2 중 하나만 가짐:
  hos_no 10956946: g3=2만 있음
  hos_no 13762164: g3=1만 있음

하지만 이것은 문제가 되지 않습니다!
```

**왜 문제가 되지 않는가?**
1. **Random effects는 다른 구조를 가짐**
   - Random effects는 평균 0인 정규분포를 따름
   - 분산 σ²_b를 데이터에서 추정
   - 설계 행렬 X에 직접 포함되지 않음

2. **Random effects는 patient-level intercept**
   - 각 patient마다 다른 랜덤 효과를 가짐
   - Patient 내부에서만 작동 (within-patient variation)
   - Patient 간 변이를 모델링 (between-patient variation)

3. **Fixed effects는 patient를 넘어서 작동**
   - g3 효과는 모든 patient에 공통
   - Patient 간 비교를 통해 추정
   - 각 patient가 g3=1 또는 g3=2 중 하나만 가져도 문제 없음

**수식으로 설명:**
```
log(μ_ij) = log(offset_ij) + X_ij^T β + b_i

where:
  X_ij: fixed effects (g3, GEM 등)
  b_i: random effect (hos_no i의 랜덤 절편)

각 patient i는:
  - g3=1 또는 g3=2 중 하나만 가짐
  - 하지만 b_i는 patient 내부의 모든 세포에 공통
  - g3 효과는 patient 간 비교를 통해 추정됨

예시:
  Patient 1 (g3=1): log(μ_1j) = log(offset_1j) + β_g3=1 + b_1
  Patient 2 (g3=2): log(μ_2j) = log(offset_2j) + β_g3=2 + b_2

Patient 1과 Patient 2를 비교하면:
  β_g3=2 - β_g3=1 = (log(μ_2j) - log(offset_2j) - b_2) - 
                     (log(μ_1j) - log(offset_1j) - b_1)

이렇게 patient 간 비교를 통해 g3 효과를 추정할 수 있음!
```

## 3. GEM을 covar_effects에 포함할지 여부

### 3.1 GEM을 포함하지 않는 경우 (권장)

```r
neb1 <- runNEBULA(
  is5s, 
  fixed_effects = c("g3"), 
  covar_effects = NULL,  # GEM 제거
  patient_col = "hos_no", 
  offset = "nCount_RNA"
)
```

**장점:**
1. **완전 분리 문제 해결**
   - GEM과 g3 간 완전 분리 문제 없음
   - 설계 행렬이 특이하지 않음
   - 분석이 안정적으로 수행됨

2. **모델 단순화**
   - 모델 복잡도 감소
   - 추정 파라미터 수 감소
   - 계산 속도 향상
   - 과적합(overfitting) 위험 감소

3. **g3 효과에 집중**
   - 주요 관심 변수인 g3 효과를 명확하게 추정
   - GEM의 혼란 효과(confounding) 제거

**단점:**
1. **GEM 효과를 고려하지 않음**
   - GEM이 유전자 발현에 영향을 미치는 경우
   - GEM 효과를 보정하지 못함

### 3.2 GEM을 포함하는 경우

```r
neb1 <- runNEBULA(
  is5s, 
  fixed_effects = c("g3"), 
  covar_effects = "GEM",  # GEM 포함
  patient_col = "hos_no", 
  offset = "nCount_RNA"
)
```

**문제점:**
1. **완전 분리 문제**
   - GEM과 g3 간 완전 분리
   - 설계 행렬이 특이해짐
   - 분석 실패

2. **모델 복잡도 증가**
   - 추가 파라미터 추정 필요
   - 계산 비용 증가
   - 과적합 위험 증가

**해결 방법:**
1. **완전 분리되지 않은 GEM만 사용**
   ```r
   gem_levels_to_keep <- c("GEM1", "GEM3", "GEM4", "GEM7")
   sobj_filtered <- is5s[, is5s@meta.data$GEM %in% gem_levels_to_keep]
   ```

2. **GEM을 fixed_effects로 사용 (g3 제거)**
   ```r
   fixed_effects = c("GEM"), 
   covar_effects = NULL
   ```

### 3.3 권장 사항

**GEM을 covar_effects에 포함하지 않는 것을 권장합니다.**

**이유:**
1. **완전 분리 문제 해결**
   - 분석이 안정적으로 수행됨
   - 오류 없이 결과를 얻을 수 있음

2. **모델 단순화**
   - 해석이 용이함
   - 계산 속도 향상
   - 과적합 위험 감소

3. **g3 효과에 집중**
   - 주요 관심 변수인 g3 효과를 명확하게 추정
   - GEM의 혼란 효과 제거

**GEM 효과가 중요한 경우:**
- 완전 분리되지 않은 GEM만 사용
- 또는 GEM을 fixed_effects로 사용 (g3 제거)
- 또는 별도의 분석으로 GEM 효과 분석

## 4. NEBULA 모델의 가정

### 4.1 주요 가정

1. **Negative Binomial 분포**
   - Count 데이터가 Negative Binomial 분포를 따름
   - Overdispersion을 모델링할 수 있음
   - 단일 세포 데이터에 적합

2. **Log-linear 모델**
   - log(μ_ij) = log(offset_ij) + X_ij^T β + b_i
   - 선형 관계 가정
   - Multiplicative effects

3. **Random Effects 가정**
   - b_i ~ N(0, σ²_b): 정규분포
   - Patient 간 독립성
   - Patient 내부 상관관계 고려

4. **Fixed Effects 가정**
   - 설계 행렬 X가 full rank (특이하지 않음)
   - 완전 분리 문제 없음
   - 선형 독립성

### 4.2 가정 위반 시 문제

1. **설계 행렬 특이성**
   - Fixed effects 간 완전 분리
   - 계수 추정 불가
   - NA/NaN 발생

2. **Random Effects 가정 위반**
   - Patient 간 상관관계
   - 정규분포 가정 위반
   - 결과 해석 어려움

3. **Negative Binomial 가정 위반**
   - 다른 분포 (예: Zero-inflated)
   - 모델 적합도 저하
   - 결과 신뢰도 감소

## 5. 모델 복잡성

### 5.1 복잡도 요소

1. **Fixed Effects 수**
   - 더 많은 fixed effects → 더 복잡한 모델
   - 추정 파라미터 수 증가
   - 계산 비용 증가

2. **Random Effects 수**
   - 더 많은 patient → 더 많은 random effects
   - 분산 추정 복잡도 증가
   - 계산 비용 증가

3. **유전자 수**
   - 더 많은 유전자 → 더 많은 분석
   - 계산 비용 선형 증가
   - 메모리 사용량 증가

### 5.2 복잡도 관리

1. **Fixed Effects 최소화**
   - 필요한 변수만 포함
   - 완전 분리 문제 확인
   - 변수 선택 신중히

2. **Random Effects 구조**
   - Patient-level random intercept만 사용
   - 복잡한 random effects 구조 피함
   - 계산 비용 고려

3. **유전자 필터링**
   - min_count로 필터링
   - 저발현 유전자 제거
   - 계산 비용 감소

## 6. 실전 권장 사항

### 6.1 모델 설정

```r
# 권장 설정
neb1 <- runNEBULA(
  is5s, 
  fixed_effects = c("g3"),           # 주요 관심 변수만
  covar_effects = NULL,              # 완전 분리 문제 있는 변수 제외
  patient_col = "hos_no",            # Random effect
  offset = "nCount_RNA",             # Sequencing depth
  min_count = 20                     # 적절한 필터링
)
```

### 6.2 문제 해결 체크리스트

1. **완전 분리 확인**
   - Fixed effects 간 contingency table 확인
   - 완전 분리 발견 시 변수 제외

2. **설계 행렬 특이성 확인**
   - rank < columns 확인
   - 특이 행렬 발견 시 변수 조정

3. **Random Effects 확인**
   - Patient 수 확인 (최소 10개 이상 권장)
   - Patient별 세포 수 확인

4. **유전자 필터링**
   - min_count 적절히 설정
   - 저발현 유전자 제거
   - 계산 비용 고려

### 6.3 결과 해석

1. **Fixed Effects 해석**
   - β_g3: g3 효과
   - Patient 간 비교를 통한 추정
   - 통계적 유의성 확인

2. **Random Effects 해석**
   - σ²_b: Patient 간 변이
   - 큰 σ²_b: Patient 간 변이가 큼
   - 작은 σ²_b: Patient 간 변이가 작음

3. **모델 적합도**
   - Residual 확인
   - 모델 가정 검증
   - 결과 신뢰도 평가

## 7. 요약

### 7.1 핵심 포인트

1. **완전 분리 문제**
   - Fixed effects 간 완전 분리: 문제!
   - Fixed effects와 random effects 간 완전 분리: 문제 없음

2. **GEM을 covar_effects에 포함할지 여부**
   - 권장: 포함하지 않음
   - 이유: 완전 분리 문제 해결, 모델 단순화

3. **hos_no와 g3 간 완전 분리**
   - 문제 없음
   - Random effects와 fixed effects 간 완전 분리는 정상
   - Patient 간 비교를 통해 g3 효과 추정 가능

### 7.2 권장 설정

```r
# 가장 안정적인 설정
neb1 <- runNEBULA(
  is5s, 
  fixed_effects = c("g3"), 
  covar_effects = NULL,  # GEM 제외
  patient_col = "hos_no", 
  offset = "nCount_RNA",
  min_count = 20
)
```

### 7.3 참고 자료

- NEBULA 패키지 문서
- Negative Binomial mixed-effects model
- Complete separation 문제
- Mixed-effects model 가정

