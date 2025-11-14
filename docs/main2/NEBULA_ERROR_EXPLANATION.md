# NEBULA 오류 원인 및 해결 방법

## 오류 메시지

```
Warning in nlminb(para[1], pql_nbm_gamma_ll, gamma = para[2], betas = betae,  :
  NA/NaN function evaluation

Error in value[[3L]](cond) : 
  NEBULA failed: task 7 failed - "objective in x0 returns NA". 
  Try reducing number of genes or checking data quality.
```

## 문제 원인

### 1. 완전 분리(Complete Separation) 문제

**GEM과 g3 변수 사이에 완전 분리가 발생했습니다:**

- **GEM2, GEM5, GEM8**: g3=2만 있고, g3=1이 없음
- **GEM6**: g3=1만 있고, g3=2가 없음

이로 인해 설계 행렬(`model.matrix(~ g3 + GEM)`)이 **특이(singular)**해집니다:
- 일부 GEM 레벨에서 g3의 모든 레벨이 없음
- 설계 행렬의 rank < columns
- NEBULA의 최적화 과정에서 NA/NaN 발생

### 2. 데이터 분포 확인

```r
# GEM x g3 contingency table
        g3
GEM     1   2
  GEM1 105 193  # ✓ 양쪽 모두 있음
  GEM2   0 327  # ✗ g3=1 없음
  GEM3 105  52  # ✓ 양쪽 모두 있음
  GEM4 209 208  # ✓ 양쪽 모두 있음
  GEM5   0 192  # ✗ g3=1 없음
  GEM6 342   0  # ✗ g3=2 없음
  GEM7  97 197  # ✓ 양쪽 모두 있음
  GEM8   0 320  # ✗ g3=1 없음
```

## 해결 방법

### 방법 1: GEM을 covar_effects에서 제거 (권장)

```r
neb1 <- runNEBULA2_v1(
  is5s, 
  fixed_effects = c("g3"), 
  covar_effects = NULL,  # GEM 제거
  patient_col = "hos_no", 
  offset = "nCount_RNA",
  min_count = 20  # 더 높은 min_count로 필터링 강화
)
```

**장점:**
- 간단하고 빠름
- 완전 분리 문제 해결
- g3의 효과만 분석

**단점:**
- GEM의 효과를 고려하지 않음

### 방법 2: 완전 분리된 GEM 레벨 제거

```r
# 완전 분리되지 않은 GEM만 사용
gem_levels_to_keep <- c("GEM1", "GEM3", "GEM4", "GEM7")
sobj_filtered <- is5s[, is5s@meta.data$GEM %in% gem_levels_to_keep & 
                       !is.na(is5s@meta.data$g3) & 
                       !is.na(is5s@meta.data$hos_no)]

neb1 <- runNEBULA2_v1(
  sobj_filtered, 
  fixed_effects = c("g3"), 
  covar_effects = "GEM",
  patient_col = "hos_no", 
  offset = "nCount_RNA",
  min_count = 20
)
```

**장점:**
- GEM의 효과를 고려할 수 있음
- 완전 분리 문제 해결

**단점:**
- 일부 GEM 레벨(GEM2, GEM5, GEM6, GEM8) 제거
- 샘플 수 감소

### 방법 3: GEM을 fixed_effects로 사용

```r
neb1 <- runNEBULA2_v1(
  is5s, 
  fixed_effects = c("GEM"),  # GEM을 fixed_effects로
  covar_effects = NULL,
  patient_col = "hos_no", 
  offset = "nCount_RNA",
  min_count = 20
)
```

**장점:**
- GEM의 효과를 분석할 수 있음
- 완전 분리 문제 없음

**단점:**
- g3의 효과를 고려하지 않음

### 방법 4: min_count를 높여서 유전자 수 제한

```r
neb1 <- runNEBULA2_v1(
  is5s, 
  fixed_effects = c("g3"), 
  covar_effects = "GEM",
  patient_col = "hos_no", 
  offset = "nCount_RNA",
  min_count = 50  # 더 높은 min_count
)
```

**장점:**
- 간단함

**단점:**
- 완전 분리 문제는 해결되지 않음
- 많은 유전자가 제외됨

## 개선 사항

`runNEBULA2_v1` 함수에 다음 기능을 추가했습니다:

1. **설계 행렬 특이성 확인**: `qr()`을 사용하여 rank 확인
2. **완전 분리 감지**: 변수 간 contingency table 확인
3. **자동 유전자 제한**: 설계 행렬이 특이하면 유전자 수를 1000개로 제한
4. **상세한 오류 메시지**: 문제 원인과 해결 방법 제시

## 권장 사항

1. **방법 1 사용** (GEM 제거): 가장 간단하고 안정적
2. **GEM의 효과가 중요하면 방법 2 사용** (완전 분리된 레벨 제거)
3. **분석 전에 데이터 분포 확인**: `table(sobj@meta.data$GEM, sobj@meta.data$g3)`

## 테스트 스크립트

문제를 진단하고 해결 방법을 테스트하려면:

```r
source("/home/user3/data_user3/git_repo/mylit-main2/test_nebula_issue.R")
```

## 참고

- **완전 분리(Complete Separation)**: 범주형 변수 간에 일부 조합이 완전히 없을 때 발생
- **설계 행렬 특이성(Singularity)**: 설계 행렬의 rank가 columns보다 작을 때 발생
- **NEBULA**: Negative Binomial mixed-effects model로, 설계 행렬이 특이하면 최적화가 실패할 수 있음

