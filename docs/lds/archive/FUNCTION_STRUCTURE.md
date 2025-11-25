# LDS 함수 구조 및 기능 상세 분석

## 함수 개요

**LDS (Limma-Dream-SVA)**는 세 가지 주요 기술을 결합한 통합 파이프라인입니다:
1. **SVA (Surrogate Variable Analysis)**: 숨겨진 공변량 탐지
2. **Dream (variancePartition)**: 다중 임의 효과를 포함한 LMM 피팅
3. **Limma**: Empirical Bayes 조정 및 차등 발현 분석

## 함수 시그니처

```r
LDS(
  sobj,                    # 입력 데이터
  formula,                 # LMM 포뮬러
  meta.data = NULL,        # 메타데이터
  layer = "counts",        # Seurat layer
  n_sv = NULL,             # SV 개수 (NULL=자동)
  sv_var_cutoff = 0.5,     # SV 분산 설명 비율
  n_cores = ...            # 병렬 처리 코어 수
)
```

## 7단계 파이프라인 상세 분석

### 단계 1: 입력값 검증 및 데이터 추출

**목적**: 다양한 입력 형식을 통일된 형식으로 변환

**처리 과정**:
```r
1. 입력 타입 확인:
   - Seurat 객체 → GetAssayData(sobj, layer = "counts")
   - DGEList → sobj$counts
   - Matrix → sobj (그대로 사용)

2. 메타데이터 추출:
   - Seurat: sobj@meta.data
   - DGEList: sobj$samples
   - Matrix: meta.data 파라미터 필수

3. 일치성 검증:
   - ncol(counts_matrix) == nrow(meta.data)
```

**핵심 포인트**:
- Seurat의 경우 `layer` 파라미터로 count 데이터 위치 지정 가능
- Matrix 입력 시 `meta.data` 필수 제공

### 단계 2: 포뮬러 파싱

**목적**: LMM 포뮬러에서 고정 효과와 임의 효과 분리

**처리 과정**:
```r
1. 고정 효과 추출:
   fixed_effects_formula <- lme4::nobars(formula)
   # 예: ~ treatment + (1|patient) → ~ treatment

2. 고정 효과가 없는 경우 처리:
   if (is.null(fixed_effects_formula)) {
     fixed_effects_formula <- ~ 1
   }
```

**이유**:
- SVA는 고정 효과 모델(`mod`)과 null 모델(`mod0`)이 필요
- `filterByExpr`도 고정 효과 기반 디자인 행렬 필요
- Dream은 전체 포뮬러(고정 + 임의 + SV) 사용

**예시**:
```r
입력: ~ g3_clean + (1|hos_no) + (1|GEM)
고정 효과: ~ g3_clean
임의 효과: (1|hos_no) + (1|GEM)
```

### 단계 3: DGEList 생성, 필터링, 정규화

**목적**: RNA-seq 데이터 전처리

**처리 과정**:
```r
1. DGEList 생성:
   dge <- edgeR::DGEList(counts_matrix, samples = meta.data)

2. 유전자 필터링:
   design_for_filter <- model.matrix(fixed_effects_formula, data = meta.data)
   keep_genes <- edgeR::filterByExpr(dge, design = design_for_filter)
   dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

3. 정규화:
   dge <- edgeR::calcNormFactors(dge)
```

**filterByExpr의 역할**:
- 너무 낮게 발현되는 유전자 제거
- 모든 그룹에서 발현되지 않는 유전자 제거
- 고정 효과 모델 기반으로 필터링 (임의 효과는 고려하지 않음)

**결과**:
- 예: 51,795개 유전자 → 55개 유전자 (다운샘플링 데이터 기준)

### 단계 4: SVA 실행

**목적**: 숨겨진 공변량(배치 효과, 숨겨진 생물학적 요인) 탐지

**처리 과정**:

#### 4-1. Voom 변환
```r
mod_sva <- model.matrix(fixed_effects_formula, data = meta.data)
mod0_sva <- model.matrix(~ 1, data = meta.data)
v_sva <- limma::voom(dge, mod_sva, plot = FALSE)
```

**이유**: SVA는 log2 변환된 데이터에 실행

#### 4-2. SVA 실행
```r
sva_obj <- sva::sva(
  v_sva$E,           # voom 변환된 발현 데이터
  mod = mod_sva,      # 고정 효과 모델
  mod0 = mod0_sva,    # null 모델 (절편만)
  n.sv = NULL         # 최대 SV 개수 자동 탐지
)
```

**SVA의 원리**:
1. 고정 효과 모델로 피팅: `Y ~ mod`
2. 잔차 계산: `residuals = Y - Y_predicted`
3. 잔차의 SVD: 잔차에서 주요 변동성 방향 추출
4. SV 추출: 주요 성분을 Surrogate Variable로 사용

#### 4-3. SV 개수 결정

**방법 1: 사용자 지정**
```r
if (!is.null(n_sv)) {
  n_sv_final <- min(n_sv, n_sv_max)
}
```

**방법 2: 자동 결정 (기본값)**
```r
# 잔차 행렬 계산
res_matrix <- t(resid(lm.fit(mod_sva, t(v_sva$E))))

# SVD
svd_res <- svd(res_matrix)

# 분산 비율 계산
percent_var_explained <- (svd_res$d^2) / sum(svd_res$d^2)
cumulative_var <- cumsum(percent_var_explained)

# Cutoff 만족하는 최소 SV 개수
n_sv_auto <- which(cumulative_var >= sv_var_cutoff)[1]
```

**예시**:
- SV 1개: 30% 설명
- SV 2개: 45% 설명
- SV 3개: 54.2% 설명
- `sv_var_cutoff = 0.5` → SV 3개 선택

#### 4-4. SV 메타데이터에 추가
```r
svs_final <- sva_obj$sv[, 1:n_sv_final, drop = FALSE]
colnames(svs_final) <- paste0("SV", 1:n_sv_final)
meta.data <- cbind(meta.data, svs_final)
```

### 단계 5: 최종 포뮬러 생성

**목적**: 원본 포뮬러에 SV 추가

**처리 과정**:
```r
original_formula_str <- paste(deparse(formula), collapse = "")

if (n_sv_final > 0) {
  sv_str <- paste(colnames(svs_final), collapse = " + ")
  final_formula_str <- paste(original_formula_str, sv_str, sep = " + ")
} else {
  final_formula_str <- original_formula_str
}

final_formula <- as.formula(final_formula_str)
```

**예시**:
```
입력: ~ g3_clean + (1|hos_no) + (1|GEM)
출력: ~ g3_clean + (1|hos_no) + (1|GEM) + SV1 + SV2 + SV3
```

**이유**: Dream 모델에서 SV를 공변량으로 포함하여 숨겨진 변동성 보정

### 단계 6: Limma-Dream 파이프라인

**목적**: 다중 임의 효과를 포함한 LMM 피팅

#### 6-1. VoomWithDreamWeights
```r
v_dream <- variancePartition::voomWithDreamWeights(
  dge, 
  final_formula, 
  meta.data, 
  BPPARAM = BPPARAM_SETUP
)
```

**역할**:
- 일반 `voom`과 달리 LMM 구조를 고려한 가중치 계산
- 임의 효과의 변동성을 반영한 가중치

**차이점**:
- `limma::voom`: 고정 효과만 고려
- `variancePartition::voomWithDreamWeights`: 임의 효과도 고려

#### 6-2. Dream (LMM 피팅)
```r
fit_dream <- variancePartition::dream(
  v_dream, 
  final_formula, 
  meta.data, 
  BPPARAM = BPPARAM_SETUP
)
```

**역할**:
- 각 유전자에 대해 LMM 피팅
- 포뮬러: `log2(expression) ~ g3_clean + SV1 + SV2 + SV3 + (1|hos_no) + (1|GEM)`
- 병렬 처리로 속도 향상

**모델 구조**:
```
Y_ij = β₀ + β₁*g3_clean + β₂*SV1 + β₃*SV2 + β₄*SV3 + 
       b_patient[i] + b_batch[i] + ε_ij

where:
  b_patient[i] ~ N(0, σ²_patient)
  b_batch[i] ~ N(0, σ²_batch)
  ε_ij ~ N(0, σ²_residual)
```

#### 6-3. Empirical Bayes 조정
```r
fit_ebayes <- limma::eBayes(fit_dream)
```

**역할**:
- 유전자 간 분산 정보 공유
- 작은 샘플 크기에서도 안정적인 통계량 추정
- Moderated t-statistics 계산

### 단계 7: 결과 반환

**반환 구조**:
```r
list(
  fit = fit_ebayes,           # MArrayLM2 객체
  voom = v_dream,              # EList 객체
  sva_obj = sva_obj,          # SVA 결과
  svs_used = svs_final,       # 사용된 SV 매트릭스
  final_formula = final_formula, # 최종 포뮬러
  dge = dge                   # DGEList
)
```

**사용 예시**:
```r
# 차등 발현 유전자 추출
top_genes <- limma::topTable(result$fit, number = 100)

# 특정 contrast
top_genes <- limma::topTable(
  result$fit, 
  coef = "g3_clean",  # g3_clean의 계수
  number = 100
)
```

## 핵심 설계 원칙

### 1. 유연한 입력 형식
- Seurat, DGEList, Matrix 모두 지원
- 사용자의 데이터 형식에 맞춰 자동 변환

### 2. 자동화된 SV 탐지
- 사용자가 SV 개수를 지정하지 않아도 자동 결정
- 잔차 분산 기반으로 통계적으로 의미있는 SV만 선택

### 3. 다중 임의 효과 지원
- 환자, 배치 등 여러 임의 효과 동시 모델링
- 반복 측정 설계에 적합

### 4. 숨겨진 변동성 보정
- SVA로 알려지지 않은 공변량 탐지
- Dream 모델에 SV를 공변량으로 포함하여 보정

## 데이터 흐름도

```
Raw Count Matrix (Seurat/DGEList/Matrix)
    ↓
[1] 데이터 추출
    ↓
Count Matrix + Meta.data
    ↓
[2] 포뮬러 파싱
    ↓
Fixed Effects + Random Effects
    ↓
[3] DGEList 생성
    ↓
Filtered & Normalized DGEList
    ↓
[4] Voom 변환 → SVA 실행
    ↓
Surrogate Variables (SV1, SV2, SV3, ...)
    ↓
[5] 포뮬러 확장
    ↓
Original Formula + SV
    ↓
[6] VoomWithDreamWeights → Dream → eBayes
    ↓
MArrayLM2 객체 (피팅 결과)
    ↓
[7] 결과 반환
```

## 주요 의존성

| 패키지 | 역할 |
|--------|------|
| `edgeR` | DGEList, 필터링, 정규화 |
| `limma` | voom 변환, eBayes 조정 |
| `sva` | Surrogate Variable Analysis |
| `variancePartition` | dream, voomWithDreamWeights |
| `lme4` | 포뮬러 파싱 |
| `BiocParallel` | 병렬 처리 |

## 사용 시나리오

### 시나리오 1: GeoMx 데이터
- 샘플 수가 적음 (10-50개)
- 여러 환자에서 반복 측정
- 배치 효과 존재
- **→ LDS 적합**: SVA로 숨겨진 변동성 보정, Dream으로 환자/배치 모델링

### 시나리오 2: 단일세포 데이터 (pseudobulk)
- 셀을 클러스터별로 집계
- 환자별, 배치별 구조
- **→ LDS 적합**: 임의 효과로 환자/배치 모델링

### 시나리오 3: 단순한 비교
- 두 그룹 비교만 필요
- 임의 효과 없음
- **→ LDS 불필요**: 일반 limma 사용

## 주의사항

1. **샘플 수**: 너무 적으면 SVA가 SV를 찾지 못할 수 있음
2. **유전자 필터링**: `filterByExpr`로 대부분의 유전자가 필터링될 수 있음
3. **메모리**: 큰 데이터셋에서는 메모리 사용량이 많을 수 있음
4. **실행 시간**: Dream은 계산 집약적 (병렬 처리로 완화)

## 성능 최적화

1. **병렬 처리**: `n_cores` 조정
2. **유전자 사전 필터링**: 분석 전에 저발현 유전자 제거
3. **SV 개수 제한**: `n_sv`로 SV 개수 제한
4. **다운샘플링**: 테스트 시 작은 데이터셋 사용

