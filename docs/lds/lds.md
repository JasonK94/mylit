# Limma-Dream-SVA (LDS) 파이프라인

## 개요

LDS (Limma-Dream-SVA)는 단일세포 RNA-seq 또는 GeoMx Digital Spatial Profiling 데이터에서 다중 임의 효과(예: 환자, 배치)를 포함한 선형혼합모형(LMM)을 피팅하는 통합 파이프라인입니다. SVA(Surrogate Variable Analysis)로 숨겨진 공변량을 탐지하고, limma의 `dream` 함수를 사용하여 차등 발현 분석을 수행합니다.

**주요 특징**:
- GeoMx 데이터와 같이 샘플 수가 많지 않은 경우에 적합
- 다중 임의 효과 모델링 (예: `(1|patient) + (1|batch)`)
- SVA를 통한 숨겨진 변동성 자동 탐지 및 보정
- 잔차 분산 기반 SV 개수 자동 결정

## 핵심 개념

### SVA (Surrogate Variable Analysis)

SVA는 알려지지 않은 공변량(예: 기술적 배치 효과, 숨겨진 생물학적 요인)을 탐지하는 방법입니다:

1. **고정 효과 모델 피팅**: 사용자가 지정한 고정 효과(예: treatment)로 모델 피팅
2. **잔차 계산**: 고정 효과를 제거한 잔차에서 숨겨진 패턴 탐지
3. **SV 추출**: SVD(Singular Value Decomposition)를 통해 잔차의 주요 변동성 방향 추출
4. **모델 보정**: 추출된 SV를 공변량으로 추가하여 모델 보정

### Limma-Dream

`dream`은 limma 패키지의 확장 기능으로, 다중 임의 효과를 포함한 LMM을 피팅할 수 있습니다:

- **voomWithDreamWeights**: 가중치를 고려한 voom 변환
- **dream**: LMM 피팅 (lme4 수식 지원)
- **eBayes**: Empirical Bayes 조정

### 파이프라인 워크플로우

```
입력 데이터 (Seurat/DGEList/Matrix)
    ↓
[1] 데이터 추출 및 검증
    ↓
[2] 포뮬러 파싱 (고정 효과 추출)
    ↓
[3] DGEList 생성, 필터링, 정규화
    ↓
[4] SVA 실행 (voom 변환 후)
    - 최대 SV 개수 탐지
    - 잔차 분산 기반 SV 개수 자동 결정 (또는 사용자 지정)
    ↓
[5] 최종 포뮬러 생성 (원본 + SV)
    ↓
[6] limma-dream 파이프라인
    - voomWithDreamWeights
    - dream (LMM 피팅)
    - eBayes (Empirical Bayes)
    ↓
최종 결과 (MArrayLM 객체)
```

## 주요 함수

### `LDS()`

전체 LDS 파이프라인을 실행하는 메인 함수입니다.

#### 파라미터

| 파라미터 | 설명 | 기본값 |
|---------|------|--------|
| `sobj` | Seurat 객체, DGEList 객체, 또는 Raw Count Matrix | 필수 |
| `formula` | 모델 포뮬러 (lme4 수식) | 필수 |
| `meta.data` | 메타데이터 (Matrix 입력 시 필수) | `NULL` |
| `layer` | Seurat 객체의 count layer | `"counts"` |
| `n_sv` | 사용할 SV 개수 | `NULL` (자동) |
| `sv_var_cutoff` | SV가 설명해야 할 잔차 분산 비율 | `0.5` |
| `n_cores` | 병렬 처리 CPU 코어 수 | `max(1, detectCores()-2)` |

#### 포뮬러 예시

```r
# 기본 예시: treatment 효과, 환자와 배치를 임의 효과로
~ treatment + (1|patient) + (1|batch)

# 복잡한 예시: 여러 고정 효과와 임의 효과
~ response + age + sex + (1|emrid) + (1|set)

# 임의 효과만 (고정 효과 없음)
~ (1|patient) + (1|batch)
```

#### 반환값

리스트 객체:
- `fit`: MArrayLM 객체 (limma의 `topTable()` 사용 가능)
- `voom`: 가중치가 포함된 EList 객체
- `sva_obj`: 원본 SVA 객체
- `svs_used`: 모델에 실제 사용된 SV 매트릭스
- `final_formula`: 최종 사용된 포뮬러 (문자열)
- `dge`: 필터링/정규화된 DGEList

#### 사용 예시

```r
# Seurat 객체에서 실행
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient) + (1|batch),
  n_sv = NULL,  # 자동 결정
  sv_var_cutoff = 0.5
)

# 결과 확인
top_genes <- limma::topTable(result$fit, number = 100)
head(top_genes)

# 사용된 SV 확인
head(result$svs_used)
result$final_formula
```

## SV 개수 결정 방법

### 1. 사용자 지정 (`n_sv` 제공)

```r
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient),
  n_sv = 5  # 정확히 5개 SV 사용
)
```

### 2. 자동 결정 (`n_sv = NULL`, 기본값)

잔차 분산의 누적 비율을 기반으로 자동 결정:

1. SVA로 최대 SV 개수 탐지
2. 각 SV가 설명하는 잔차 분산 비율 계산
3. 누적 분산이 `sv_var_cutoff` (기본 50%)를 넘는 최소 SV 개수 선택

```r
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient),
  n_sv = NULL,  # 자동 결정
  sv_var_cutoff = 0.5  # 잔차 분산의 50% 설명
)
```

**예시**: SV 1개가 30%, SV 2개가 55%를 설명하면 → 2개 선택

## 결과 해석

### `topTable()` 사용

```r
# 전체 결과 추출
all_results <- limma::topTable(
  result$fit,
  number = Inf,
  sort.by = "P"
)

# 특정 contrast에 대한 결과
# (dream의 경우 contrast를 지정해야 할 수 있음)
```

### 주요 컬럼

- `logFC`: 로그 fold change
- `AveExpr`: 평균 발현량
- `t`: t-statistic
- `P.Value`: p-value
- `adj.P.Val`: FDR 조정된 p-value
- `B`: log-odds

## 주의사항

### 1. 데이터 형식

- **Seurat 객체**: `layer = "counts"`로 raw count 추출
- **DGEList**: `$counts`와 `$samples` 사용
- **Matrix**: `meta.data` 필수 제공

### 2. 포뮬러 작성

- lme4 수식 문법 사용: `(1|group)` 형태
- 고정 효과와 임의 효과 구분
- 고정 효과만 없는 경우 (`~ (1|patient)`)도 가능하지만, SVA의 `mod0`는 `~ 1`로 설정됨

### 3. 샘플 수

- GeoMx 데이터처럼 샘플 수가 적은 경우에 적합
- 너무 적은 샘플 수에서는 SV 탐지가 어려울 수 있음

### 4. 패키지 의존성

필수 패키지:
- `limma` (dream, voomWithDreamWeights 포함)
- `edgeR` (DGEList, filterByExpr, calcNormFactors, voom)
- `lme4` (포뮬러 파싱)
- `BiocParallel` (병렬 처리)
- `sva` (SVA 실행)

### 5. 메모리 및 실행 시간

- 큰 데이터셋의 경우 메모리 사용량이 많을 수 있음
- SVA와 dream은 계산 집약적
- `n_cores` 조정으로 병렬 처리 속도 향상 가능

## 문제 해결

### SVA가 SV를 찾지 못함

```
... SVA가 유의미한 대리 변수(SV)를 찾지 못했습니다.
```

- 샘플 수가 너무 적을 수 있음
- 고정 효과 모델이 이미 대부분의 변동성을 설명할 수 있음
- 정상적인 경우일 수 있음 (SV 없이 진행)

### 모든 유전자가 필터링됨

```
모든 유전자가 필터링되었습니다.
```

- `filterByExpr` 조건이 너무 엄격함
- 데이터 품질 문제
- `design_for_filter` 확인 필요

### dream 피팅 실패

- 포뮬러 문법 확인
- 메타데이터에 필요한 변수가 모두 있는지 확인
- 샘플 수가 너무 적지 않은지 확인

## 참고 자료

- limma 패키지 문서: https://bioconductor.org/packages/limma
- dream 함수: `?limma::dream`
- SVA 논문: Leek et al. (2007) Nature Reviews Genetics
- 원본 구현: `myR/R/test_working.R` (lines 1554-1757)

