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

SVA는 알려지지 않은 공변량(예: 기술적 배치 효과, 숨겨진 생물학적 요인)을 탐지하는 방법입니다.

#### SVA의 수학적 원리

SVA는 **잔차 분산(residual variance)**을 설명하는 주성분을 찾습니다:

1. **고정 효과 모델 피팅**: 사용자가 지정한 고정 효과(예: `g3_clean + (1|hos_no) + (1|GEM)`)로 모델 피팅
   ```
   Y (genes × cells) = X (cells × covariates) × β + R (genes × cells)
   ```
   여기서 R은 잔차 행렬입니다.

2. **잔차 계산**: 고정 효과를 제거한 잔차 계산
   ```
   R = Y - X × β^T
   ```
   각 유전자별로 회귀 모델 피팅 후 잔차 계산

3. **잔차 행렬의 SVD**: 잔차 행렬 전체의 구조를 분석
   ```
   R = U × D × V^T
   ```
   - U: 유전자 공간의 주성분 (genes × SV)
   - D: 특이값 (SV의 중요도)
   - V: 셀 공간의 주성분 (cells × SV) ← 이것이 SV!

4. **SV 선택**: 각 SV가 설명하는 잔차 분산 비율 계산
   ```
   분산 설명 비율 = (D[i]^2) / sum(D^2)
   ```
   - 누적 분산이 `sv_var_cutoff`를 넘는 최소 SV 개수 선택
   - 또는 `sv_var_cutoff = NULL`이면 전체 유의미한 SV 사용

5. **모델 보정**: 추출된 SV를 공변량으로 추가하여 모델 보정

**중요한 점**:
- SVA는 **전체 291개 유전자의 총 잔차 분산**을 설명하는 것이 아니라
- **고정 효과로 설명되지 않는 잔차 분산**을 설명합니다
- 이 잔차 분산에는 배치 효과, 숨겨진 공변량 등이 포함될 수 있습니다

#### 유의미한 SV의 기준

SVA 패키지의 `sva()` 함수가 `n.sv = NULL`일 때 자동으로 결정합니다:

- **Permutation Test 기반**: `mod` (full model)과 `mod0` (null model)의 잔차를 비교하여 통계적으로 유의미한 SV 탐지
- **잔차 분산 기반**: 잔차 행렬의 SVD에서 특이값이 큰 순서로 SV 선택
- 일반적으로 **p-value < 0.05** 기준으로 유의미한 SV 결정
- 이후 `sv_var_cutoff`로 추가 필터링하여 실제 사용 개수 결정

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
| `sv_var_cutoff` | SV가 설명해야 할 잔차 분산 비율. `NULL`이면 전체 유의미한 SV 사용 | `0.5` |
| `n_cores` | 병렬 처리 CPU 코어 수 | `max(1, detectCores()-2)` |
| `remove_na` | NA 값 필터링 여부 | `TRUE` |
| `min.count` | `filterByExpr`의 `min.count` | `10` |
| `min.total.count` | `filterByExpr`의 `min.total.count` | `15` |
| `min.prop` | `filterByExpr`의 `min.prop` | `0.1` |
| `large.n` | `filterByExpr`의 `large.n` | `10` |
| `plot_sva_correlation` | SVA 상관관계 분석 및 Heatmap 생성 여부 | `TRUE` |

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
- `fit`: MArrayLM 객체 (limma의 `topTable()` 사용 가능, p-value 포함)
- `voom`: 가중치가 포함된 EList 객체
- `sva_obj`: 원본 SVA 객체
- `svs_used`: 모델에 실제 사용된 SV 매트릭스
- `final_formula`: 최종 사용된 포뮬러 (문자열)
- `dge`: 필터링/정규화된 DGEList
- `fixed_vars`: 고정 효과 변수 목록 (topTable에서 coef로 사용)
- `contrast_applied`: Contrast 적용 여부
- `sv_correlation`: SVA와 메타데이터 간 상관관계 행렬 (선택적)
- `heatmaps`: 생성된 Heatmap 파일 경로 및 상관관계 행렬 (선택적)

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

#### 2-1. `sv_var_cutoff` 지정 (기본값 0.5)

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

#### 2-2. `sv_var_cutoff = NULL` (전체 SV 사용)

모든 유의미한 SV를 사용:

```r
result <- LDS(
  sobj = seurat_obj,
  formula = ~ treatment + (1|patient),
  n_sv = NULL,
  sv_var_cutoff = NULL  # 전체 유의미한 SV 사용
)
```

## 결과 해석

### `topTable()` 사용

```r
# 고정 효과 변수 확인
result$fixed_vars

# coef를 지정하여 p-value 계산
all_results <- limma::topTable(
  result$fit,
  coef = result$fixed_vars[1],  # 첫 번째 고정 효과 변수
  number = Inf,
  sort.by = "P"
)

# 또는 coef 없이 (p-value가 이미 계산된 경우)
all_results <- limma::topTable(
  result$fit,
  number = Inf,
  sort.by = "P"
)
```

**참고**: `dream` 결과의 경우, `LDS()` 함수 내부에서 t-statistics와 df를 사용하여 p-value를 자동 계산합니다.

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

## SVA 상관관계 분석 및 시각화

### 자동 Heatmap 생성

`plot_sva_correlation = TRUE` (기본값)일 때, LDS 함수는 자동으로 3가지 Heatmap을 생성합니다:

1. **전체 상관관계 행렬** (`{prefix}_heatmap_full.png`)
   - 메타데이터와 SV를 포함한 전체 상관관계 행렬
   - 크기: (metadata + SV) × (metadata + SV)

2. **메타데이터 × SV** (`{prefix}_heatmap_metadata_x_sv.png`)
   - 모든 메타데이터 변수와 SV 간 상관관계
   - 행: 메타데이터 변수, 열: SV

3. **상위 n개 메타데이터 × SV** (`{prefix}_heatmap_top10_metadata_x_sv.png`)
   - p-value < 0.05인 메타데이터 변수 중 상위 10개
   - p-value를 별표로 표시 (* p<0.05, ** p<0.01, *** p<0.001)

### Helper 함수

#### `lds_corrplot()`

SVA와 메타데이터 간 상관관계 플롯을 생성하는 독립 함수입니다.

**차이점**: `lds_08_create_heatmaps()`는 LDS 파이프라인의 일부로 자동 실행되지만, `lds_corrplot()`은 독립적으로 사용할 수 있는 helper 함수입니다.

```r
# 독립적으로 사용
lds_corrplot(
  sva = result$svs_used,
  meta.data = meta.data,
  method = "spearman",
  display.value = TRUE,
  output_file = "correlation_plot.png"
)
```

#### `lds_corrplot_summary()`

상관관계 요약 테이블을 반환합니다:

```r
summary_table <- lds_corrplot_summary(
  sva = result$svs_used,
  meta.data = meta.data,
  top_n = 20
)
```

### 파일 위치

생성된 Heatmap 파일은 `output_dir`에 저장됩니다:
- 기본값: `tempdir()` (R 세션의 임시 디렉터리)
- `save_intermediate = TRUE`일 때: 지정한 `output_dir`

```r
# 결과에서 파일 경로 확인
result$heatmaps$heatmap_files
```

## 참고 자료

- limma 패키지 문서: https://bioconductor.org/packages/limma
- dream 함수: `?limma::dream`
- SVA 논문: Leek et al. (2007) Nature Reviews Genetics
- 모듈화된 구현: `myR/R/lds.R`
- Helper 함수: `myR/R/lds_corrplot.R`, `myR/R/lds_08_heatmaps.R`

