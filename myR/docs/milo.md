# MiloR 차등 풍부도 분석 (Differential Abundance Analysis)

## 개요

MiloR는 단일세포 RNA 시퀀싱 데이터에서 세포 유형 또는 상태의 차등 풍부도(Differential Abundance, DA)를 분석하는 도구입니다. 이 패키지는 neighborhood 기반 접근법을 사용하여 공간적으로 인접한 세포 그룹의 풍부도 변화를 검출합니다.

## 주요 개념

### Neighborhood (이웃)
- 각 세포 주변의 k-nearest neighbor 세포들로 구성된 지역적 세포 집단
- `makeNhoods()` 함수로 생성되며, 각 neighborhood는 여러 세포를 포함할 수 있음
- Neighborhood 간에는 그래프 구조로 연결되어 있어 독립적이지 않음

### Block Method (블록 방법)
Neighborhood 간의 비독립성을 고려하기 위한 블록 생성 방법:

- **`"sample"`**: Cell name에서 sample ID를 추출하여 같은 sample에 속한 neighborhoods를 하나의 block으로 묶음
  - Cell name에서 `S\d+$` 패턴 또는 처음 8자를 추출
  - 각 neighborhood에 가장 많은 세포를 포함하는 sample을 할당
  - Sample-level 구조를 보존하여 permutation test의 타당성 확보

- **`"community"`**: nhoodGraph의 community detection (Louvain 알고리즘)을 사용하여 그래프 구조상 연결된 neighborhoods를 block으로 묶음
  - Sample 정보가 없거나 추출이 불가능한 경우 사용
  - 그래프의 연결 구조를 기반으로 한 pseudo-block 생성

- **`"none"`**: Blocking 없이 전체를 하나의 block으로 처리
  - 가장 단순한 방법이지만 비독립성 문제를 완전히 해결하지 못할 수 있음

## 함수

### `run_milo_pipeline`

MiloR 차등 풍부도 분석의 전체 파이프라인을 실행하는 통합 함수입니다.

#### 기능
- Seurat 객체를 Milo 객체로 변환
- kNN 그래프 구축 및 neighborhood 생성
- Neighborhood 간 거리 계산
- 차등 풍부도 검정 (GLM 기반)
- 결과 시각화 (UMAP, nhood graph, beeswarm plot)
- 중간 결과 캐싱 지원

#### 주요 파라미터

| 파라미터 | 설명 | 기본값 |
|---------|------|--------|
| `seurat_obj` | Seurat 객체 (메모리에 로드된 경우) | `NULL` |
| `seurat_qs_path` | Seurat 객체 `.qs` 파일 경로 | `NULL` |
| `patient_var` | 환자/샘플 식별자 컬럼명 | 필수 |
| `cluster_var` | 클러스터 식별자 컬럼명 | 필수 |
| `target_var` | 비교 대상 그룹 변수 (예: treatment) | 필수 |
| `batch_var` | 배치 효과 변수 | 필수 |
| `graph_reduction` | 그래프 구축에 사용할 차원 축소 | `"integrated.scvi"` |
| `layout_reduction` | 시각화에 사용할 차원 축소 | `"umap.scvi"` |
| `k` | kNN 그래프의 k 값 | `30` |
| `d` | 차원 축소 차원 수 | `30` |
| `prop` | Neighborhood 생성 비율 | `0.1` |
| `alpha` | 유의성 임계값 | `0.1` |
| `save` | 중간 결과 저장 여부 | `TRUE` |
| `output_dir` | 결과 저장 디렉터리 | `tempdir()/milo` |
| `prefix` | 파일명 접두사 | `"milo"` |
| `suffix` | 파일명 접미사 (자동 증가) | `NULL` |
| `force_run` | 캐시 무시하고 재실행 여부 | `FALSE` |
| `plotting` | 시각화 수행 여부 | `TRUE` |
| `max_cells` | 다운샘플링할 최대 세포 수 | `NULL` |

#### 반환값
- `milo`: Milo 객체
- `da_results`: 차등 풍부도 검정 결과 데이터프레임
- `plots`: 시각화 객체 리스트 (UMAP, nhood graph, beeswarm plot)

#### 사용 예시

```r
# Seurat 객체에서 직접 실행
result <- run_milo_pipeline(
    seurat_obj = seurat_object,
    patient_var = "patient_id",
    cluster_var = "seurat_clusters",
    target_var = "treatment",
    batch_var = "batch",
    graph_reduction = "integrated.scvi",
    layout_reduction = "umap.scvi",
    k = 30,
    d = 30,
    prop = 0.1,
    alpha = 0.1,
    save = TRUE,
    output_dir = "/path/to/output",
    prefix = "milo_analysis"
)

# .qs 파일에서 로드하여 실행
result <- run_milo_pipeline(
    seurat_qs_path = "/path/to/seurat.qs",
    patient_var = "patient_id",
    cluster_var = "seurat_clusters",
    target_var = "treatment",
    batch_var = "batch"
)
```

### `test_cluster_logfc_bias`

클러스터별 logFC 편중성을 검정하는 함수입니다. Neighborhood 간의 비독립성을 고려하여 여러 통계 검정 방법을 제공합니다.

#### 기능
- 클러스터별 평균 logFC 계산
- Block permutation test (권장)
- Correlation-adjusted t-test (neff 보정)
- Mixed-effects model (선택적, lme4 필요)
- Empirical Bayes 추정 (선택적, ashr 필요)

#### 주요 파라미터

| 파라미터 | 설명 | 기본값 |
|---------|------|--------|
| `da_results` | `miloR::testNhoods()` 결과 데이터프레임 | 필수 |
| `milo` | Milo 객체 | 필수 |
| `cluster_col` | 클러스터 식별자 컬럼명 | `"major_cluster"` |
| `block_method` | Block 생성 방법 (`"sample"`, `"community"`, `"none"`) | `"sample"` |
| `test_methods` | 실행할 검정 방법 | `c("permutation", "neff")` |
| `n_perm` | Permutation test 반복 횟수 | `2000` |
| `max_nhoods` | 다운샘플링할 최대 neighborhood 수 | `NULL` |
| `seed` | 난수 시드 | `1` |
| `verbose` | 진행 메시지 출력 여부 | `TRUE` |

#### 반환값
클러스터별 통계량이 포함된 데이터프레임:
- `major_cluster`: 클러스터 식별자
- `mean_logFC`: 클러스터별 평균 logFC
- `n_nhoods`: 클러스터에 속한 neighborhood 수
- `p_perm`: Block permutation test p-value (선택적)
- `p_neff`: Correlation-adjusted t-test p-value (선택적)
- `p_lmm`: Mixed-effects model p-value (선택적)
- `eb_mean`: Empirical Bayes posterior mean (선택적)

#### 사용 예시

```r
# 기본 사용 (permutation + neff)
cluster_bias <- test_cluster_logfc_bias(
    da_results = da_results,
    milo = milo,
    cluster_col = "major_cluster",
    block_method = "sample",
    test_methods = c("permutation", "neff"),
    n_perm = 2000
)

# 다운샘플링으로 빠른 테스트
cluster_bias_fast <- test_cluster_logfc_bias(
    da_results = da_results,
    milo = milo,
    block_method = "community",
    test_methods = c("permutation"),
    n_perm = 500,
    max_nhoods = 1000
)

# 모든 검정 방법 사용
cluster_bias_full <- test_cluster_logfc_bias(
    da_results = da_results,
    milo = milo,
    test_methods = c("permutation", "neff", "lmm", "ashr")
)
```

## 검정 방법 상세

### 1. Block Permutation Test (권장)
- **원리**: 각 block 내에서만 logFC 값을 섞어서 permutation을 수행
- **장점**: Block 구조를 보존하여 비독립성 문제를 완화
- **단점**: 계산 시간이 오래 걸림 (n_perm에 비례)

### 2. Correlation-adjusted t-test (neff)
- **원리**: Neighborhood 그래프의 인접 행렬을 기반으로 effective sample size (neff)를 계산하여 t-test의 자유도를 보정
- **장점**: 빠른 계산 속도
- **단점**: 근사적 방법이므로 permutation보다 덜 보수적일 수 있음

### 3. Mixed-Effects Model (LMM)
- **원리**: `logFC ~ 1 + (1 | block_id)` 형태의 선형혼합모형으로 block-level random effect를 모델링
- **장점**: Block 구조를 명시적으로 모델링
- **단점**: lme4 패키지 필요, 작은 클러스터에서는 수렴 실패 가능

### 4. Empirical Bayes (ashr)
- **원리**: p-value를 z-score로 변환한 후 Empirical Bayes shrinkage를 적용하여 effect size를 추정
- **장점**: 작은 effect도 안정적으로 추정 가능
- **단점**: ashr 패키지 필요, 검정보다는 effect size 추정에 유용

## 주의사항

### 1. Milo 객체 구조
- **중요**: `milo@Nhoods` (대문자)가 올바른 slot 이름입니다. `milo@nhoods` (소문자)는 존재하지 않습니다.
- `buildNhoodGraph()`를 호출한 후에만 `plotNhoodGraphDA()`를 사용할 수 있습니다.

### 2. Plotting 함수
- `plotDAbeeswarm()`는 `alpha`와 `SpatialFDR` 값에 민감합니다. 모든 값이 비유의적이면 오류가 발생할 수 있으므로, `beeswarm_metric`을 `"PValue"` 또는 `"logFC_percentile"`로 변경하거나 `beeswarm_alpha`를 조정하세요.

### 3. R 환경
- R은 반드시 `st/` 디렉터리에서 실행하여 `st/start.R`이 자동으로 소싱되도록 해야 합니다.
- 이렇게 해야 `renv` 환경이 올바르게 로드되고 필요한 패키지들이 사용 가능합니다.

### 4. Caching
- `save=TRUE` (기본값)일 때 중간 결과가 `.qs` 형식으로 저장됩니다.
- `force_run=FALSE`일 때는 기존 캐시 파일이 있으면 자동으로 재사용하여 시간을 절약합니다.
- 각 단계(`nhoods`, `distances`, `testing`)별로 캐시 파일이 생성됩니다.

## 워크플로우 예시

```r
# 1. Milo 파이프라인 실행
milo_result <- run_milo_pipeline(
    seurat_qs_path = "/path/to/seurat.qs",
    patient_var = "patient_id",
    cluster_var = "seurat_clusters",
    target_var = "treatment",
    batch_var = "batch",
    save = TRUE,
    output_dir = "/path/to/output"
)

# 2. 클러스터별 편중성 검정
cluster_bias <- test_cluster_logfc_bias(
    da_results = milo_result$da_results,
    milo = milo_result$milo,
    cluster_col = "major_cluster",
    block_method = "sample",
    test_methods = c("permutation", "neff"),
    n_perm = 2000
)

# 3. 결과 확인
print(cluster_bias)
# 유의한 클러스터 필터링
significant_clusters <- cluster_bias %>%
    filter(p_perm < 0.05 | p_neff < 0.05)
```

## 참고 자료

- MiloR 공식 문서: https://www.bioconductor.org/packages/release/bioc/html/miloR.html
- 원본 논문: Dann et al. (2022) Nature Biotechnology
- 개발 과정의 주요 이슈와 해결 방법은 `myR/docs/DEVLOG_Korean.md`와 `myR/docs/context_Korean.md`를 참고하세요.

