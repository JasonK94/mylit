# Patient-Level Dimensionality Reduction and Anomaly Detection

## 개요

단일세포 RNA-seq 데이터에서 환자 수준의 차원 축소 및 이상 탐지를 위한 파이프라인입니다. 셀 레벨 데이터를 환자×클러스터 단위로 집계하여, 환자 간 유사성/차이를 분석하고 임상 변수(예: `g3`, `nih_change`)와의 연관성을 탐색합니다.

## 핵심 개념

### "뷰(View)"란?

환자 수준의 feature를 구성하는 **독립적인 정보 소스**를 의미합니다. 현재 파이프라인은 3가지 뷰를 지원하며, 각각을 선택적으로 포함할 수 있습니다:

1. **Frequency view**: 클러스터 빈도 (CLR 변환)
2. **Signature view**: 클러스터별 top marker 유전자들의 모듈 스코어
3. **Latent view**: 차원 축소 임베딩(예: scVI)의 클러스터별 평균

각 뷰는 서로 다른 측면의 정보를 담고 있으므로, 필요에 따라 일부만 선택하거나 가중치를 조정하여 결합할 수 있습니다.

### Feature 생성 과정

1. **Cluster frequency (CLR 변환)**
   - 셀 레벨 메타데이터 → 환자×클러스터 빈도 행렬
   - 조성 데이터 특성을 보존하기 위해 CLR(Centered Log-Ratio) 변환 적용
   - 기본 가중치: 1.5 (빈도 정보의 중요도 반영)

2. **Signature scores (클러스터별 top 20 유전자 모듈 스코어)**
   - `FindAllMarkers` 결과에서 클러스터별 상위 20개 유전자 선택
   - `Seurat::AddModuleScore()`로 셀 레벨 스코어 계산
   - 환자×클러스터 평균으로 집계
   - 기본 가중치: 1.0

3. **Latent embeddings (reduction의 클러스터별 평균)**
   - 기본 reduction: `integrated.scvi`
   - 각 클러스터 내 셀들의 잠재 임베딩을 환자×클러스터 평균으로 집계
   - 기본 가중치: 1.0

### 파이프라인 워크플로우

```
셀 레벨 데이터
    ↓
[뷰 선택] frequency / signatures / latent
    ↓
[블록 스케일링] 각 뷰를 독립적으로 Z-score 정규화
    ↓
[가중치 적용] frequency_weight × frequency + signature_weight × signatures + latent_weight × latent
    ↓
[결합] 환자 교집합만 유지하여 결합
    ↓
[배치 보정] (선택) limma::removeBatchEffect
    ↓
[PCA] 안전한 PC 수 선택 (elbow + rank fraction)
    ↓
[UMAP] PCA 공간에서 UMAP 임베딩
    ↓
[메타데이터 결합] g3, nih_change 등 임상 변수 병합
    ↓
최종 plot_df + embedding
```

## 주요 함수

### End-to-end 파이프라인

#### `patient_dimensionality_reduction()`

전체 워크플로우를 실행하는 메인 함수입니다.

**파라미터:**
- `include_frequency`, `include_signatures`, `include_latent`: 각 뷰 포함 여부 (기본값: 모두 TRUE)
- `frequency_weight`, `signature_weight`, `latent_weight`: 각 뷰의 가중치
- `batch_var`: 배치 보정 변수 (NULL이면 스킵)
- `reduction`: 잠재 임베딩 reduction 이름 (기본: `"integrated.scvi"`)

**반환값:**
- `embedding`: UMAP 좌표 행렬 (행=환자, 열=UMAP1/UMAP2)
- `plot_df`: 플로팅용 데이터프레임 (UMAP 좌표 + 메타데이터)
- `pca`: PCA 결과 객체
- `combined_features`: 결합된 feature 행렬 (PCA 입력)
- `view_columns`: 각 뷰의 컬럼명 리스트

### Feature 생성 함수들

#### `make_cluster_signatures()`
`FindAllMarkers` 결과를 `AddModuleScore`용 시그니처 리스트로 변환합니다.

#### `compute_patient_cluster_frequency()`
셀 레벨 메타데이터를 환자×클러스터 빈도 행렬로 변환합니다 (CLR 변환 옵션).

#### `compute_patient_signature_matrix()`
모듈 스코어를 환자×클러스터 평균으로 집계합니다.

#### `compute_patient_latent_matrix()`
잠재 임베딩을 환자×클러스터 평균으로 집계합니다.

### 유틸리티 함수들

#### `clr_transform()`
조성 데이터를 CLR 변환합니다 (pseudo-count 포함).

#### `combine_patient_feature_views()`
여러 뷰를 블록 스케일링 후 결합합니다.

#### `remove_batch_effect_if()`
선택적 배치 보정을 수행합니다 (`limma::removeBatchEffect`).

#### `choose_pca_k()`, `run_pca_safely()`
안전한 PC 수를 선택하고 PCA를 실행합니다 (데이터 크기에 따라 `irlba` 또는 `prcomp` 자동 선택).

### 시각화 및 보강

#### `plot_patient_umap()`
환자 UMAP 플로팅 함수입니다. `plot_embedding()`을 래핑하여 연속/이산 변수 모두 지원합니다.

#### `augment_plot_df()`
메타데이터를 보강합니다:
- ID 매칭 검증
- 중복/충돌 검사 (같은 ID에 서로 다른 값이 있는지)
- NA 카운트 리포트

## 사용 예시

### 기본 사용법

```r
# 1. FindAllMarkers 실행
markers_df <- FindAllMarkers(
  seurat_obj,
  group.by = "anno3.scvi",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# 2. 파이프라인 실행
res <- patient_dimensionality_reduction(
  seurat_obj = seurat_obj,
  markers_df = markers_df,
  reduction = "integrated.scvi",
  verbose = TRUE
)

# 3. 플로팅
plot_patient_umap(res$plot_df, color_by = "g3")
plot_patient_umap(res$plot_df, color_by = "nih_change")
```

### 뷰 선택 및 가중치 조정

```r
# frequency 뷰만 사용
res <- patient_dimensionality_reduction(
  seurat_obj = seurat_obj,
  markers_df = markers_df,
  include_signatures = FALSE,
  include_latent = FALSE,
  frequency_weight = 2.0  # 빈도 정보에 더 높은 가중치
)

# signatures와 latent만 사용 (frequency 제외)
res <- patient_dimensionality_reduction(
  seurat_obj = seurat_obj,
  markers_df = markers_df,
  include_frequency = FALSE
)
```

### 메타데이터 보강

```r
# 원본 Seurat 객체의 메타데이터에서 추가 변수 가져오기
aug <- augment_plot_df(
  plot_df = res$plot_df,
  metadata = is5@meta.data,
  id = "hos_no",
  add = c("g3", "nih_change", "GEM")
)

# 보강된 데이터로 플로팅
plot_patient_umap(aug$data, color_by = "g3")

# 진단 리포트 확인
aug$report  # missing IDs, NA counts, duplicates 등
```

## Unbiased Clustering 및 분리 평가

이 파이프라인은 **unsupervised** 방식으로 환자 수준의 임베딩을 생성합니다. 즉, 임상 변수(`g3`, `nih_change` 등)를 사용하지 않고 순수하게 feature 기반으로 차원 축소를 수행합니다.

### 분리 평가 방법

생성된 임베딩이 임상 변수와 얼마나 잘 연관되어 있는지 평가하려면:

1. **PERMANOVA (Permutational Multivariate Analysis of Variance)**
   - `vegan::adonis2()` 사용
   - `g3` 그룹 간 거리 차이의 통계적 유의성 검정
   - 예시:
     ```r
     library(vegan)
     dist_mat <- dist(res$embedding)
     adonis2(dist_mat ~ g3, data = res$plot_df, permutations = 9999)
     ```

2. **LISI (Local Inverse Simpson's Index)**
   - 각 환자 주변의 이웃들 중 `g3` 그룹 다양성 측정
   - 낮은 LISI = 좋은 분리 (이웃들이 같은 그룹)
   - 높은 LISI = 나쁜 분리 (이웃들이 여러 그룹에 섞임)
   - 예시:
     ```r
     library(lisi)
     lisi_scores <- compute_lisi(res$embedding, res$plot_df, c("g3"))
     # 낮은 LISI 값이 좋은 분리를 의미
     ```

### 해석

- **PERMANOVA p-value < 0.05**: `g3` 그룹 간 거리 차이가 통계적으로 유의함 → feature/embedding이 `g3`를 어느 정도 분리함
- **PERMANOVA p-value ≥ 0.05**: feature/embedding만으로는 `g3`를 분리하기 어려움 → 추가 feature나 다른 접근 필요
- **낮은 LISI**: 임베딩 공간에서 `g3` 그룹이 잘 분리됨
- **높은 LISI**: 임베딩 공간에서 `g3` 그룹이 섞여 있음

## 설계 원칙

1. **다중 뷰 통합**: 빈도/시그니처/잠재축을 독립적으로 스케일링 후 가중 결합
2. **조성 데이터 처리**: 빈도는 CLR 변환으로 처리하여 조성 특성 보존
3. **안전한 차원 축소**: 데이터 크기에 따라 PCA 백엔드 자동 선택
4. **검증 및 진단**: 메타데이터 보강 시 충돌/NA 자동 체크

## 참고

- 원본 구현: `st/KDW_251110.Rmd` (lines 4303-4428)
- 관련 함수: `myR/R/patient_dim_reduction.R`
- 테스트 스크립트: `scripts/patient_dim_reduction_downsample.R`

