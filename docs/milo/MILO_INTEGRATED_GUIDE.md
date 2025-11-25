# MiloR Differential Abundance Analysis

## Overview

MiloR is a tool for analyzing differential abundance (DA) of cell types or states in single-cell RNA sequencing data. This package uses a neighborhood-based approach to detect abundance changes in spatially adjacent cell groups.

## Key Concepts

### Neighborhood
- Local cell groups consisting of k-nearest neighbor cells around each cell
- Created with `makeNhoods()` function, each neighborhood can contain multiple cells
- Neighborhoods are connected via graph structure, so they are not independent

### Block Method
Methods for creating blocks to account for non-independence between neighborhoods:

- **`"sample"`**: Groups neighborhoods belonging to the same block based on Block ID
  - **Priority 1**: If `block_var` parameter is provided, use directly from `colData(milo)` (recommended)
    - Examples: `patient_var`, `batch_var`, `GEM`, `set`, etc.
  - **Priority 2**: Extract from cell name (fallback)
    - `S\d+$` pattern (e.g., "S1", "S2") or first 8 characters, or part before first separator
    - **Note**: GEM well suffix (`-1`, `-2`, etc.) is not removed (to prevent duplicates)
  - Assigns the block containing the most cells to each neighborhood
  - Preserves block-level structure to ensure validity of permutation test

- **`"community"`**: Uses community detection (Louvain algorithm) of nhoodGraph to group graph-connected neighborhoods into blocks
  - Used when sample information is unavailable or cannot be extracted
  - Creates pseudo-blocks based on graph connection structure

- **`"none"`**: Processes entire dataset as one block without blocking
  - Simplest method but may not fully resolve non-independence issues

## Functions

### `run_milo_pipeline`

Integrated function that executes the entire MiloR differential abundance analysis pipeline.

#### Features
- Converts Seurat object to Milo object
- Builds kNN graph and creates neighborhoods
- Calculates distances between neighborhoods
- Performs differential abundance testing (GLM-based)
- Visualizes results (UMAP, nhood graph, beeswarm plot)
- Supports caching of intermediate results

#### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `seurat_obj` | Seurat object (if loaded in memory) | `NULL` |
| `seurat_qs_path` | Path to Seurat object `.qs` file | `NULL` |
| `patient_var` | Patient/sample identifier column name | Required |
| `cluster_var` | Cluster identifier column name | Required |
| `target_var` | Group variable for comparison (e.g., treatment) | Required |
| `batch_var` | Batch effect variable | Required |
| `graph_reduction` | Dimensionality reduction for graph construction | `"integrated.scvi"` |
| `layout_reduction` | Dimensionality reduction for visualization | `"umap.scvi"` |
| `k` | k value for kNN graph | `30` |
| `d` | Number of dimensions for dimensionality reduction | `30` |
| `prop` | Proportion for neighborhood creation | `0.1` |
| `alpha` | Significance threshold | `0.1` |
| `target_include` | Groups to include from `target_var` values (e.g., `c("low","high")`) | `NULL` |
| `target_levels` | Vector forcing `target_var` factor order (`c("low","high")`) | `NULL` |
| `save` | Whether to save intermediate results | `TRUE` |
| `output_dir` | Directory to save results | `tempdir()/milo` |
| `prefix` | Filename prefix | `"milo"` |
| `suffix` | Filename suffix (auto-increment) | `NULL` |
| `cache_files` | Reuse only specified stages as cache (`nhoods`, `distances`, `da_milo`, `da_results`, `plots`, etc.) | `NULL` |
| `force_run` | Whether to ignore cache and re-run | `FALSE` |
| `plotting` | Whether to perform visualization | `TRUE` |
| `max_cells` | Maximum number of cells to downsample | `NULL` |

#### Return Value
- `milo`: Milo object
- `da_results`: Differential abundance test results (includes existing columns + `comparison_reference`, `comparison_test`, `comparison_label`, `enriched_in` for immediate logFC direction checking)
- `plots`: List of visualization objects (UMAP, nhood graph, beeswarm plot)

#### Usage Example

```r
# Run directly on Seurat object
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

# Load from .qs file and run
result <- run_milo_pipeline(
    seurat_qs_path = "/path/to/seurat.qs",
    patient_var = "patient_id",
    cluster_var = "seurat_clusters",
    target_var = "treatment",
    batch_var = "batch"
)
```

### `test_cluster_logfc_bias`

Function that tests logFC bias by cluster. Provides multiple statistical test methods accounting for non-independence between neighborhoods.

#### Features
- Calculates mean logFC by cluster
- Block permutation test (recommended)
- Correlation-adjusted t-test (neff correction)
- Mixed-effects model (optional, requires lme4)
- Empirical Bayes estimation (optional, requires ashr)

#### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `da_results` | Data frame of `miloR::testNhoods()` results | Required |
| `milo` | Milo object | Required |
| `cluster_col` | Cluster identifier column name | `"major_cluster"` |
| `block_method` | Block creation method (`"sample"`, `"community"`, `"none"`) | `"sample"` |
| `test_methods` | Test methods to execute | `c("permutation", "neff")` |
| `n_perm` | Number of permutation test iterations | `2000` |
| `max_nhoods` | Maximum number of neighborhoods to downsample | `NULL` |
| `block_var` | Block ID column name in `colData(milo)` (e.g., `patient_var`, `batch_var`, `GEM`, `set`) | `NULL` |
| `seed` | Random seed | `1` |
| `verbose` | Whether to output progress messages | `TRUE` |

#### Return Value
Data frame containing statistics by cluster:
- `major_cluster`: Cluster identifier
- `mean_logFC`: Mean logFC by cluster
- `n_nhoods`: Number of neighborhoods in cluster
- `p_perm`: Block permutation test p-value (optional)
- `p_neff`: Correlation-adjusted t-test p-value (optional)
- `p_lmm`: Mixed-effects model p-value (optional)
- `eb_mean`: Empirical Bayes posterior mean (optional)

#### 사용 예시

```r
# 기본 사용 (permutation + neff, block_var 사용 권장)
cluster_bias <- test_cluster_logfc_bias(
    da_results = da_results,
    milo = milo,
    cluster_col = "major_cluster",
    block_method = "sample",
    block_var = "patient_id",  # 또는 "batch_var", "GEM", "set" 등
    test_methods = c("permutation", "neff"),
    n_perm = 2000
)

# batch_var 사용 예시
cluster_bias_batch <- test_cluster_logfc_bias(
    da_results = da_results,
    milo = milo,
    block_method = "sample",
    block_var = "GEM",  # GEM well 기준으로 blocking
    test_methods = c("permutation", "neff")
)

# block_var 없이 cell name에서 추출 (fallback)
cluster_bias_auto <- test_cluster_logfc_bias(
    da_results = da_results,
    milo = milo,
    block_method = "sample",  # cell name에서 자동 추출 시도
    test_methods = c("permutation", "neff")
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

### logFC 방향 해석과 레벨 제어
- `target_include`: `target_var` 값 중 비교에 사용할 그룹만 남깁니다. 예를 들어 `target_include = c("low","high")`로 지정하면 "middle" 샘플은 countCells 이전에 제거됩니다.
- `target_levels`: factor 순서를 고정합니다. `target_levels = c("low","high")`이면 `logFC > 0` 은 `high`가 `low`보다 풍부하다는 뜻이며, `comparison_reference = "low"`, `comparison_test = "high"`가 `da_results`에 기록됩니다.
- `testNhoods()`는 기본적으로 factor의 첫 번째 레벨을 reference로, 두 번째 레벨을 contrast로 사용합니다. 레벨이 3개 이상이면 R의 더미 코딩 규칙을 따르므로, 관심 있는 두 수준만 비교하려면 `target_include` 또는 `target_levels`로 명시하는 것이 좋습니다.
- `da_results$enriched_in`은 `logFC` 부호를 기준으로 어느 그룹이 풍부한지를 텍스트로 보여주며, `milo$commands$.milo_run_da` 로그에도 동일 정보가 저장됩니다.

### 캐시 관리와 command 로그
- `cache_files = list(nhoods = "/abs/path/01_nhoods.qs", ...)`처럼 전달하면 해당 단계만 캐시를 재사용하고, 명시하지 않은 단계는 항상 재계산합니다. 부분 지정 시 `save = TRUE` && `suffix = NULL`이면 타임스탬프 기반 suffix가 자동 부여되어 새로운 `.qs` 출력이 기존 캐시를 덮어쓰지 않습니다.
- 모든 `.milo_run_da()` 호출은 `milo$commands`에 기록됩니다. 여기에는 사용한 캐시 경로, `target_include`, `target_levels`, FDR 가중치, 비교 참조/테스트 그룹 등 복원에 필요한 정보가 담깁니다.
- Plot 번들은 이제 `.rds`가 아닌 `.qs`로 저장되며, `cache_files$plots`를 지정해 재사용할 수 있습니다.

## 검정 방법 상세

### 1. Block Permutation Test (권장)

**원리**:
1. 각 클러스터별로 독립적으로 검정 수행
2. 관찰된 평균 logFC 계산: `obs = mean(logFC)`
3. Permutation 과정 (n_perm번 반복):
   - 각 **block 내에서만** logFC 값을 무작위로 섞음 (block 구조 보존)
   - Block별 평균 logFC 계산
   - Block 평균들의 전체 평균 계산 → `perm_mean`
4. p-value 계산: `(sum(abs(perm_means) >= abs(obs)) + 1) / (n_perm + 1)`
   - 관찰값보다 극단적인 permutation 결과의 비율

**수식**:
```
obs = mean(logFC) for cluster
For each permutation i:
  For each block b:
    logFC_perm[b] = sample(logFC[b])  # block 내에서만 섞기
  perm_mean[i] = mean(mean(logFC_perm[b] for all blocks))
p = (count(|perm_mean| >= |obs|) + 1) / (n_perm + 1)
```

**장점**:
- Block 구조를 완전히 보존하여 비독립성 문제를 완화
- Non-parametric 방법으로 분포 가정 불필요
- 가장 신뢰할 수 있는 방법

**단점**:
- 계산 시간이 오래 걸림 (n_perm에 비례, 기본값 2000)
- 작은 클러스터에서는 permutation 수가 제한적일 수 있음

### 2. Correlation-adjusted t-test (neff)

**원리**:
1. Neighborhood 그래프의 인접 행렬(adjacency matrix) 계산
2. 인접 행렬의 고유값(eigenvalue) 계산
3. Effective sample size (neff) 추정: `neff = sum(eigenvalues > 0.1)`
   - 고유값이 클수록 독립적인 정보를 많이 포함
   - 0.1보다 작은 고유값은 상관관계로 인한 중복 정보로 간주
4. 클러스터별 neff 배분: `neff_cluster = neff * (n_cluster / n_total)`
5. t-test 수행 (자유도 = neff_cluster - 1):
   ```
   t_stat = mean_logFC / (sd_logFC / sqrt(neff_cluster))
   p = 2 * (1 - pt(|t_stat|, df = neff_cluster - 1))
   ```

**수식**:
```
G = nhoodGraph adjacency matrix
eigenvals = eigenvalues(G)
neff = sum(eigenvals > threshold)  # threshold = 0.1
neff_cluster = neff * (n_cluster / n_total)
t = mean_logFC / (sd_logFC / sqrt(neff_cluster))
p = 2 * P(T > |t|) where T ~ t(neff_cluster - 1)
```

**장점**:
- 매우 빠른 계산 속도
- 그래프 구조를 직접 활용하여 상관관계를 고려

**단점**:
- 근사적 방법이므로 permutation보다 덜 보수적일 수 있음
- 고유값 계산이 큰 그래프에서는 계산 비용이 클 수 있음
- 작은 클러스터에서는 neff_cluster가 너무 작아질 수 있음

### 3. Mixed-Effects Model (LMM)

**원리**:
1. 각 클러스터별로 선형혼합모형 피팅
2. 현재 구현: `logFC ~ 1 + (1 | block_id)`
   - 고정 효과: 절편만 (평균 logFC 추정)
   - 랜덤 효과: `block_id` (block-level 변동성 모델링)
3. 절편의 유의성 검정 (H0: mean_logFC = 0)

**수식**:
```
logFC_ij = μ + b_i + ε_ij
where:
  μ = overall mean (fixed effect)
  b_i ~ N(0, σ²_block) (random effect for block i)
  ε_ij ~ N(0, σ²_residual) (residual error)
  
Test: H0: μ = 0
```

**주의사항**:
- **현재 구현은 `batch_var`를 보정하지 않습니다**
- Block-level random effect만 모델링
- `batch_var`를 보정하려면 모델을 `logFC ~ 1 + batch_var + (1 | block_id)` 형태로 확장해야 함

**장점**:
- Block 구조를 명시적으로 모델링
- Block-level 변동성을 분리하여 추정

**단점**:
- lme4 패키지 필요
- 작은 클러스터나 block 수가 적을 때 수렴 실패 가능
- 현재는 batch 효과를 보정하지 않음

### 4. Empirical Bayes (ashr)

**원리**:
1. p-value를 z-score로 변환:
   ```
   z = qnorm(1 - PValue/2) * sign(logFC)
   ```
2. Standard error 추정:
   ```
   se = |logFC| / max(|z|, 1e-6)
   ```
3. Empirical Bayes shrinkage 적용:
   - `ashr::ash()` 함수 사용
   - Prior distribution을 데이터에서 추정
   - Posterior distribution 계산
4. Posterior mean 계산 (effect size 추정):
   ```
   eb_mean = mean(posterior_mean for all neighborhoods in cluster)
   ```

**수식**:
```
z_i = qnorm(1 - PValue_i/2) * sign(logFC_i)
se_i = |logFC_i| / max(|z_i|, 1e-6)

# Empirical Bayes
posterior_i = ash(logFC_i, se_i, method="fdr")
eb_mean = mean(posterior_mean_i for all i in cluster)
```

**장점**:
- 작은 effect size도 안정적으로 추정 가능
- Shrinkage를 통해 과도한 추정을 방지
- Multiple testing 문제를 완화

**단점**:
- ashr 패키지 필요
- 검정(test)보다는 effect size 추정에 유용
- p-value가 아닌 posterior mean을 제공

**용도**:
- Effect size의 방향성과 크기를 추정할 때 유용
- 검정 결과를 보완하는 정보로 활용

## Warnings

### 1. Milo Object Structure
- **Important**: `milo@Nhoods` (uppercase) is the correct slot name. `milo@nhoods` (lowercase) does not exist.
- `plotNhoodGraphDA()` can only be used after calling `buildNhoodGraph()`.

### 2. Plotting Functions
- `plotDAbeeswarm()` is sensitive to `alpha` and `SpatialFDR` values. If all values are non-significant, errors may occur. Change `beeswarm_metric` to `"PValue"` or `"logFC_percentile"`, or adjust `beeswarm_alpha`.

### 3. R Environment
- R must be run from the `st/` directory so that `st/start.R` is automatically sourced.
- This ensures the `renv` environment is correctly loaded and required packages are available.

### 4. Caching
- When `save=TRUE` (default), intermediate results are saved in `.qs` format.
- When `force_run=FALSE`, existing cache files are automatically reused to save time.
- Cache files are created for each stage (`nhoods`, `distances`, `testing`).

## Workflow Example

```r
# 1. Run Milo pipeline
milo_result <- run_milo_pipeline(
    seurat_qs_path = "/path/to/seurat.qs",
    patient_var = "patient_id",
    cluster_var = "seurat_clusters",
    target_var = "treatment",
    batch_var = "batch",
    save = TRUE,
    output_dir = "/path/to/output"
)

# 2. Test cluster bias
cluster_bias <- test_cluster_logfc_bias(
    da_results = milo_result$da_results,
    milo = milo_result$milo,
    cluster_col = "major_cluster",
    block_method = "sample",
    test_methods = c("permutation", "neff"),
    n_perm = 2000
)

# 3. Check results
print(cluster_bias)
# Filter significant clusters
significant_clusters <- cluster_bias %>%
    filter(p_perm < 0.05 | p_neff < 0.05)
```

## References

- MiloR official documentation: https://www.bioconductor.org/packages/release/bioc/html/miloR.html
- Original paper: Dann et al. (2022) Nature Biotechnology
- For major issues and solutions during development, refer to `myR/docs/DEVLOG_Korean.md` and `myR/docs/context_Korean.md`.

