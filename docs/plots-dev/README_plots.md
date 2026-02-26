# Plot Functions Development

이 디렉토리는 scRNAseq/GeoMx 데이터 플롯 함수들의 표준화 및 통합 개발을 위한 worktree입니다.

## 구조

### Core Utilities

- **`utils_data.R`**: 데이터 변환 및 검증
  - `.sobj_to_df()`: Seurat 객체를 dataframe으로 변환
  - `.validate_df()`: 데이터프레임 무결성 검증
  - `.prepare_plot_data()`: Seurat/df 입력을 통합 처리

- **`utils_aggregation.R`**: Cell aggregation 로직
  - `.aggregate_cells()`: 그룹 변수별 셀 집계
  - `.prepare_data_with_aggregation()`: 유연한 aggregation 지원
  - `.reshape_for_heatmap()`: Heatmap용 데이터 reshape

### Plot Functions

- **`plots_scatter.R`**: Scatter plot 계열 통합
  - `plot_scatter()`: Expression vs covariate scatter plot with aggregation
  
- **`plots_volcano.R`**: Volcano plot 계열 통합
  - `plot_volcano()`: Effect size vs p-value volcano plot
  
- **`plots_heatmap.R`**: Heatmap 계열 통합
  - `plot_heatmap_genes()`: Gene expression heatmap
  - `plot_heatmap_genesets()`: Gene set (metadata) heatmap
  
- **`plots_umap_density.R`**: UMAP 빈도/커널 밀도 플롯
  - `plot_umap_density()`: 2D KDE 기반 밀도, zebra/color 스타일 및 ident 차이 맵
  
- **`plots_box.R`**: Boxplot 계열 (예정 - 기존 mybox 함수들)

## 설계 원칙

### 1. 입력 데이터 유연성
- Seurat 객체와 data.frame 모두 지원
- `metadata_df`로 추가 메타데이터 join 가능
- `match_id`로 다양한 단위(cell/patient) 매칭 지원

### 2. Aggregation 자유도
- `group.by`, `split.by`: Seurat 스타일의 그룹핑
- `sample_col`: Patient/sample 단위 aggregation
- `aggregate_by`: 다단계 aggregation 지원
  - 예: `cell -> patient -> group`

### 3. Gene Set 처리
- Gene set은 metadata로 처리 (별도 aggregation 불필요)
- Seurat의 `FetchData()`로 metadata와 feature 동시 추출 가능

### 4. 표준화
- 모든 함수는 `.prepare_plot_data()` 또는 `.prepare_data_with_aggregation()` 사용
- 일관된 파라미터 네이밍 (`group.by`, `split.by`, `sample_col`)
- `utils_plots.R`의 유틸리티 함수 활용

## 사용 예시

### Scatter Plot
```r
# Basic scatter: expression vs covariate
p <- plot_scatter(sobj, feature = "CD3D", x_var = "nih_change", 
                  group.by = "sample", aggregate = TRUE)

# With split.by and smoothing
p <- plot_scatter(sobj, feature = "CD3D", x_var = "nih_change",
                  split.by = "condition", fitted_line = "linear")
```

### Volcano Plot
```r
# From statistical results
p <- plot_volcano(lmm_results, x_col = "estimate", y_col = "p_value",
                  effect_threshold = 0.5, p_threshold = 0.05, label_top = 20)

# With filtering
p <- plot_volcano(lmm_results, filter_pattern = "drug.*timepoint")
```

### Heatmap
```r
# Gene expression heatmap
p <- plot_heatmap_genes(sobj, features = c("CD3D", "CD4", "CD8A"),
                        group.by = "seurat_clusters", normalize = TRUE)

# Gene set heatmap (metadata columns)
p <- plot_heatmap_genesets(sobj, features = c("T_cell_score", "B_cell_score"),
                           group.by = "cluster")
```

### Data Preparation
```r
# Cell-level data
plot_df <- .prepare_plot_data(sobj, features = "CD3D", 
                               group.by = "cell_type", split.by = "condition")

# Patient-level aggregation
plot_df <- .prepare_data_with_aggregation(sobj, features = "CD3D",
                                           sample_col = "patient", aggregate = TRUE)
```

### UMAP Density
```r
# 기본 UMAP 밀도 (zebra)
plot_umap_density(sobj, reduction = "umap.scvi")

# condition별 facet + color heatmap
plot_umap_density(
  sobj,
  reduction = "umap",
  split.by = "condition",
  style = "color",
  palette = "turbo"
)

# 두 그룹 간 signed density difference + 클러스터 overlay
plot_umap_density(
  sobj,
  ident.by = "condition",
  ident.1 = "responder",
  ident.2 = "non_responder",
  group.by = "anno3",
  overlay = TRUE,
  overlay_type = "both",          # contour + grid
  overlay_label = TRUE,
  overlay_grid_min_density = 0.1,
  style = "color"
)
```

## 개발 상태

- [x] `utils_data.R` - 데이터 변환/검증
- [x] `utils_aggregation.R` - Aggregation 로직
- [x] `plots_scatter.R` - Scatter plot 통합
- [x] `plots_volcano.R` - Volcano plot 통합
- [x] `plots_heatmap.R` - Heatmap 통합
- [x] `plots_umap_density.R` - UMAP KDE/빈도 플롯
- [ ] `plots_box.R` - Boxplot 통합 (예정)

## 참고 함수

기존 코드베이스의 다음 함수들을 참고:
- `scatter_smooth_colored`, `scatter_smooth_colored3` (test.R)
- `plot_volcano` (test.R, test_claude.R)
- `myhm_genes4`, `myhm_genesets4` (plots.R)
- `mybox_cell`, `mybox_pb` (plots.R)

