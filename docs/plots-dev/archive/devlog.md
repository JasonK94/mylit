# Plot Functions Development Log

## 2024-11-12: Initial Standardization

### 목표
scRNAseq/GeoMx 데이터 플롯 함수들의 통폐합 및 표준화

### 완료된 작업

#### 1. Core Utilities 생성
- **`utils_data.R`**: 데이터 변환 및 검증 함수
  - `.sobj_to_df()`: Seurat 객체를 dataframe으로 변환 (metadata join 지원)
  - `.validate_df()`: 데이터프레임 무결성 검증
  - `.prepare_plot_data()`: Seurat/df 입력을 통합 처리

- **`utils_aggregation.R`**: Cell aggregation 로직
  - `.aggregate_cells()`: 그룹 변수별 셀 집계 (다단계 지원)
  - `.prepare_data_with_aggregation()`: 유연한 aggregation 파이프라인
  - `.reshape_for_heatmap()`: Heatmap용 데이터 reshape

#### 2. 표준화된 Plot 함수 생성
- **`plots_scatter.R`**
  - `plot_scatter()`: Expression vs covariate scatter plot
  - Seurat/df 입력, aggregation, split.by, color_by, fitted_line 지원

- **`plots_volcano.R`**
  - `plot_volcano()`: Effect size vs p-value volcano plot
  - 유연한 column naming, filtering, labeling 지원

- **`plots_heatmap.R`**
  - `plot_heatmap_genes()`: Gene expression heatmap
  - `plot_heatmap_genesets()`: Gene set (metadata) heatmap
  - Z-score normalization, clustering 지원

#### 3. 설계 원칙
1. **입력 데이터 유연성**: Seurat 객체와 data.frame 모두 지원
2. **Aggregation 자유도**: `group.by`, `split.by`, `sample_col`로 다단계 aggregation
3. **Gene Set 처리**: metadata column으로 처리 (별도 aggregation 불필요)
4. **표준화**: 일관된 파라미터 네이밍 및 공통 유틸리티 사용

#### 4. Worktree 설정
- `/data/user3/git_repo/_wt/plots`를 `plots-dev` 브랜치의 worktree로 설정
- 메인 브랜치와 분리된 개발 환경

### 2024-11-12: Boxplot 통합 및 테스트

#### 완료된 작업
- **`plots_box.R`**: Boxplot 함수 통합
  - `plot_box()`: Cell-level 및 aggregated boxplot 통합
  - `mybox_cell`, `mybox_pb` 기능 통합
  - Seurat/df 지원, aggregation 지원

- **배경 문제 수정**: 모든 plot 함수에 `panel.background`, `plot.background` 명시적 설정
  - `plots_scatter.R`: theme_bw() + white background
  - `plots_volcano.R`: theme_bw() + white background  
  - `plots_heatmap.R`: theme_minimal() + white background

- **테스트 스크립트**: `test_plots.R` 생성
  - TXNIP 기준 테스트
  - 3계층 grouping 지원 (cell → patient → group)

- **함수 문제점 분석**: `function_issues.md` 작성
  - `upset_gene_lists`, `vln_p`, `cmb`, `acmb`, `cml`, `cdf`, `cdf_multi` 분석
  - 공통 문제점 및 개선 제안 정리

### 다음 단계
- [ ] 실제 데이터로 테스트 실행
- [ ] 테스트 결과 확인 및 버그 수정
- [ ] 기존 함수들과의 호환성 확인

## 2024-12-XX: Enhanced Plot Customization Options

### Completed Features

1. **NA Value Handling**
   - Added `remove_na = FALSE` option for both scatter and heatmap plots
   - NA values are now properly handled with `na.value = "grey90"` in heatmap scales
   - NA labels are displayed in scatter plots when `label = TRUE`

2. **Normalization Control**
   - Fixed `normalize_transpose` logic: `FALSE` = column normalization (sample/group level), `TRUE` = row normalization (gene level)
   - Default behavior: `normalize_by = "row"` for gene-level normalization (better for group comparisons)

3. **Scatter Plot Enhancements**
   - `each_fit = TRUE`: Fit separate regression lines for each `color_by` group
   - `show_stats = TRUE`: Display regression statistics (y, p-value) for each fit
   - `stats_in_legend = FALSE`: Control whether stats appear in legend or as annotations
   - `shape.group.by`: Add point shape grouping for additional visual distinction
   - `label = TRUE`: Add group labels with automatic coloring and legend support
   - Labels automatically use `color_by` or `split.by` for coloring when available
   - NA values in labels are properly displayed as "NA"

4. **Heatmap Enhancements**
   - `show_group_separator = TRUE`: Add vertical dashed lines between groups
   - Enhanced facet separation: `strip.background` with black border and reduced `panel.spacing` to clearly separate g3=1 and g3=2 panels
   - Column normalization option via `normalize_transpose = FALSE`

5. **Group-level Scatter Plot**
   - Default `label = TRUE` for group-level plots
   - Automatic coloring and legend for labels
   - `each_fit = TRUE` support for separate fits per group

### Testing
- All 6 plot types successfully generated:
  1. Cluster-level scatter (with each_fit)
  2. Cluster-level heatmap (column normalization)
  3. Patient-level scatter (with each_fit)
  4. Patient-level heatmap (split by g3, column normalization, facet borders)
  5. Group-level scatter (with labels, coloring, legend, each_fit)
  6. Group-level heatmap (with NA removal, vertical separators)

