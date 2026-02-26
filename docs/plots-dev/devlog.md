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

## 2025-12-18: Bar Plots (cmb/acmb) & Boxplot Enhancements

### 완료된 작업

#### 1. `cmb` (Proportional) & `acmb` (Absolute Count) 함수 대규모 업데이트
- **다중 정렬 시스템 (Sequential Sorting)**: `sort.by`에 벡터를 입력받아 여러 기준으로 순차적 정렬 지원.
  - 우선순위: Identity 클러스터 빈도 > Metadata 평균 > Feature(Gene) 발현량 평균.
- **그룹 분할 및 자동 구분선 (Split by Group)**: 
  - `split.by` 변수를 기준으로 1차 그룹화 후 내부 정렬 수행.
  - `split.by` 값이 변하는 지점에 자동적으로 수직 점선(`vlines`)을 삽입하여 시각적 가독성 향상.
- **샘플 라벨 컬러링 (label.by)**:
  - x축 샘플 이름 영역에 메타데이터 기반의 색상 박스(rect) 표시.
  - `label.numeric`: 연속형(Gradient)과 범주형(Discrete)을 명시적으로 구분하거나 자동 판별하는 로직 구현.
- **정렬 보조 그래프 (sort.graph)**:
  - `sort.graph.overlay = TRUE`: 메인 Bar plot 위에 Secondary Y-axis를 사용하여 선 그래프를 중첩.
  - `sort.graph.overlay = FALSE`: `patchwork`를 이용해 메인 그래프 위쪽에 독립적인 선 그래프 배치 (좌측 Y축 사용).
  - `sort.graph.iqr`: 정렬 기준이 되는 데이터의 IQR(1/3 사분위수) 영역을 음영으로 표시 지원.
- **스태킹 순서 제어 (cluster_order)**: 누적 바 내부의 클러스터 쌓이는 순서를 사용자 정의 가능하도록 개선.
- **상세 정보 표시 (legend_detail)**: identity, group.by, sort.by 등 사용된 주요 파라미터를 플롯 하단 캡션에 자동 생성.

#### 2. `plot_box` (Boxplot) 시각화 개선
- **통계 라벨 위치 최적화**: `ggpubr`의 p-value 라벨이 클러스터 이름과 겹치지 않도록 `vjust`, `bracket.nudge.y` 등 파라미터 조정.
- **독립적 Y축 스케일 (Free Y)**: `facet_wrap` 시 `scales = "free_y"`를 적용하여 세포 수가 적은 희귀 클러스터의 분포가 뭉개지지 않고 잘 보이도록 수정.

#### 3. 내부 로직 및 안정성 강화
- **`.calculate_sort_values` 헬퍼 함수**: 다중 정렬 기준 병합, IQR 계산, Seurat 객체 데이터 추출 로직을 하나로 통합 및 고도화.
- **레이아웃 제어**: `coord_cartesian(clip = "off")`와 동적 여백 설정을 통해 라벨 박스나 보조 그래프가 플롯 영역 밖에서 잘리는 문제 해결.

### 다음 단계
- [ ] 정렬 및 보조 그래프 로직의 엣지 케이스(데이터 부족 등) 테스트
- [ ] `cmb`/`acmb` 함수를 `plots_bar.R`로 분리 및 표준화 모듈에 포함
