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

### 다음 단계
- [ ] Boxplot 함수 통합 (`mybox_cell`, `mybox_pb` → `plots_box.R`)
- [ ] 테스트 및 버그 수정
- [ ] 기존 함수들과의 호환성 확인

