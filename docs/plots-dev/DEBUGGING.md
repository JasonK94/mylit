# Debugging Guide

## Common Issues and Solutions

### 1. `rlang::%||%` Error
**Problem**: `rlang` 패키지가 로드되지 않아 `%||%` 연산자 사용 불가

**Solution**: 모든 함수에서 `%||%`를 로컬로 정의하도록 수정
```r
`%||%` <- function(x, y) if (is.null(x)) y else x
```

### 2. `match.arg()` with NULL
**Problem**: `match.arg()`는 NULL을 처리할 수 없음

**Solution**: NULL 체크 후 match.arg() 호출
```r
if (!is.null(fitted_line)) {
  fitted_line <- match.arg(fitted_line, choices = c("linear", "loess", "lasso"))
}
```

### 3. Feature Not Found
**Problem**: `feature[1]`이 "nih_change"인데 이것이 gene이 아니라 metadata일 수 있음

**Solution**: 
- Gene과 metadata를 분리하여 처리
- Scatter plot은 gene만 사용 (metadata는 x_var로 사용 가능)
- Heatmap은 gene + metadata 모두 사용 가능

### 4. Missing Dependencies
**Problem**: 필수 패키지가 로드되지 않음

**Solution**: `test_plots.R`에서 사전 체크 추가
```r
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "rlang", "tidyr", "viridisLite", "RColorBrewer")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
```

### 5. File Path Issues
**Problem**: 상대 경로로 인한 source() 실패

**Solution**: 절대 경로 사용 또는 `here::here()` 사용

### 6. Background Color Issues
**Problem**: 플롯 배경이 검은색

**Solution**: 명시적으로 white background 설정
```r
theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA)
)
```

### 7. `plot_umap_density()` Issues
**Problem**:
- `"Reduction 'umap.scvi' not found"` 또는 `ident.by column ... not found` 에러.
- `Unknown colour name: inferno` 에러.
- `object 'ndensity' not found` 또는 `Binned scales only support continuous data`.
- `overlay = TRUE`인데 contour가 일부 클러스터에만 그려지거나, grid가 거울 대칭으로 찍힘.

**Solution**:
- **Reduction/Ident**: Seurat 객체에 해당 리덕션/컬럼이 있는지 확인. `ident.1/2` 사용 시 `ident.by` 명시 필수.
- **Palette**: `inferno` 같은 viridis 옵션은 이제 내부적으로 자동 처리됩니다(패치 완료).
- **ggplot2 호환성**: `ndensity` 대신 `as.numeric(level)`을 사용하여 버전 간 호환성과 연속/이산형 스케일 충돌을 해결했습니다.
- **Overlay/Contour**: 그룹별로 `geom_density_2d` 레이어를 루프로 분리해 모든 클러스터가 강제로 그려지도록 수정했습니다.
- **Grid 좌표**: `MASS::kde2d` 결과의 x/y 매핑 순서를 `expand.grid`에서 올바르게 맞춰 해결했습니다.


## Debugging Steps

1. **Check function loading**:
   ```r
   exists(".sobj_to_df")
   exists("plot_scatter")
   ```

2. **Test data preparation**:
   ```r
   df <- .sobj_to_df(sobj, features = NULL, metadata_only = TRUE)
   ```

3. **Test with minimal example**:
   ```r
   p <- plot_scatter(sobj, feature = "TXNIP", x_var = "nFeature_RNA", 
                     group.by = "anno3.scvi", aggregate = FALSE)
   ```

4. **Check error messages**: 모든 tryCatch에서 stack trace 출력

## Test Scripts

- `test_plots.R`: 기본 테스트 스크립트
- `test_plots_debug.R`: 향상된 디버깅 기능 포함

