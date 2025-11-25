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

