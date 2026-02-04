# Testing Guide

## Function Loading Test

모든 함수가 정상적으로 로드되는지 확인:

```bash
cd /data/user3/git_repo/_wt/plots
Rscript -e "
source('myR/R/utils_data.R')
source('myR/R/utils_aggregation.R')
source('myR/R/plots_scatter.R')
source('myR/R/plots_volcano.R')
source('myR/R/plots_heatmap.R')
cat('All functions loaded successfully!\n')
"
```

## Running Full Tests

### Option 1: With Seurat Object File

```bash
Rscript execute_test.R /path/to/your/sobj.rds
```

또는 RData 파일:
```bash
Rscript execute_test.R /path/to/your/sobj.RData
```

### Option 2: In R Environment

```r
# Load Seurat object
sobj <- readRDS("/path/to/sobj.rds")

# Source test script
source("/data/user3/git_repo/_wt/plots/test_plots.R")

# Run tests
test_plots(
  sobj = sobj,
  feature = c("TXNIP", "DDIT4", "UTY", "S100B", "XIST", "HLA-B", "CCL4", "HLA-C"),
  output_dir = "test_output",
  group.by = "anno3.scvi",
  sample_col = "hos_no",
  split.by = "g3"
)
```

### Option 3: Automated Search

`run_tests.R`는 다음 위치에서 Seurat 객체를 자동으로 찾습니다:
- `sobj.RData` / `sobj.rds` (현재 디렉토리)
- `../sobj.RData` / `../sobj.rds`
- `../../sobj.RData` / `../../sobj.rds`

또는 환경 변수 `sobj`가 설정되어 있으면 사용합니다.

```bash
Rscript run_tests.R
```

## Expected Output

테스트가 성공하면 다음 파일들이 생성됩니다:

```
test_output/
├── 01_scatter_cell.png      # Cell-level scatter plot
├── 02_heatmap_cell.png      # Cell-level heatmap
├── 03_scatter_patient.png   # Patient-level scatter plot
├── 04_heatmap_patient.png   # Patient-level heatmap
├── 05_scatter_group.png     # Group-level scatter plot
└── 06_heatmap_group.png     # Group-level heatmap
```

## Troubleshooting

### Missing Columns

필요한 컬럼이 없으면 `execute_test.R`이 자동으로 대체 컬럼을 찾습니다:
- `anno3.scvi` → `cluster`, `anno`, `cell_type` 등
- `hos_no` → `sample`, `patient`, `hos`, `id` 등
- `g3` → `group`, `condition`, `treatment`, `g[0-9]` 등

### Missing Features

테스트 feature가 assay에 없으면 첫 5개 유전자를 사용합니다.

### Error Messages

에러가 발생하면:
1. Stack trace 확인
2. `DEBUGGING.md` 참조
3. 함수별 에러 메시지 확인

## Manual Testing

개별 함수 테스트:

```r
# Load functions
source("myR/R/utils_data.R")
source("myR/R/utils_aggregation.R")
source("myR/R/plots_scatter.R")

# Test data conversion
df <- .sobj_to_df(sobj, features = NULL, metadata_only = TRUE)

# Test scatter plot
p <- plot_scatter(
  data = sobj,
  feature = "TXNIP",
  x_var = "nFeature_RNA",
  group.by = "anno3.scvi",
  aggregate = FALSE
)
print(p)
```

