# 빠른 시작 가이드: 테스트 및 분석 실행

## 준비 사항

1. **R 세션 시작**: alias st를 치면 들어가는 디렉터리에서 R 시작
   ```bash
   cd /home/user3/GJC_KDW_250721
   R
   ```

2. **필요한 패키지**: 다음 패키지들이 설치되어 있어야 합니다
   - Seurat
   - muscat
   - nebula
   - SingleCellExperiment
   - SummarizedExperiment
   - qs
   - dplyr
   - Matrix

## 빠른 테스트 실행

### 방법 1: 대화형 테스트 스크립트 (권장)

```r
# R 세션에서 실행
source("/home/user3/data_user3/git_repo/_wt/main2/test_interactive.R")
```

이 스크립트는:
1. 데이터 로드 및 확인
2. runMUSCAT2_v1 테스트 실행
3. 결과 저장

### 방법 2: 수동 실행

```r
# 1. 환경 설정
source("st/start.R")
source("/home/user3/data_user3/git_repo/_wt/main2/myR/R/test_analysis.R")

# 2. 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# 3. runMUSCAT2_v1 테스트
res_muscat2 <- runMUSCAT2_v1(
  sobj = sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  contrast = "IS - SAH",  # 또는 사용 가능한 그룹
  method = "edgeR",
  remove_na_groups = TRUE,
  keep_clusters = c("0", "1", "2")  # 처음 3개 클러스터만
)

# 4. 결과 확인
head(res_muscat2)
nrow(res_muscat2)

# 5. 결과 저장
qs::qsave(res_muscat2, "/data/user3/sobj/test_muscat2_v1_result.qs")
```

## 테스트 순서

### 1단계: 데이터 확인
```r
source("/home/user3/data_user3/git_repo/_wt/main2/test_simple.R")
```

### 2단계: runMUSCAT2_v1 테스트
```r
source("/home/user3/data_user3/git_repo/_wt/main2/test_interactive.R")
```

### 3단계: runNEBULA2_v1 테스트 (선택)
```r
# 작은 서브셋으로 테스트
sobj_sub <- sobj[1:min(1000, nrow(sobj)), ]

res_nebula2 <- runNEBULA2_v1(
  sobj = sobj_sub,
  layer = "counts",
  fixed_effects = "g3",
  patient_col = "hos_no",
  offset = "nCount_RNA",
  remove_na_cells = TRUE
)

# 결과 확인
str(res_nebula2, max.level = 2)

# 결과 저장
qs::qsave(res_nebula2, "/data/user3/sobj/test_nebula2_v1_result.qs")
```

### 4단계: runNEBULA2_v1_with_pseudobulk 테스트 (선택)
```r
res_nebula2_pb <- runNEBULA2_v1_with_pseudobulk(
  sobj = sobj,
  layer = "counts",
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  fixed_effects = "g3",
  patient_col = "hos_no",
  offset_method = "sum",
  remove_na_cells = TRUE,
  keep_clusters = c("0", "1", "2")
)

# 결과 확인
str(res_nebula2_pb, max.level = 2)

# 결과 저장
qs::qsave(res_nebula2_pb, "/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs")
```

## 결과 확인

### 결과 파일 위치
- `/data/user3/sobj/test_muscat2_v1_result.qs`
- `/data/user3/sobj/test_nebula2_v1_result.qs`
- `/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs`

### 결과 로드
```r
# 결과 로드
res_muscat2 <- qs::qread("/data/user3/sobj/test_muscat2_v1_result.qs")

# 결과 확인
head(res_muscat2)
summary(res_muscat2)
```

## 문제 해결

### 패키지 오류
패키지가 설치되지 않은 경우:
```r
# BiocManager 설치
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor 패키지 설치
BiocManager::install(c("muscat", "SingleCellExperiment", "SummarizedExperiment"))

# CRAN 패키지 설치
install.packages(c("Seurat", "nebula", "qs", "dplyr", "Matrix"))
```

### 함수가 로드되지 않음
```r
# 함수 소스 직접 로드
source("/home/user3/data_user3/git_repo/_wt/main2/myR/R/test_analysis.R")

# 또는 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/main2/myR")
```

### 데이터 로드 오류
```r
# 데이터 파일 확인
file.exists("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# qs 패키지 확인
requireNamespace("qs", quietly = TRUE)

# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
```

## 실행 시간 예상

- **runMUSCAT2_v1**: 몇 분 ~ 수십 분 (클러스터 수와 데이터 크기에 따라)
- **runNEBULA2_v1**: 수십 분 ~ 수 시간 (유전자 수와 데이터 크기에 따라)
- **runNEBULA2_v1_with_pseudobulk**: 수십 분 ~ 수 시간 (클러스터 수와 데이터 크기에 따라)

## 다음 단계

1. 테스트 성공 후 전체 데이터로 확장
2. 결과 분석 및 시각화
3. 다른 함수와의 통합 테스트
4. 성능 최적화

## 상세 문서

- `TEST_INSTRUCTIONS.md`: 상세한 테스트 지침
- `DEVELOPMENT_SUMMARY.md`: 개발 요약 및 함수 설명

