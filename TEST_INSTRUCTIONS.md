# 테스트 및 분석 실행 방법

## 환경 설정

### 1. R 세션 시작
```bash
# alias st를 치면 들어가는 디렉터리에서 R 시작
cd /home/user3/GJC_KDW_250721
R
```

### 2. R 세션에서 실행
```r
# start.R 자동 실행 (Rprofile.R에서 설정되어 있을 수 있음)
# 또는 수동으로 실행
source("st/start.R")

# 함수 소스 로드
source("/home/user3/data_user3/git_repo/_wt/main2/myR/R/test_analysis.R")

# 또는 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/main2/myR")
```

## 테스트 실행

### 1. 데이터 확인
```r
# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# 메타데이터 확인
colnames(sobj@meta.data)
table(sobj@meta.data$type, useNA = "ifany")
table(sobj@meta.data$g3, useNA = "ifany")
table(sobj@meta.data$seurat_clusters, useNA = "ifany")
```

### 2. runMUSCAT2_v1 테스트
```r
# 처음 3개 클러스터만 테스트
res_muscat2 <- runMUSCAT2_v1(
  sobj = sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  contrast = "IS - SAH",  # 또는 사용 가능한 그룹으로 변경
  method = "edgeR",
  remove_na_groups = TRUE,
  keep_clusters = c("0", "1", "2")  # 처음 3개 클러스터만
)

# 결과 확인
head(res_muscat2)
nrow(res_muscat2)
colnames(res_muscat2)

# 결과 저장
qs::qsave(res_muscat2, "/data/user3/sobj/test_muscat2_v1_result.qs")
```

### 3. runNEBULA2_v1 테스트
```r
# 작은 서브셋으로 테스트 (처음 1000개 유전자만)
sobj_sub <- sobj[1:min(1000, nrow(sobj)), ]

res_nebula2 <- runNEBULA2_v1(
  sobj = sobj_sub,
  layer = "counts",
  fixed_effects = "g3",
  covar_effects = NULL,
  patient_col = "hos_no",
  offset = "nCount_RNA",
  min_count = 10,
  remove_na_cells = TRUE
)

# 결과 확인
str(res_nebula2, max.level = 2)

# 결과 저장
qs::qsave(res_nebula2, "/data/user3/sobj/test_nebula2_v1_result.qs")
```

### 4. runNEBULA2_v1_with_pseudobulk 테스트
```r
# 처음 3개 클러스터만 테스트
res_nebula2_pb <- runNEBULA2_v1_with_pseudobulk(
  sobj = sobj,
  layer = "counts",
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  fixed_effects = "g3",
  covar_effects = NULL,
  patient_col = "hos_no",
  offset_method = "sum",
  min_count = 10,
  min_cells_per_pb = 3,
  remove_na_cells = TRUE,
  keep_clusters = c("0", "1", "2")  # 처음 3개 클러스터만
)

# 결과 확인
str(res_nebula2_pb, max.level = 2)
nrow(res_nebula2_pb$pseudobulk_meta)
nrow(res_nebula2_pb$pseudobulk_counts)

# 결과 저장
qs::qsave(res_nebula2_pb, "/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs")
```

## 전체 테스트 스크립트

### test_functions.R 실행
```r
# test_functions.R을 열어서 테스트 플래그를 TRUE로 변경
# test_muscat2 <- TRUE
# test_nebula2 <- TRUE
# test_nebula2_pb <- TRUE

# 스크립트 실행
source("/home/user3/data_user3/git_repo/_wt/main2/test_functions.R")
```

## 결과 확인

### 결과 파일 위치
- `test_muscat2_v1_result.qs`: runMUSCAT2_v1 결과
- `test_nebula2_v1_result.qs`: runNEBULA2_v1 결과
- `test_nebula2_v1_with_pseudobulk_result.qs`: runNEBULA2_v1_with_pseudobulk 결과

### 결과 로드
```r
# 결과 로드
res_muscat2 <- qs::qread("/data/user3/sobj/test_muscat2_v1_result.qs")
res_nebula2 <- qs::qread("/data/user3/sobj/test_nebula2_v1_result.qs")
res_nebula2_pb <- qs::qread("/data/user3/sobj/test_nebula2_v1_with_pseudobulk_result.qs")

# 결과 확인
head(res_muscat2)
str(res_nebula2, max.level = 2)
str(res_nebula2_pb, max.level = 2)
```

## 주의사항

1. **패키지 설치**: 필요한 패키지가 설치되어 있어야 합니다
   - Seurat
   - muscat
   - nebula
   - SingleCellExperiment
   - SummarizedExperiment
   - 기타

2. **메모리**: 큰 데이터셋의 경우 메모리 사용량이 많을 수 있습니다

3. **실행 시간**: 
   - runMUSCAT2_v1: 수 분 ~ 수십 분
   - runNEBULA2_v1: 수십 분 ~ 수 시간
   - runNEBULA2_v1_with_pseudobulk: 수십 분 ~ 수 시간

4. **g3 결측치**: g3 변수에 결측치가 있을 수 있으므로 `remove_na_cells=TRUE` 사용 권장

## 문제 해결

### 패키지 오류
```r
# 패키지 설치 확인
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("muscat", quietly = TRUE)) {
  BiocManager::install("muscat")
}
if (!requireNamespace("nebula", quietly = TRUE)) {
  install.packages("nebula")
}
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

# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
```

## 다음 단계

1. 테스트 성공 후 전체 데이터로 확장
2. 결과 분석 및 시각화
3. 다른 함수와의 통합 테스트
4. 성능 최적화

