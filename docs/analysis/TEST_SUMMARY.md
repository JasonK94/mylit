# 테스트 및 분석 진행 요약

## 완료된 작업

### 1. 함수 개발
- ✅ `runMUSCAT`: MUSCAT 분석 함수 (결측치 처리 개선)
- ✅ `runNEBULA`: NEBULA 분석 함수 (결측치 처리 개선)
- ✅ `runNEBULA (pseudobulk mode)`: Pseudobulk와 NEBULA 결합 함수

### 2. 테스트 스크립트 작성
- ✅ `test_interactive.R`: 대화형 테스트 스크립트 (권장)
- ✅ `test_simple.R`: 기본 데이터 확인 스크립트
- ✅ `test_step1_data_check.R`: 데이터 로드 및 검증 스크립트
- ✅ `test_run_muscat2.R`: runMUSCAT 전용 테스트 스크립트
- ✅ `test_functions.R`: 전체 함수 테스트 스크립트
- ✅ `test_functions_quick.R`: 빠른 함수 검증 스크립트

### 3. 문서 작성
- ✅ `QUICK_START.md`: 빠른 시작 가이드
- ✅ `TEST_INSTRUCTIONS.md`: 상세한 테스트 지침
- ✅ `DEVELOPMENT_SUMMARY.md`: 개발 요약 및 함수 설명
- ✅ `TEST_SUMMARY.md`: 테스트 요약 (본 문서)

### 4. 실행 스크립트
- ✅ `run_test_muscat2.sh`: 배치 테스트 실행 스크립트

## 테스트 실행 방법

### 방법 1: 대화형 R 세션 (권장)

```bash
# 1. 작업 디렉터리로 이동
cd /home/user3/GJC_KDW_250721

# 2. R 세션 시작
R

# 3. R 세션에서 실행
source("/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_interactive.R")
```

### 방법 2: 수동 실행

```r
# 1. 환경 설정
source("st/start.R")
source("/home/user3/data_user3/git_repo/_wt/analysis/myR/R/test_analysis.R")

# 2. 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# 3. runMUSCAT 테스트
res_muscat2 <- runMUSCAT(
  sobj = sobj,
  cluster_id = "seurat_clusters",
  sample_id = "hos_no",
  group_id = "type",
  contrast = "IS - SAH",
  method = "edgeR",
  remove_na_groups = TRUE,
  keep_clusters = c("0", "1", "2")
)

# 4. 결과 확인 및 저장
head(res_muscat2)
qs::qsave(res_muscat2, "/data/user3/sobj/test_muscat2_v1_result.qs")
```

## 테스트 순서

### 1단계: 데이터 확인
```r
source("/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_simple.R")
```

### 2단계: runMUSCAT 테스트
```r
source("/home/user3/data_user3/git_repo/_wt/analysis/scripts/analysis/test_interactive.R")
```

### 3단계: runNEBULA 테스트 (선택)
```r
# 작은 서브셋으로 테스트
sobj_sub <- sobj[1:min(1000, nrow(sobj)), ]

res_nebula2 <- runNEBULA(
  sobj = sobj_sub,
  layer = "counts",
  fixed_effects = "g3",
  patient_col = "hos_no",
  offset = "nCount_RNA",
  remove_na_cells = TRUE
)

qs::qsave(res_nebula2, "/data/user3/sobj/test_nebula2_v1_result.qs")
```

### 4단계: runNEBULA (pseudobulk mode) 테스트 (선택)
```r
res_nebula2_pb <- runNEBULA (pseudobulk mode)(
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
nrow(res_muscat2)
colnames(res_muscat2)
```

## 예상 실행 시간

- **runMUSCAT**: 몇 분 ~ 수십 분 (클러스터 수와 데이터 크기에 따라)
- **runNEBULA**: 수십 분 ~ 수 시간 (유전자 수와 데이터 크기에 따라)
- **runNEBULA (pseudobulk mode)**: 수십 분 ~ 수 시간 (클러스터 수와 데이터 크기에 따라)

## 주요 파일

### 함수 정의
- `myR/R/test_analysis.R`: 함수 정의 (1,208줄)

### 테스트 스크립트
- `test_interactive.R`: 대화형 테스트 (권장)
- `test_simple.R`: 기본 데이터 확인
- `test_run_muscat2.R`: runMUSCAT 전용 테스트
- `test_functions.R`: 전체 함수 테스트

### 문서
- `QUICK_START.md`: 빠른 시작 가이드
- `TEST_INSTRUCTIONS.md`: 상세한 테스트 지침
- `DEVELOPMENT_SUMMARY.md`: 개발 요약

## 다음 단계

1. **테스트 실행**: R 세션에서 `test_interactive.R` 실행
2. **결과 확인**: 결과 파일 확인 및 분석
3. **전체 데이터 확장**: 테스트 성공 시 전체 데이터로 확장
4. **결과 분석**: 결과 분석 및 시각화
5. **성능 최적화**: 필요시 성능 최적화

## 문제 해결

### 패키지 오류
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
source("/home/user3/data_user3/git_repo/_wt/analysis/myR/R/test_analysis.R")
```

### 데이터 로드 오류
```r
# 데이터 파일 확인
file.exists("/data/user3/sobj/IS_scvi_251107_ds2500.qs")

# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs")
```

## Git 상태

- **브랜치**: main2
- **커밋 수**: 4개
- **작업 트리**: clean (모든 변경사항 커밋됨)
- **최신 커밋**: b2d8bfd

## 참고사항

1. **실제 R 세션에서 실행**: Rscript로는 패키지 환경 문제가 있을 수 있으므로, 실제 R 세션에서 실행하는 것을 권장합니다.
2. **패키지 설치**: 필요한 패키지가 모두 설치되어 있어야 합니다.
3. **메모리**: 큰 데이터셋의 경우 메모리 사용량이 많을 수 있습니다.
4. **실행 시간**: 분석에 시간이 오래 걸릴 수 있으므로, 작은 서브셋으로 먼저 테스트하는 것을 권장합니다.

