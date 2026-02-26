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
# start.R 자동 실행 (Rprofile.R에서 설정되어 있을 수 있음; "/home/user3/GJC_KDW_250721" 에서는 자동 정상 실행됨)
# 또는 수동으로 실행
# source("st/start.R")
#
패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
# 함수 소스 로드
source("/home/user3/data_user3/git_repo/_wt/<your_worktree>/myR/scripts/test_scripts.R") #example.

```

## 테스트 실행

### 1. 데이터 확인
```r
# 데이터 로드
library(qs)
is5s <- qs::qread("/data/user3/sobj/IS_scvi_251107_ds2500.qs") #for test; 2500 cells downsampled
is5 = qs::qread("/data/user3/sobj/IS_scvi_251107.qs") #full dataset; seurat object.

# 메타데이터 확인
colnames(is5s@meta.data)
table(is5s@meta.data$hos_no, useNA = "ifany") #sample ids
meta_clinical=is5@meta.data%>%distinct(hos_no,.keep_all=T)

table(is5s@meta.data$hos_no, useNA = "ifany") #sample ids

table(is5s@meta.data$sex, useNA = "ifany") #can be used as biological cofactor
table(is5s@meta.data$g3, useNA = "ifany") #target group variable; values are "NA", 1, 2, so basically charcater but can be recognized as numeric after removing "NA"; thus caution needed.
table(is5s@meta.data$anno3.scvi, useNA = "ifany") #Cluster_annotated (scvi integrated)


table(meta_clinical$hos_no, meta_clinical$GEM) #patients are in a specific GEM. be caustious or use nested formula.
table(meta_clinical$hos_no, meta_clinical$g3) #patients are in a specific group.

table(meta_clinical$GEM, meta_clinical$g3) #some are perfectly separated; be cautious when setting formula

sample_key="hos_no"
batch_key="GEM"
group_key="g3"
covar_key="sex"
cluster_key="anno3.scvi"
```

### 2. runMUSCAT 테스트
```r
# 처음 3개 클러스터만 테스트
res_muscat2 <- runMUSCAT(
  sobj = is5s,
  cluster_id = cluster_key,
  sample_id = sample_key,
  group_id = group_key,
  contrast = "2 - 1",  # group 2: g3==2, bad prognosis.
  method = "edgeR",
  remove_na_groups = TRUE,
  keep_clusters = unique(is5s$anno3.scvi)c(1:3)  # 일부 클러스터만.
)

# 결과 확인
head(res_muscat2)
nrow(res_muscat2)
colnames(res_muscat2)

# 결과 저장
qs::qsave(res_muscat2, "/data/user3/sobj/test_muscat2_v1_result.qs")
```

### 3. 전체 클러스터 테스트
```r
# 전체 클러스터 테스트
res_muscat2 <- runMUSCAT(
  sobj = is5s,
  cluster_id = cluster_key,
  sample_id = sample_key,
  group_id = group_key,
  contrast = "2 - 1",  # group 2: g3==2, bad prognosis.
  method = "edgeR",
  remove_na_groups = TRUE
)

# 결과 확인
head(res_muscat2)
nrow(res_muscat2)
colnames(res_muscat2)

# 결과 저장
qs::qsave(res_muscat2, "/data/user3/sobj/test_muscat2_v1_result.qs")
```

### 4. 원본데이터 테스트
```r
# 원본 데이터로 테스트
res_muscat3 <- runMUSCAT(
  sobj = is5,
  cluster_id = cluster_key,
  sample_id = sample_key,
  group_id = group_key,
  contrast = "2 - 1",  # group 2: g3==2, bad prognosis.
  method = "edgeR",
  remove_na_groups = TRUE,
  keep_clusters = unique(is5$anno3.scvi)c(1:3)  # 일부 클러스터만.
)

# 결과 확인
head(res_muscat3)
nrow(res_muscat3)
colnames(res_muscat3)

# 결과 저장
qs::qsave(res_muscat3, "/data/user3/sobj/test_muscat2_v1_result.qs")
```

## 전체 테스트 스크립트

### test_functions.R 실행
```r
# test_functions.R을 열어서 테스트 플래그를 TRUE로 변경
# test_muscat2 <- TRUE
# test_nebula2 <- TRUE
# test_nebula2_pb <- TRUE

# 스크립트 실행
source("/home/user3/data_user3/git_repo/mylit/scripts/test_functions.R")
```

## 결과 확인

### 결과 파일 위치
- `test_muscat2_v1_result.qs`: runMUSCAT 결과
- `test_nebula2_v1_result.qs`: runNEBULA 결과
- `test_nebula2_v1_with_pseudobulk_result.qs`: runNEBULA (pseudobulk mode) 결과

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
   - runMUSCAT: 수 분 ~ 수십 분
   - runNEBULA: 수십 분 ~ 수 시간
   - runNEBULA (pseudobulk mode): 수십 분 ~ 수 시간

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
source("/home/user3/data_user3/git_repo/_wt/analysis/myR/R/test_analysis.R")

# 또는 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/analysis/myR")
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

