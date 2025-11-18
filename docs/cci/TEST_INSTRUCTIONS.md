# CCI 분석 도구 테스트 가이드

## 환경 설정

### 1. R 세션 시작
```bash
# alias st를 치면 들어가는 디렉터리에서 R 시작
cd /home/user3/GJC_KDW_250721
R
```

### 2. R 세션에서 실행
```r
# 패키지 로드
devtools::load_all("/home/user3/data_user3/git_repo/_wt/cci/myR")

# 또는 함수 소스 직접 로드
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
cci_core_worktree <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
cci_core_mainrepo <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
} else {
  stop("CCI.R not found in worktree or main repository.")
}
```

## 데이터 준비

### 1. 데이터 로드
```r
library(qs)

# 테스트 데이터 (다운샘플)
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 원본 데이터
sobj_full <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")
```

### 2. 메타데이터 확인
```r
# 메타데이터 컬럼 확인
colnames(sobj_test@meta.data)

# 클러스터 확인
table(sobj_test@meta.data$anno3.scvi, useNA = "ifany")

# 조건 변수 확인
table(sobj_test@meta.data$g3, useNA = "ifany")

# 샘플 ID 확인
table(sobj_test@meta.data$hos_no, useNA = "ifany")
```

### 3. DEG 리스트 준비
```r
# 예시 1: 기존 DEG 분석 결과 사용
# (예: muscat, edgeR, DESeq2 등으로 분석한 결과)

# 예시 2: 간단한 DEG 리스트 생성 (테스트용)
deg_df <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
  cluster = rep("Cluster1", 5),
  avg_log2FC = c(1.5, 2.0, -1.2, 0.8, -0.5),
  p_val_adj = c(0.001, 0.0001, 0.01, 0.05, 0.1)
)
```

### 4. 사전 계산된 receiver DEG 재사용 (`receiver_de_table`)
`cci_nichenet_results_*.qs` 또는 사용자가 별도로 저장한 `.qs` 파일에는 receiver DEG 테이블이 그대로 들어 있습니다. 동일한 receiver를 다시 분석할 때는 아래처럼 바로 다시 불러와 `run_nichenet_analysis()`에 전달하면 `FindMarkers()`가 재실행되지 않습니다.

```r
library(qs)

receiver_deg_qs <- "/data/user3/sobj/receiver_CD4_deg_table.qs"  # 필요 시 갱신
if (file.exists(receiver_deg_qs)) {
  receiver_deg_table <- qs::qread(receiver_deg_qs)
} else {
  stop("Receiver DEG table (.qs) not found. Run run_cci_analysis once and save the table.")
}

nichenet_rerun <- run_nichenet_analysis(
  seurat_obj = sobj_test,
  species = "human",
  sender_celltypes = c("Cluster2", "Cluster3"),
  receiver_celltype = "Cluster1",
  assay_name = "SCT",
  cluster_col = "anno3.scvi",
  receiver_DE_ident1 = "2",
  receiver_DE_ident2 = "1",
  receiver_DE_group_by = "g3",
  receiver_de_table = receiver_deg_table,
  receiver_gene_col = "gene",        # 필요 시 "symbol" 등으로 변경 가능
  receiver_logfc_col = "avg_log2FC", # 외부 DEG면 "logFC"/"avg_logFC"도 허용
  receiver_pval_col = "p_val_adj",
  p_val_adj_cutoff = 1.1,
  logfc_cutoff = 0.05,
  verbose = TRUE
)
```

- `run_cci_analysis()`는 위와 동일한 방식으로 자동 전달하므로, 테스트 시 메시지에 `Precomputed receiver DE tables...` 로그와 함께 FindMarkers 재호출이 생략되었는지 확인합니다.
- receiver DEG를 새로 생성해야 할 때는 `save_cci_intermediate()`가 만든 `.qs` 파일에서 `receiver_degs` 요소를 추출하여 재사용하는 것이 가장 빠릅니다.

## 테스트 실행 방법

> 참고: `run_cci_analysis()` Step 6 로그에는 sender/receiver 수, receiver DEG 수, 그리고 `Estimated runtime` 메시지가 포함됩니다. 테스트 시 해당 로그가 표준 출력에 나타나는지 확인하십시오.

### 방법 1: 대화형 R 세션에서 실행 (권장)

```r
# R 세션 시작
cd /home/user3/GJC_KDW_250721
R

# R 세션에서:
source("/home/user3/data_user3/git_repo/_wt/cci/scripts/cci/test_cci_interactive.R")
```

### 방법 2: 단계별 수동 실행

```r
# 1. 함수 로드
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
cci_core_worktree <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
cci_core_mainrepo <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
} else {
  stop("CCI.R not found in worktree or main repository.")
}

# 2. 데이터 로드
library(qs)
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 3. DEG 리스트 준비
clusters <- unique(sobj_test@meta.data$anno3.scvi)
clusters <- clusters[!is.na(clusters)]
receiver_cluster_test <- clusters[1]

deg_df_example <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  cluster = rep(receiver_cluster_test, 3),
  avg_log2FC = c(1.5, 2.0, -1.2),
  p_val_adj = c(0.001, 0.0001, 0.01),
  stringsAsFactors = FALSE
)

# 4. CCI 분석 실행
results <- run_cci_analysis(
  sobj = sobj_test,
  cluster_col = "anno3.scvi",
  deg_df = deg_df_example,
  receiver_cluster = receiver_cluster_test,
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  verbose = TRUE
)
```

## 테스트 실행

### Test 1: 기본 기능 테스트 (다운샘플 데이터)
```r
# 기본 파라미터로 실행
results_test <- run_cci_analysis(
  sobj = sobj_test,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Cluster1",  # 실제 클러스터 ID로 변경
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  output_dir = "/data/user3/sobj/cci_test_output"
)

# 결과 확인
names(results_test)
head(results_test$nichenet_results$ligand_activities)
head(results_test$nichenet_results$best_upstream_ligands)
```

### Test 2: Sender 클러스터 지정
```r
# 특정 sender 클러스터 지정
results_test2 <- run_cci_analysis(
  sobj = sobj_test,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Cluster1",
  sender_clusters = c("Cluster2", "Cluster3"),  # sender 지정
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human"
)
```

### Test 3: 자동 Sender 식별
```r
# sender_clusters를 NULL로 설정하여 자동 식별
results_test3 <- run_cci_analysis(
  sobj = sobj_test,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Cluster1",
  sender_clusters = NULL,  # 자동 식별
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human"
)
```

### Test 4: 전체 데이터 테스트
```r
# 원본 데이터로 테스트 (시간이 오래 걸릴 수 있음)
results_full <- run_cci_analysis(
  sobj = sobj_full,
  cluster_col = "anno3.scvi",
  deg_df = deg_df_full,  # 전체 데이터의 DEG 리스트
  receiver_cluster = "Cluster1",
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  output_dir = "/data/user3/sobj/cci_full_output"
)
```

## 결과 확인

### 1. 결과 구조 확인
```r
# 결과 리스트 구조
str(results_test, max.level = 2)

# NicheNet 결과 확인
names(results_test$nichenet_results)

# Ligand activity 확인
head(results_test$nichenet_results$ligand_activities)

# Top ligands 확인
results_test$nichenet_results$best_upstream_ligands

# Sender-receiver 매핑 확인
results_test$sender_receiver_map

# DEG 요약 확인
results_test$deg_summary
```

### 2. 시각화 확인
```r
# Ligand-target 네트워크 히트맵
results_test$nichenet_results$plot_ligand_target_network

# Ligand-receptor 네트워크 히트맵
results_test$nichenet_results$plot_ligand_receptor_network

# Circos plot (있는 경우)
if (!is.null(results_test$nichenet_results$plot_circos)) {
  plot(results_test$nichenet_results$plot_circos)
}
```

### 3. 결과 저장 확인
```r
# 저장된 파일 확인
list.files("/data/user3/sobj", pattern = "cci.*\\.qs$")

# 저장된 결과 로드
library(qs)
saved_results <- qs::qread("/data/user3/sobj/cci_analysis_results_YYYYMMDD_HHMMSS.qs")
```

## 결과 저장

### 자동 저장
각 단계마다 자동으로 `/data/user3/sobj`에 저장됩니다:
- `cci_prepared_data_{timestamp}.qs`: 데이터 준비 결과
- `cci_nichenet_results_{timestamp}.qs`: NicheNet 분석 결과
- `cci_analysis_results_{timestamp}.qs`: 최종 분석 결과

### 수동 저장
```r
# 결과 수동 저장
qs::qsave(results_test, "/data/user3/sobj/my_cci_results.qs")

# 특정 결과만 저장
qs::qsave(
  results_test$nichenet_results$ligand_activities,
  "/data/user3/sobj/ligand_activities.qs"
)
```

## 문제 해결

### 패키지 오류
```r
# 필수 패키지 설치 확인
if (!requireNamespace("nichenetr", quietly = TRUE)) {
  BiocManager::install("nichenetr")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
```

### 데이터 로드 오류
```r
# 데이터 파일 확인
file.exists("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 데이터 로드
library(qs)
sobj <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")
```

### DEG 리스트 형식 오류
```r
# 필수 컬럼 확인
required_cols <- c("gene", "cluster", "logFC", "p_val_adj")
# 또는
required_cols <- c("gene", "cluster", "avg_log2FC", "p_val_adj")

# 컬럼명 확인
colnames(deg_df)

# 컬럼명 변경 (필요시)
deg_df <- deg_df %>%
  rename(
    gene = GENE,
    cluster = CLUSTER,
    avg_log2FC = logFC,
    p_val_adj = FDR
  )
```

### NicheNet 데이터 다운로드 실패
```r
# 수동으로 데이터 디렉터리 지정
results <- run_cci_analysis(
  ...,
  nichenet_data_dir = "/path/to/nichenet/data"
)
```

## 다음 단계

1. 기본 테스트 성공 후 다양한 클러스터 조합 테스트
2. 실제 DEG 분석 결과와 통합
3. 결과 시각화 및 해석
4. 다른 CCI 도구와 비교 (CellChat 등)

