# CCI 분석 도구 빠른 시작 가이드

## R 세션에서 실행

### 1. R 세션 시작
```bash
cd /home/user3/GJC_KDW_250721
R
```

### 2. 함수 로드
```r
# CCI 모듈 함수들
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")

# 기존 NicheNet 함수 (워크트리 우선)
cci_core_worktree <- "/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R"
cci_core_mainrepo <- "/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R"
if (file.exists(cci_core_worktree)) {
  source(cci_core_worktree)
} else if (file.exists(cci_core_mainrepo)) {
  source(cci_core_mainrepo)
} else {
  stop("CCI.R not found in worktree or main repository.")
}

# 필요한 패키지
library(Seurat)
library(dplyr)
library(qs)
```

### 3. 데이터 로드
```r
# 테스트 데이터
sobj_test <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 메타데이터 확인
colnames(sobj_test@meta.data)
table(sobj_test@meta.data$anno3.scvi, useNA = "ifany")
table(sobj_test@meta.data$g3, useNA = "ifany")
```

### 4. DEG 리스트 준비
```r
# 클러스터 확인
clusters <- unique(sobj_test@meta.data$anno3.scvi)
clusters <- clusters[!is.na(clusters)]
receiver_cluster_test <- clusters[1]  # 첫 번째 클러스터 사용

# 예시 DEG 데이터프레임 (실제로는 분석 결과 사용)
deg_df_example <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
  cluster = rep(receiver_cluster_test, 5),
  avg_log2FC = c(1.5, 2.0, -1.2, 0.8, -0.5),
  p_val_adj = c(0.001, 0.0001, 0.01, 0.05, 0.1),
  stringsAsFactors = FALSE
)
```

### 5. CCI 분석 실행
```r
results <- run_cci_analysis(
  sobj = sobj_test,
  cluster_col = "anno3.scvi",
  deg_df = deg_df_example,
  receiver_cluster = receiver_cluster_test,
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  species = "human",
  verbose = TRUE,
  auto_save = TRUE
)
```

### 5-b. `run_nichenet_analysis()`를 직접 호출하며 `receiver_de_table` 재사용
다운샘플 테스트에서 추출한 receiver DEG 테이블을 `.qs`로 저장했다면 동일한 데이터를 다시 계산하지 않고 그대로 사용할 수 있습니다.

```r
receiver_deg_qs <- "/data/user3/sobj/receiver_CD4_deg_table.qs"
receiver_deg_table <- qs::qread(receiver_deg_qs)

nichenet_results <- run_nichenet_analysis(
  seurat_obj = sobj_test,
  species = "human",
  sender_celltypes = c("Cluster2", "Cluster3"),
  receiver_celltype = receiver_cluster_test,
  assay_name = "SCT",
  cluster_col = "anno3.scvi",
  receiver_DE_ident1 = "2",
  receiver_DE_ident2 = "1",
  receiver_DE_group_by = "g3",
  receiver_de_table = receiver_deg_table,
  receiver_gene_col = "gene",
  receiver_logfc_col = "avg_log2FC",
  receiver_pval_col = "p_val_adj",
  verbose = TRUE
)
```
`receiver_logfc_col`과 `receiver_pval_col`을 통해 외부 DEG 분석 엔진에서 생성한 컬럼명을 그대로 지정할 수 있습니다.

### 6. 결과 확인
```r
# 결과 구조 확인
str(results, max.level = 2)

# Top ligands 확인
results$nichenet_results$best_upstream_ligands

# Ligand activities 확인
head(results$nichenet_results$ligand_activities)

# 저장된 경로 확인
results$saved_path
```

## 자동 테스트 실행

대화형 테스트 스크립트 실행:
```r
source("/home/user3/data_user3/git_repo/_wt/cci/scripts/cci/test_cci_interactive.R")
```

## 문제 해결

### 오류: "Cannot find run_nichenet_analysis function"
```r
# CCI.R 파일을 직접 소스
source("/home/user3/data_user3/git_repo/mylit/myR/R/CCI.R")
```

### 오류: "No DEGs found"
- `deg_df`의 `cluster` 컬럼 값이 `receiver_cluster`와 정확히 일치하는지 확인
- `deg_df`에 실제 유전자명이 있는지 확인

### 오류: "No potential ligands found"
- Sender 클러스터에 ligand가 발현되는지 확인
- `min_pct_expressed` 값을 낮춰보기 (예: 0.05)

## 다음 단계

1. 실제 DEG 분석 결과 사용
2. 다양한 클러스터 조합 테스트
3. 결과 시각화 및 해석

