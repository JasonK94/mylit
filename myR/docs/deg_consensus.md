# Multi-Model DEG Consensus (deg-consensus)

## 1. 개요

- **목적**: 여러 DEG 방법론(limma, edgeR, DESeq2, muscat, nebula 등)을 **같은 데이터**에 적용해서  
  - 방법론별 DEG 결과를 수집하고  
  - 형식을 표준화한 뒤  
  - 방법론 간 일치도/클러스터링을 기반으로 **Consensus DEG signature**를 만드는 엔진.
- **입력**: Seurat object (예: `IS6_sex_added_251110.qs`, `IS6_sex_added_251110_ds2500.qs`)
- **핵심 단계**:
  1. 여러 DEG 방법 실행 (`run_deg_consensus`)
  2. 결과 표준화 (`standardize_deg_results`)
  3. gene × method 행렬 구성 (`build_deg_matrices`)
  4. 방법론 간 일치도/클러스터링 (`compute_agreement_scores`, `perform_deg_pca`, `cluster_deg_methods`)
  5. 최종 Consensus DEG 리스트 (`compute_consensus_scores`, `generate_consensus_deg_list`)

---

## 2. 주요 코드 파일 (어디에 무엇이 있는지 한눈에)

> 경로는 repo root 기준: `myR/R/deg_consensus/*.R`, 스크립트는 `_wt/deg-consensus/scripts/*.R`

### 2.1 DEG 메인/핵심 함수들

- `myR/R/deg_consensus/run_deg_consensus.R`
  - **`run_deg_consensus()`**
  - 여러 DEG 방법을 한 번에 실행하는 **통합 엔진**.
  - 지원 방법:
    - `muscat-edgeR`, `muscat-DESeq2`, `muscat-limma-voom`, `muscat-limma-trend`
    - `nebula`
    - `limma-voom`, `limma-trend`, (`limma-wt`, `dream`는 향후 확장)
    - `edgeR-LRT`, `edgeR-QLF`, (`edgeR-robust` 향후)
    - `DESeq2-Wald`, `DESeq2-LRT`
  - 공통 인자:  
    - `sobj`: Seurat object  
    - `contrast`: 예 `"2 - 1"`  
    - `cluster_id`, `sample_id`, `group_id`, `batch_id`  
    - `pb_min_cells`: pseudobulk 최소 셀 수  
    - `remove_na_groups`: NA 그룹 제거 여부
  - 반환:
    - `results`: 방법론별 raw DEG 결과 (리스트)
    - `methods_run`, `methods_failed`, `errors`

- `myR/R/test_analysis.R`
  - **`runMUSCAT2_v1()`**, **`runNEBULA2_v1()`**
  - 원래 있던 muscat/nebula 기반 DEG pipeline.
  - `run_deg_consensus()`의 muscat/nebula 래퍼에서 호출.

---

### 2.2 방법론별 구현 파일

- `myR/R/deg_consensus/deg_methods_limma.R`
  - **`runLIMMA_voom_v1()`**, **`runLIMMA_trend_v1()`**
  - `runMUSCAT2_v1`와 동일한 인터페이스로 만든 limma 계열.
  - **전처리/로직**:
    - `Seurat -> SCE -> muscat::prepSCE`
    - `muscat::aggregateData(sce, by = c("cluster_id","sample_id"))`로 **클러스터×샘플 pseudobulk** 생성
    - batch confounding (group와 batch가 완전히 겹치는 경우) 자동 체크 후 batch 제외
    - library size 0인 샘플 제거 후 **다시 design/contrast matrix 생성**
  - 클러스터별로 독립적인 limma-voom / limma-trend 분석 수행, cluster별 DEG 결과를 하나의 data.frame으로 반환.

- `myR/R/deg_consensus/deg_methods_edger.R`
  - **`runEDGER_LRT_v1()`**, **`runEDGER_QLF_v1()`** (향후 `runEDGER_robust_v1`)
  - limma 버전과 동일한 인터페이스/전처리, 대신 edgeR GLM(LRT/QLF) 사용.
  - **전처리**는 limma 버전과 거의 동일:
    - NA 제거, `prepSCE`, `aggregateData`로 pseudobulk
    - library size 0인 샘플 제거 + batch confounding 재평가

- `myR/R/deg_consensus/deg_methods_deseq2.R`
  - **`runDESEQ2_Wald_v1()`**, **`runDESEQ2_LRT_v1()`**
  - 동일한 pseudobulk 전처리를 거친 뒤 클러스터별로 DESeq2를 적용.
  - 주의: DESeq2는 **정수 카운트**가 필요하므로, pseudobulk가 비정수인 경우 rounding 등 전처리가 필요할 수 있음(현재 데이터 특성 상 제약이 있을 수 있음).

- `myR/R/deg_consensus/deg_methods_base.R`
  - 기본/공통 유틸 및 NEBULA pseudobulk 변형 등 일부 실험적 함수가 포함.
  - 주로 내부 개발용/참고용으로 사용 (현재 메인 파이프라인은 위 3개 파일 중심).

---

### 2.3 결과 표준화 및 Consensus 분석

- `myR/R/deg_consensus/deg_standardize.R`
  - **`standardize_deg_results()`**
    - 각 방법론의 결과 형식이 제각각인 것을 **공통 포맷**으로 변환.
    - 공통 컬럼:
      - `gene`, `logFC`, `pvalue`, `pvalue_adj`, `statistic`, (선택) `se`, `cluster_id`, `method`.
    - p-value, adj p-value 컬럼명 자동 인식:
      - p: `pvalue`, `p_value`, `PValue`, `p_val` 등
      - padj: `pvalue_adj`, `p_val_adj`, `padj`, `FDR`,  
        특히 muscat: **`p_adj.loc`(local FDR) > `p_adj.glb`(global FDR)** 우선 사용.
  - **`build_deg_matrices()`**
    - 표준화된 결과들을 합쳐서 **gene × method 행렬** 구성:
      - `beta`, `pvalue`, `logp`, `significance`, `statistic`, `se`.
    - 같은 유전자가 여러 클러스터에 나타나는 경우, **가장 유의한 결과(최소 FDR 또는 p)**를 대표로 사용.

- `myR/R/deg_consensus/deg_consensus_analysis.R`
  - **`compute_agreement_scores(significance_matrix)`**
    - 각 유전자에 대해 **방법론 간 유의 여부 일치도** (0~1) 계산.
  - **`perform_deg_pca(deg_matrices)`**
    - 방법론별 특성(평균 |beta|, 평균 logp, 유의 비율, 유의한 유전자 수)을 사용해 **PCA로 방법론 공간 분석**.
  - **`cluster_deg_methods(deg_matrices, method = "hierarchical")`**
    - 같은 특성들로 방법론 간 **계층적 클러스터링 또는 k-means** 수행.
  - **`compute_consensus_scores(deg_matrices, agreement_scores, ...)`**
    - 각 유전자에 대해:
      - `mean_beta`, `mean_logp`, `agreement`, `n_significant`, `n_methods`
      - **weighted beta** (방법론별 가중치 반영 가능)
      - `consensus_score = agreement * |weighted_beta|`
  - **`generate_consensus_deg_list(consensus_scores, fdr_threshold, agreement_threshold, min_methods)`**
    - 최종 Consensus DEG 필터링:
      - `agreement >= agreement_threshold`
      - `n_significant >= min_methods` (최소 몇 개 방법에서 유의해야 하는지)
      - (fdr_threshold는 주로 upstream FDR 기준에 맞춰 해석용)

---

## 3. 실행/테스트용 스크립트

> 위치: `/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts`

- `run_consensus_simple.R`
  - **R 콘솔에서 바로 실행 가능한 최소 버전**.
  - 전제: Seurat 객체(`is5` 등)가 이미 메모리에 있음.
  - 하는 일:
    1. `devtools::load_all(".../mylit/myR")` + deg_consensus 관련 R 파일들 `source()`
    2. 적당한 방법론 셋 (muscat + limma + edgeR + DESeq2) 실행 (`run_deg_consensus`)
    3. `standardize_deg_results` → `build_deg_matrices`
    4. `compute_agreement_scores` → `compute_consensus_scores` → `generate_consensus_deg_list`
    5. `fdr_threshold = 0.1`, `agreement_threshold = 0.3`, `min_methods = 2` 기본
    6. 최종 결과를 `"/data/user3/sobj/deg_consensus_final_result.qs"`로 저장.

- `run_consensus_analysis.R`
  - **조금 더 자세한 전체 파이프라인 스크립트**.
  - PCA/클러스터링 결과까지 포함해서 요약 출력.
  - FDR/threshold는 기본적으로 조금 더 엄격 (`fdr_threshold = 0.05`, `agreement_threshold = 0.5`, `min_methods ≈ 절반 방법 수`).

- `test_step_by_step.R`
  - **디버깅/검증용 단계별 스크립트**.
  - 순서:
    1. 함수 로드 (각 파일 `source`)  
    2. 데이터 로드 (downsampled `IS6_sex_added_251110_ds2500.qs`)  
    3. 개별 방법론 테스트 (예: `muscat-edgeR`, `limma-voom`, `edgeR-LRT`)  
    4. 여러 방법론 동시에 실행 (`run_deg_consensus`)  
    5. 표준화/행렬/consensus 단계별 확인  
    6. 필요 시 `traceback()` 호출 (interactive 환경에서)

- 그 외:
  - `test_phase2_limma.R`, `test_phase2_limma_interactive.R`: limma 계열 단독 테스트용
  - `test_full_pipeline.R`: 선택된 일부 방법론으로 full pipeline 테스트용

---

## 4. R 콘솔에서 직접 테스트하는 최소 예제

### 4.1 다운샘플 데이터로 빠른 테스트

# 0. 패키지
library(qs)
library(devtools)

# 1. Seurat 객체 로드 (downsampled)
is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 2. 간단 스크립트 실행 (모든 함수와 파이프라인 한 번에)
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/run_consensus_simple.R")
# -> 콘솔에 방법론별 요약 + Consensus DEG 개수 출력
# -> 결과 객체는 /data/user3/sobj/deg_consensus_final_result.qs 로 저장됨### 4.2 본 데이터(전체)로 실행

library(qs)
library(devtools)

# 1. 전체 데이터 로드
is5 <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# 2. 동일하게 실행
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/run_consensus_simple.R")### 4.3 파라미터(FDR, min_methods 등) 조금 만져보고 싶을 때

- `run_consensus_simple.R` 안에서:
  - `fdr_threshold` (기본 0.1): 낮추면 더 엄격, 올리면 더 많은 DEG.
  - `agreement_threshold` (기본 0.3): 방법론 간 일치도 기준.
  - `min_methods` (기본 2): 최소 몇 개 방법에서 유의해야 하는지.
- 또는, R 콘솔에서 직접 실행하고 싶으면:

# 1) 함수 로드
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 2) 메타데이터 키 설정 (is5가 이미 로드되어 있다고 가정)
cluster_key <- if ("anno3.scvi" %in% colnames(is5@meta.data)) "anno3.scvi" else "seurat_clusters"
sample_key  <- "hos_no"
group_key   <- "g3"
batch_key   <- if ("GEM" %in% colnames(is5@meta.data)) "GEM" else NULL
contrast_str <- "2 - 1"

# 3) 여러 방법론 실행
methods_to_run <- c(
  "muscat-edgeR",
  "muscat-DESeq2",
  "muscat-limma-voom",
  "muscat-limma-trend",
  "limma-voom",
  "limma-trend",
  "edgeR-LRT",
  "edgeR-QLF",
  "DESeq2-Wald",
  "DESeq2-LRT"
)

result_consensus <- run_deg_consensus(
  sobj = is5,
  contrast = contrast_str,
  methods = methods_to_run,
  cluster_id = cluster_key,
  sample_id = sample_key,
  group_id  = group_key,
  batch_id  = batch_key,
  remove_na_groups = TRUE,
  verbose = TRUE
)

# 4) 표준화 + 행렬 + consensus
standardized_results <- lapply(result_consensus$methods_run, function(m) {
  standardize_deg_results(result_consensus$results[[m]], m)
})
names(standardized_results) <- result_consensus$methods_run

fdr_threshold <- 0.1
deg_matrices <- build_deg_matrices(standardized_results, fdr_threshold = fdr_threshold)

agreement_scores <- compute_agreement_scores(deg_matrices$significance)
consensus_scores <- compute_consensus_scores(deg_matrices, agreement_scores)

consensus_deg_list <- generate_consensus_deg_list(
  consensus_scores,
  fdr_threshold      = fdr_threshold,
  agreement_threshold = 0.3,
  min_methods         = 2
)

# 5) 결과 간단 확인
nrow(consensus_deg_list)
head(consensus_deg_list)---

## 5. 상태 메모

- Phase 1~5 (통합 함수, 방법론별 구현, 결과 표준화, consensus 분석)는 구현 및 기본 테스트 완료.
- **FDR/유의성 문제**:
  - muscat의 `p_adj.loc`(local FDR)을 우선 사용하도록 수정하여 유의 유전자 수를 늘림.
  - 클러스터별 결과가 겹치는 유전자에 대해서는 가장 유의한 결과만 사용하도록 `build_deg_matrices()` 수정.
  - `run_consensus_simple.R`에서 기본 `fdr_threshold = 0.1`, `agreement_threshold = 0.3`, `min_methods = 2`로 설정해 조금 더 완화된 기준 사용.
- Phase 6 (시각화/결과 출력 전용 함수)는 추후 별도 구현 예정.