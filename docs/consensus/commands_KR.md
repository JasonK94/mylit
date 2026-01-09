
# DEG Consensus 방법론 사용 가이드

이 문서는 DEG Consensus 파이프라인의 개별 방법론 함수, 통합 Wrapper 함수, 그리고 CLI 스크립트 사용 방법을 설명합니다.

## 1. 주요 파라미터 설명

### `cluster_col` (클러스터 ID)
분석을 수행할 세포 집단(Subpopulation)을 정의하는 컬럼입니다. 분석은 **각 클러스터별로** 수행됩니다.
- **배경**: DEG 분석은 보통 특정 세포 유형(예: "B cells") 내에서 조건 간의 차이를 찾기 위해 수행됩니다.
- **Pseudobulk**: Pseudobulk 기반 방법론에서는 `(sample_id, cluster_id)` 단위로 세포 데이터를 합산(aggregation)합니다.
- **예시**: `anno3` (세포 유형), `seurat_clusters` (숫자 클러스터).

### `group_col` (그룹 ID)
발현 차이를 분석할 주 관심 변수입니다 (예: 환자군 vs 대조군).
- 최소 2개 이상의 레벨이 있어야 합니다 (대비 분석을 위해 범주형 변수 권장).
- **예시**: `g3` (뇌졸중 유형), `diagnosis` (질병 유무).

### `fixed_effects` vs `random_effects` (고정 효과 vs 변량 효과)
- **Fixed Effects (고정 효과)**: 모델에서 보정하고자 하는 공변량입니다 (예: `age`, `sex`, `batch`).
    - *범주형/연속형 모두 사용 가능한가요?*: 네. `model.matrix`가 자동으로 처리합니다.
    - **Batch**: `GEM`과 같은 배치 정보가 있다면 고정 효과 모델(edgeR, DESeq2, limma)에서는 여기에 포함시킵니다.
- **Random Effects (변량 효과)**: 혼합 모형(mixed models: dream, nebula, MAST-glmer)에서 피험자(Subject) 내의 상관성을 보정하기 위해 사용합니다.
    - **표준**: `(1|sample_id)` 또는 `(1|hos_no)`.
    - **주의**: 피험자에 대한 변량 효과를 사용하는 경우, 배치(Batch)가 피험자와 완전히 겹친다면(Collinear, 즉 한 피험자는 하나의 배치에만 속함) 배치를 별도로 넣지 않는 것이 일반적입니다. (Mixed Model에서는 피험자 효과가 배치를 포함할 수 있음).
    - **권장 사항**: `dream`/`nebula` 사용 시 `random_effect="(1|hos_no)"`를 사용하고, `GEM`은 고정 효과에서 제외하십시오.

## 2. R 함수 사용 예시

### 설정 (Setup)
```{r}
library(Seurat)
source("myR/R/deg_consensus/run_deg_consensus.R")

sobj <- qs::qread("/data/user3/sobj/is2_IS_3_1_plots.qs")

# 설정
cluster_col <- "anno3big"
group_col <- "g3"
contrast <- "2 - 1" # 그룹 2 vs 그룹 1 (g3 컬럼의 레벨)
fixed_covars <- c("age", "sex", "GEM") 
mixed_covars <- c("age", "sex")       
random_var <- "hos_no"
```

### A. Pseudobulk 방법론 (고정 효과)
**지원**: edgeR, DESeq2, limma-voom, limma-trend.

```{r}
# 1. edgeR (LRT)
res_edger <- runEDGER_LRT_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_covars, contrast=contrast
)

# 2. DESeq2 (LRT)
res_deseq2 <- runDESEQ2_LRT_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_covars, contrast=contrast
)

# 3. limma-voom
res_voom <- runLIMMA_voom_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_covars, contrast=contrast
)

# 4. limma-trend (muscat 스타일)
res_trend <- runLIMMA_trend_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_covars, contrast=contrast
)
```

### B. 혼합 모형 방법론 (변량 효과)
**지원**: dream, nebula, MAST.

```{r}
# 5. dream (Pseudobulk LMM)
# 참고: runDREAM_v1은 sample_id를 사용하여 pseudobulking을 수행하며, 
# (1|sample_id) 효과를 모델에 포함합니다.
res_dream <- runDREAM_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  sample_id=random_var,          # Aggregation 및 Random Effect ID
  covar_effects=mixed_covars,    # 고정 효과: age, sex (GEM 제외)
  batch_id=NULL,                 
  contrast=contrast
)

# 6. Nebula (Single-Cell NBMM)
res_nebula <- runNebula_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  sample_id=random_var,          # Random effect 그룹 변수
  covar_effects=mixed_covars,
  contrast=contrast
)

# 7. MAST (Single-Cell Hurdle)
res_mast <- runMAST_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  sample_id=random_var,
  covar_effects=mixed_covars,
  contrast=contrast,
  random_effect=random_var       # 이 값이 설정되면 glmer(혼합모형) 모드로 작동 (느림)
)
```

## 3. 통합 Wrapper 함수 (All-in-One)

`run_deg_consensus` 함수는 여러 방법론을 한 번에 실행하고 결과 형식을 표준화합니다.

```{r}
res_cons <- run_deg_consensus(
  sobj,
  methods = c("edgeR-LRT", "limma-voom", "nebula"), # 실행할 방법론 선택
  cluster_id = cluster_col,
  group_id = group_col,
  sample_id = random_var,
  covar_effects = fixed_covars, # 주의: 혼합 모형과 고정 효과 모형을 섞어 쓸 때 공변량 처리에 주의 필요.
                                # 현재 Wrapper는 입력된 covar_effects를 모든 방법론에 동일하게 적용함.
                                # 정밀한 제어가 필요하면 CLI를 사용하여 분리 실행 권장.
  contrast = contrast,
  n_cores = 4
)
```

## 4. CLI 사용법

터미널에서 `scripts/consensus/run_deg_consensus_cli.R` 스크립트를 통해 분석을 실행할 수 있습니다.

### 문법 (Syntax)
```bash
Rscript scripts/consensus/run_deg_consensus_cli.R \
  --input [PATH_TO_SEURAT_QS] \
  --output [OUTPUT_DIR] \
  --cluster [CLUSTER_COL] \
  --group [GROUP_COL] \
  --covariates [COVAR_STRING] \
  --methods [METHODS_STRING] \
  --contrast [CONTRAST] \
  --cores [N] \
  --sampleid [SAMPLE_ID]
```

### 예시 1: 고정 효과 모형 실행 (edgeR, DESeq2 etc.)
```bash
Rscript scripts/consensus/run_deg_consensus_cli.R \
  --input data/sobj.qs \
  --output results/fixed_run \
  --cluster anno3big \
  --group g3 \
  --covariates "sex,age,GEM" \
  --methods "edgeR-LRT,DESeq2-Wald,limma-voom" \
  --contrast "2 - 1" \
  --cores 8
```

### 예시 2: 혼합 모형 실행 (Nebula, Dream)
```bash
Rscript scripts/consensus/run_deg_consensus_cli.R \
  --input data/sobj.qs \
  --output results/mixed_run \
  --cluster anno3big \
  --group g3 \
  --covariates "sex,age" \
  --sampleid "hos_no" \
  --methods "nebula,dream" \
  --contrast "2 - 1" \
  --cores 8
```
