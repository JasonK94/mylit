# myR 패키지

## 개요

`myR`는 단일세포 RNA 시퀀싱(scRNA-seq) 분석을 위한 종합 R 패키지로, GeoMx 공간 전사체 데이터도 전문적으로 지원합니다. 이 패키지는 차등 발현 분석, 의사시간 추론, 경로 풍부도 분석, 세포 간 상호작용 분석, 시각화 도구 등을 제공합니다.

## 설치

```r
# devtools가 설치되어 있지 않다면 먼저 설치합니다
install.packages("devtools")

# 로컬 디렉터리에서 패키지 설치
devtools::install("path/to/myR", dependencies = TRUE)
```

## 주요 기능

- **차등 발현 분석**: 의사벌크 및 LMM 기반 DEG 분석
- **경로 풍부도 분석**: GO, KEGG, GSEA 분석
- **의사시간 분석**: Slingshot과 Monocle3를 이용한 궤적 추론
- **세포 간 상호작용**: NicheNet 기반 상호작용 분석
- **모듈 스코어링**: PROGENy 및 기타 시그니처 스코어링
- **GeoMx 분석**: NanoString GeoMx 데이터를 위한 특화 도구

## 주요 컴포넌트

### 핵심 분석 모듈

- **pseudobulk_deg.R**: 의사벌크 차등 발현 분석
- **pathway_analysis.R**: 경로 풍부도 분석(GO, KEGG, GSEA)
- **pseudotime.R**: 의사시간 및 궤적 분석
- **CCI.R**: 세포 간 상호작용 분석
- **deg.R**: 모듈 스코어링 및 DEG 분석

### 유틸리티 모듈

- **plots.R**: 포괄적인 시각화 함수
- **utils.R**: 데이터 처리를 위한 유틸리티 함수
- **gene_list.R**: 유전자 리스트 관리 및 조작
- **cluster_frequency.R**: 클러스터 비율 분석

## 함수 비교 표

### 핵심 분석 함수

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| `prepare_geomx_data` | count 파일, 메타데이터 파일 | 리스트(발현 행렬, 메타데이터, gene_info) | GeoMx 데이터를 분석용 형식으로 전처리 |
| `q3_normalize` | expr_matrix, scaling_factor | 정규화 행렬(log2) | 3사분위(Q3) 기반 발현 정규화 |
| `find_deg_geomx` | norm_expr, metadata, group_var | DEG 결과 데이터프레임 | limma/edgeR를 이용한 DEG 탐색 |
| `run_lmm_multiple_genes` | seurat_obj, genes, config | LMM 결과 | 여러 유전자에 대한 선형혼합모형 분석 |
| `find_response_differential_genes` | lmm_summary, pval_threshold | 유의 유전자 데이터프레임 | 처치 반응성이 유의한 유전자 탐색 |

### 의사벌크 및 DEG 함수

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| `run_pseudobulk_deg` | seurat_obj, analysis_level, cluster_group, condition_col | DEG 결과 | 의사벌크 차등 발현 분석 수행 |
| `prepare_pseudobulk_edgeR` | seurat_obj, cluster_group, sample_col, counts_assay | 의사벌크 행렬 | edgeR 분석을 위한 의사벌크 데이터 준비 |
| `cluster_pseudobulk_deg` | sobj, cluster_group, condition_col, genes, ... | 클러스터별 DEG 결과 | 클러스터 수준의 의사벌크 DEG 분석 |
| `pseudobulk_linear_fit` | sobj, genes, sample_col, numeric_predictor, ... | 선형 피팅 결과 | 의사벌크 데이터의 선형 모델 피팅 |

### 경로 분석 함수

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| `myGO` | DEG 데이터프레임, pathway_set, analysis_type | 경로 결과 리스트 | GO/GSEA 기반 경로 풍부도 분석 |
| `run_go_analysis` | genes, gene_type, ont, ... | GO 결과 | Gene Ontology 풍부도 분석 |
| `run_kegg_analysis` | genes, background, pval_cutoff | KEGG 결과 | KEGG 경로 풍부도 분석 |
| `run_gsea_analysis` | ranked_genes, pathway_set, ... | GSEA 결과 | 유전자 집합 풍부도 분석 |

### 모듈 스코어링 및 시그니처

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| `AddMultipleModuleScores` | seurat_object, gene_modules | 모듈 스코어가 추가된 Seurat 객체 | 여러 모듈 스코어를 한 번에 계산 |
| `add_progeny_scores` | seurat_obj, organism, topn | PROGENy 스코어가 추가된 Seurat 객체 | PROGENy 경로 활성도 스코어 추가 |
| `add_signature_enrichit` | seurat_obj, gene_sets, ... | 시그니처 스코어가 추가된 Seurat 객체 | EnrichIt 기반 시그니처 스코어 추가 |
| `find_gene_signature` | data, signature_type, ... | 유전자 시그니처 리스트 | 데이터에서 시그니처를 추출 |

### 시각화 함수

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| `plot_volcano` | lmm_summary, ... | 화산도 | LMM 결과를 시각화 |
| `PlotModuleScoreHeatmap` | seurat_object, features, assays | 히트맵 | 모듈 스코어 히트맵 생성 |
| `plot_cluster_fractions` | sobj_metadata, cluster_col, ... | 클러스터 비율 그래프 | 클러스터 비율 시각화 |
| `plot_interaction_for_gene` | sobj, gene, patient, treatment, timepoint | 상호작용 플롯 | 유전자 발현 상호작용 플롯 |
| `myhm_genes4` | sobj, features, group.by, ... | 히트맵 | 유전자 발현 히트맵 생성 |
| `cdf` | data, probability_col, ratio_col | 누적 분포 그래프 | 누적 분포 시각화 |

### 의사시간 분석 함수

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| `run_slingshot_from_seurat` | seurat_obj, cluster_col, reduced_dim_name | SingleCellExperiment 객체 | Slingshot 의사시간 추론 실행 |
| `analyze_gene_dynamics` | gene_id, cds_obj, condition_col_name, ... | 분석 결과 | 의사시간 축에서 유전자 동역학 분석 |
| `process_gene_list_dynamics` | gene_list, cds_obj, condition_col_name, ... | 결과 리스트 | 여러 유전자의 동역학 분석 |
| `analyze_gene_dynamics_tradeSeq` | gene_id, cds_obj, condition_col_name, ... | tradeSeq 결과 | tradeSeq 기반 의사시간 분석 |

### 세포 간 상호작용 함수

| 함수명 | 입력 | 출력 | 설명 |
|---|---|---|---|
| CCI.R 내 함수들 | 다중 | NicheNet 분석 결과 | NicheNet 기반 상호작용 분석 |
| `run_nichenet_analysis` | sobj, sender, receiver, condition_col | 상호작용 결과 | 리간드-수용체 분석 실행 |

## 사용 예시

### 의사벌크 DEG 분석

```r
# 패키지 로드
library(myR)

# 의사벌크 DEG 분석 실행
deg_results <- run_pseudobulk_deg(
  seurat_obj = seurat_object,
  analysis_level = "per_cluster",
  cluster_group = "cell_type",
  condition_col = "treatment"
)
```

### 경로 풍부도 분석

```r
# 경로 풍부도 분석 수행
pathway_results <- myGO(
  DEG = deg_results,
  analysis_type = "ALL",
  fc_threshold = 0.25,
  pval_threshold = 0.05
)
```

### 모듈 스코어링

```r
# PROGENy 경로 스코어 추가
seurat_obj <- add_progeny_scores(
  seurat_obj = seurat_object,
  organism = "Human",
  topn = 100
)
```

## 의존성

- Seurat
- edgeR
- limma
- slingshot
- monocle3
- NicheNet
- PROGENy
- clusterProfiler
- msigdbr
- dplyr
- ggplot2

## 문서

각 함수의 상세 문서는 R 도움말에서 확인할 수 있습니다:

```r
?run_pseudobulk_deg
?myGO
?add_progeny_scores
```

## 개발자

- Development team

## 라이선스

자세한 내용은 LICENSE 파일을 참고하세요.


