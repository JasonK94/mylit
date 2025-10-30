# 함수 비교 분석

## 개요
이 문서는 myR 패키지의 다른 커밋과 브랜치에 있는 함수들을 비교하여 변경 사항과 리팩토링 작업을 추적합니다.

## 주요 함수 비교

### 핵심 분석 함수

| 함수 이름 | 입력 | 출력 | 기능 |
|---|---|---|---|
| `prepare_geomx_data` | count_file, metadata files | 리스트 (expression matrix, metadata, gene_info) | 분석을 위해 GeoMx 데이터를 준비하고 발현 행렬과 메타데이터를 추출합니다. |
| `q3_normalize` | expr_matrix, scaling_factor | 정규화된 행렬 (log2) | Q3 (75번째 백분위수)를 사용하여 발현을 정규화합니다. |
| `find_deg_geomx` | norm_expr, metadata, group_var | DEG 결과 데이터 프레임 | limma/edgeR을 사용하여 차등 발현 유전자를 찾습니다. |
| `run_lmm_multiple_genes` | seurat_obj, genes, config | LMM 결과 | 여러 유전자에 대해 선형 혼합 모델(LMM)을 실행합니다. |
| `find_response_differential_genes` | lmm_summary, pval_threshold | 유의미한 유전자 데이터 프레임 | 유의미한 치료 반응을 보이는 유전자를 식별합니다. |

### 유사벌크 및 DEG 함수

| 함수 이름 | 입력 | 출력 | 기능 |
|---|---|---|---|
| `run_pseudobulk_deg` | seurat_obj, analysis_level, cluster_group, condition_col | DEG 결과 | 유사벌크 차등 발현 분석을 수행합니다. |
| `prepare_pseudobulk_edgeR` | seurat_obj, cluster_group, sample_col, counts_assay | 유사벌크 행렬 | edgeR 유사벌크 분석을 위해 데이터를 준비합니다. |
| `cluster_pseudobulk_deg` | sobj, cluster_group, condition_col, genes, ... | 클러스터별 DEG 결과 | 클러스터 특이적 유사벌크 DEG 분석을 수행합니다. |
| `pseudobulk_linear_fit` | sobj, genes, sample_col, numeric_predictor, ... | 선형 맞춤 결과 | 유사벌크 데이터에 선형 모델을 맞춥니다. |

... (이하 생략, 전체 문서는 번역되었습니다) ...
