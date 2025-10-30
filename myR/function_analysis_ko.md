# myR 함수 분석

## 디렉토리: `core`

### 파일: `myR/R/core/data_preparation.R`

#### 함수: `.get_feature_vector`
- **설명:** Seurat 객체에서 특정 피처(유전자 또는 메타데이터 컬럼)를 추출하는 내부 헬퍼 함수입니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `feature`: 피처 이름 (유전자 또는 메타데이터 컬럼)
  - `assay`: 사용할 Assay (기본값: DefaultAssay)
  - `slot`: 사용할 Slot (기본값: "data")
  - `cells`: 추출할 세포의 일부 (기본값: 모든 세포)
- **반환값:** 피처 값으로 이루어진 숫자 벡터

#### 함수: `get_feature_vec`
- **설명:** Seurat 객체에서 특정 피처(유전자 또는 메타데이터 컬럼)를 추출합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `feature`: 피처 이름 (유전자 또는 메타데이터 컬럼)
  - `assay`: 사용할 Assay (기본값: DefaultAssay)
  - `slot`: 사용할 Slot (기본값: "data")
  - `cells`: 추출할 세포의 일부 (기본값: 모든 세포)
- **반환값:** 피처 값으로 이루어진 숫자 벡터

#### 함수: `get_feature_matrix`
- **설명:** Seurat 객체에서 여러 피처를 추출하여 행렬로 반환합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `features`: 피처 이름으로 이루어진 문자 벡터
  - `assay`: 사용할 Assay (기본값: DefaultAssay)
  - `slot`: 사용할 Slot (기본값: "data")
  - `cells`: 추출할 세포의 일부 (기본값: 모든 세포)
- **반환값:** 피처가 행, 세포가 열인 행렬

#### 함수: `prepare_metadata_table`
- **설명:** Seurat 객체에서 메타데이터를 추출하고 선택적으로 필터링합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `columns`: 추출할 특정 컬럼 (NULL = 모든 컬럼)
  - `cells`: 일부 세포 (NULL = 모든 세포)
  - `drop_na`: NA 값이 있는 행을 제거할지 여부 (기본값: FALSE)
- **반환값:** 메타데이터 데이터 프레임

#### 함수: `aggregate_expression_by_group`
- **설명:** 그룹(예: 클러스터, 샘플) 내의 세포들에서 유전자 발현을 집계합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `features`: 집계할 피처(유전자)
  - `group_by`: 그룹화할 기준이 되는 메타데이터 컬럼
  - `method`: 집계 방법: "mean", "median", "sum" (기본값: "mean")
  - `assay`: 사용할 Assay (기본값: DefaultAssay)
  - `slot`: 사용할 Slot (기본값: "data")
- **반환값:** 피처가 행, 그룹이 열인 행렬

#### 함수: `prepare_count_matrix`
- **설명:** Seurat 객체에서 카운트 행렬을 추출하며, 선택적으로 필터링을 적용합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `assay`: 사용할 Assay (기본값: "RNA")
  - `slot`: 사용할 Slot (기본값: "counts")
  - `features`: 포함할 피처 (NULL = 모두)
  - `cells`: 포함할 세포 (NULL = 모두)
  - `min_cells`: 피처가 발현되어야 하는 최소 세포 수 (기본값: 0)
  - `min_features`: 세포에서 발현되어야 하는 최소 피처 수 (기본값: 0)
- **반환값:** 희소(sparse) 또는 밀집(dense) 카운트 행렬

#### 함수: `convert_to_long_format`
- **설명:** 발현 데이터를 ggplot2에 적합한 긴 형식(tidy format)으로 변환합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `features`: 포함할 피처
  - `metadata_cols`: 포함할 메타데이터 컬럼
  - `assay`: 사용할 Assay (기본값: DefaultAssay)
  - `slot`: 사용할 Slot (기본값: "data")
  - `cells`: 포함할 세포 (NULL = 모두)
- **반환값:** cell, feature, expression, 그리고 메타데이터 컬럼을 포함하는 긴 형식의 데이터 프레임

#### 함수: `check_data_quality`
- **설명:** 발현 데이터에 대한 기본적인 품질 검사를 수행합니다.
- **파라미터:**
  - `data`: 발현 데이터의 행렬 또는 데이터 프레임
  - `check_finite`: 무한대 값 확인 (기본값: TRUE)
  - `check_na`: NA 값 확인 (기본값: TRUE)
  - `check_negative`: 음수 값 확인 (기본값: TRUE)
- **반환값:** 발견된 문제의 논리 플래그와 개수를 포함하는 리스트

### 파일: `myR/R/core/validation.R`

#### 함수: `validate_seurat`
- **설명:** 입력값이 예상되는 구성 요소를 가진 유효한 Seurat 객체인지 확인합니다.
- **파라미터:**
  - `obj`: 검증할 객체
  - `assay`: 확인할 선택적 Assay 이름
  - `reduction`: 확인할 선택적 차원 축소 이름
  - `min_cells`: 요구되는 최소 세포 수
  - `min_features`: 요구되는 최소 피처 수
- **반환값:** 유효하면 TRUE, 그렇지 않으면 오류 메시지와 함께 중단

#### 함수: `validate_metadata_column`
- **설명:** 메타데이터 컬럼이 존재하는지 확인하고 선택적으로 타입을 검증합니다.
- **파라미터:**
  - `obj`: Seurat 객체 또는 데이터 프레임
  - `column_name`: 검증할 컬럼 이름
  - `required_type`: 요구되는 선택적 타입 ("numeric", "factor", "character")
  - `allow_na`: NA 값을 허용할지 여부
- **반환값:** 유효하면 TRUE, 그렇지 않으면 오류 메시지와 함께 중단

#### 함수: `validate_genes`
- **설명:** 객체에 유전자가 존재하는지 확인하고 선택적으로 유효한 유전자로 필터링합니다.
- **파라미터:**
  - `obj`: Seurat 객체 또는 사용 가능한 유전자의 문자 벡터
  - `genes`: 검증할 유전자의 문자 벡터
  - `min_present`: 반드시 존재해야 하는 최소 유전자 수 (기본값: 모두)
  - `assay`: 유전자를 확인할 Assay (Seurat 객체용)
  - `warn_missing`: 누락된 유전자에 대해 경고할지 여부
- **반환값:** 객체에 존재하는 유효한 유전자의 문자 벡터

#### 함수: `validate_numeric_range`
- **설명:** 숫자 파라미터가 허용 가능한 범위 내에 있는지 확인합니다.
- **파라미터:**
  - `value`: 검증할 숫자 값
  - `param_name`: 파라미터 이름 (오류 메시지용)
  - `min`: 허용되는 최소값 (포함)
  - `max`: 허용되는 최대값 (포함)
  - `allow_na`: NA를 허용할지 여부
- **반환값:** 유효하면 TRUE, 그렇지 않으면 오류 메시지와 함께 중단

#### 함수: `validate_choice`
- **설명:** 값이 허용된 선택지 중 하나인지 검증합니다 (match.arg보다 더 많은 정보 제공).
- **파라미터:**
  - `value`: 검증할 값
  - `param_name`: 파라미터 이름 (오류 메시지용)
  - `choices`: 허용된 값의 벡터
  - `multiple`: 여러 선택지를 허용할지 여부
- **반환값:** 검증된 값 (multiple=TRUE인 경우 값들)

#### 함수: `validate_path`
- **설명:** 파일 경로가 존재하는지 확인하고 선택적으로 확장자를 검증합니다.
- **파라미터:**
  - `path`: 검증할 파일 경로
  - `must_exist`: 파일이 반드시 존재해야 하는지 여부
  - `extensions`: 허용된 확장자의 선택적 벡터 (예: c("csv", "txt"))
  - `type`: 경로 타입 ("file" 또는 "directory")
- **반환값:** 유효하면 정규화된 경로, 그렇지 않으면 오류 메시지와 함께 중단

#### 함수: `create_error_message`
- **설명:** 일관성 있고 유용한 오류 메시지를 생성하는 헬퍼 함수입니다.
- **파라미터:**
  - `context`: 오류가 발생한 컨텍스트 (예: 함수 이름)
  - `message`: 주 오류 메시지
  - `suggestion`: 오류 해결을 위한 선택적 제안
- **반환값:** 형식화된 오류 메시지

#### 함수: `check_packages`
- **설명:** 필요한 패키지가 설치되어 있는지 확인하고 선택적으로 로드합니다.
- **파라미터:**
  - `packages`: 패키지 이름의 문자 벡터
  - `load`: 패키지를 로드할지 여부 (기본값: FALSE)
- **반환값:** 모든 패키지를 사용할 수 있으면 TRUE, 그렇지 않으면 오류 메시지와 함께 중단

## 디렉토리: `analysis`

### 파일: `myR/R/analysis/cell_communication/nichenet_analysis.R`

#### 함수: `ligand_to_target`
- **설명:** NicheNet 데이터베이스에서 리간드-타겟 유전자 조절 잠재력을 검색합니다.
- **파라미터:**
  - `ligand`: 리간드 유전자의 문자 벡터
  - `target`: 타겟 유전자의 문자 벡터
  - `lr_network`: 리간드-수용체 네트워크 (NicheNet 제공)
  - `sig_network`: 신호 전달 네트워크 (NicheNet 제공)
  - `gr_network`: 유전자 조절 네트워크 (NicheNet 제공)
  - `top_n_targets`: 리간드당 반환할 상위 타겟 수 (기본값: 250)
- **반환값:** 리간드-타겟 조절 잠재력 점수를 포함하는 데이터 프레임

#### 함수: `run_nichenet_analysis`
- **설명:** 리간드 활성 예측, 수용체 추론 및 시각화를 포함한 포괄적인 NicheNet 세포-세포 상호작용 분석 워크플로우를 수행합니다.
- **파라미터:**
  - `seurat_obj`: Seurat 객체
  - `sender_cells`: 발신 세포 정체성의 벡터 또는 논리 벡터
  - `receiver_cells`: 수신 세포 정체성의 벡터 또는 논리 벡터
  - `condition_oi`: 수신 세포 DE 분석을 위한 관심 조건
  - `condition_ref`: 수신 세포 DE 분석을 위한 참조 조건
  - `ident_col`: 메타데이터의 정체성 컬럼 (기본값: "seurat_clusters")
  - `condition_col`: 메타데이터의 조건 컬럼 (DE에 필요)
  - `lr_network`: 리간드-수용체 네트워크 (NicheNet 제공)
  - `sig_network`: 신호 전달 네트워크 (NicheNet 제공)
  - `gr_network`: 유전자 조절 네트워크 (NicheNet 제공)
  - `ligand_target_matrix`: 사전 계산된 리간드-타겟 행렬 (선택 사항)
  - `expressed_pct`: 발현 비율 임계값 (기본값: 0.10)
  - `top_n_ligands`: 분석할 상위 리간드 수 (기본값: 20)
  - `top_n_targets`: 리간드당 상위 타겟 수 (기본값: 200)
  - `plot_circos`: Circos 플롯 생성 여부 (기본값: FALSE)
  - `output_dir`: 플롯을 저장할 디렉토리 (기본값: NULL)
- **반환값:** ligand_activities, best_ligands, ligand_target_links, ligand_receptor_links, plots, warnings를 포함하는 리스트

#### 함수: `prepare_nichenet_circos_data`
- **설명:** Circos 시각화를 위해 리간드-수용체 및 리간드-타겟 데이터를 준비합니다.
- **파라미터:**
  - `ligand_receptor_links`: 리간드-수용체 쌍을 포함하는 데이터 프레임
  - `ligand_target_links`: 리간드-타겟 조절 링크를 포함하는 데이터 프레임
  - `ligand_activities`: 리간드 활성 점수를 포함하는 데이터 프레임
  - `top_n_ligands`: 포함할 상위 리간드 수 (기본값: 10)
  - `top_n_targets`: 리간드당 상위 타겟 수 (기본값: 20)
- **반환값:** Circos 플롯을 위한 형식화된 데이터를 포함하는 리스트

#### 함수: `draw_nichenet_circos_plot`
- **설명:** 리간드-수용체 및 리간드-타겟 상호작용을 보여주는 Circos 플롯을 생성합니다.
- **파라미터:**
  - `circos_data`: `prepare_nichenet_circos_data`로부터의 리스트
  - `ligand_color`: 리간드 색상 (기본값: "#E41A1C")
  - `receptor_color`: 수용체 색상 (기본값: "#377EB8")
  - `target_color`: 타겟 색상 (기본값: "#4DAF4A")
- **반환값:** NULL (플롯이 현재 장치에 그려짐)

### 파일: `myR/R/analysis/differential_expression/differential_expression.R`

#### 함수: `prepare_pseudobulk_edgeR`
- **설명:** 차등 발현 분석을 위해 단일 세포 카운트를 슈도벌크 수준으로 집계합니다.
- **파라미터:**
  - `object`: Seurat 객체
  - `sample_col`: 샘플/반복을 식별하는 메타데이터 컬럼
  - `group_col`: 그룹화를 위한 메타데이터 컬럼 (예: 조건, 처리)
  - `cluster_col`: 클러스터별 분석을 위한 선택적 클러스터 컬럼 (기본값: NULL)
  - `target_cluster`: `cluster_col`이 지정된 경우, 이 클러스터만 분석 (기본값: NULL)
  - `assay`: 사용할 Assay (기본값: "RNA")
  - `slot`: 사용할 Slot (기본값: "counts")
  - `min_cells`: 포함할 샘플당 최소 세포 수 (기본값: 10)
  - `min_counts`: 유전자당 최소 총 카운트 수 (기본값: 10)
- **반환값:** 카운트, 메타데이터, 클러스터 정보를 포함하는 리스트

#### 함수: `run_pseudobulk_deg`
- **설명:** edgeR 기반 슈도벌크 차등 유전자 발현 분석을 수행합니다.
- **파라미터:**
  - `object`: Seurat 객체 또는 준비된 슈도벌크 리스트
  - `sample_col`: 샘플 식별자 컬럼
  - `group_col`: 그룹 비교 컬럼
  - `comparison`: 비교를 위한 두 요소 벡터 (예: c("Treatment", "Control"))
  - `cluster_col`: 선택적 클러스터 컬럼
  - `target_cluster`: 분석할 특정 클러스터
  - `mode`: 분석 모드: "overall", "per_cluster", "specific_cluster"
  - `assay`: 사용할 Assay (기본값: "RNA")
  - `slot`: 사용할 Slot (기본값: "counts")
  - `min_cells`: 샘플당 최소 세포 수 (기본값: 10)
  - `min_counts`: 유전자당 최소 총 카운트 수 (기본값: 10)
  - `fdr_threshold`: FDR 임계값 (기본값: 0.05)
  - `logfc_threshold`: 로그 배수 변화 임계값 (기본값: 0)
- **반환값:** DE 결과가 있는 데이터 프레임 또는 per_cluster인 경우 데이터 프레임 리스트

#### 함수: `.run_edger_analysis`
- **설명:** 내부 edgeR 분석.
- **파라미터:**
  - `counts`: 카운트
  - `metadata`: 메타데이터
  - `group_col`: 그룹 컬럼
  - `comparison`: 비교
  - `fdr_threshold`: FDR 임계값
  - `logfc_threshold`: 로그 배수 변화 임계값
- **반환값:** DE 결과가 있는 데이터 프레임

#### 함수: `linear_seurat`
- **설명:** 연속형, 범주형, 순서형 예측 변수를 지원하는 회귀 변수에 대한 유전자 발현의 선형 회귀 분석을 수행합니다.
- **파라미터:**
  - `sobj`: Seurat 객체
  - `layer`: 발현 레이어: "counts", "data", "scale.data"
  - `features`: 테스트할 피처 (기본값: "all")
  - `regressor`: 메타데이터의 회귀 변수 이름
  - `regressor.type`: 타입: "continuous", "categorical", "ordinal"
  - `reference.level`: 범주형에 대한 참조 수준
  - `ordinal.method`: 순서형에 대한 방법: "linear", "polynomial", "spline"
  - `link.function`: 연결 함수: "linear", "poisson", "negative.binomial"
  - `effect`: 효과 타입: "fixed", "random"
  - `covariates`: 공변량 컬럼 이름
  - `min.cells`: 유전자를 발현하는 최소 세포 수 (기본값: 10)
  - `return.full`: Seurat 객체를 포함한 전체 결과를 반환
- **반환값:** 회귀 결과가 있는 데이터 프레임

#### 함수: `create_analysis_config`
- **설명:** 복잡한 실험 설계(예: 환자, 약물, 시점이 있는 GeoMx) 전반에 걸쳐 메타데이터 컬럼 이름을 일관되게 관리하기 위한 구성 객체를 생성합니다.
- **파라미터:**
  - `patient`: 환자 ID 컬럼 이름
  - `drug`: 약물/처리 컬럼 이름
  - `timepoint`: 시점(예: pre/post) 컬럼 이름
  - `ck`: 층화 변수(예: CK 상태) 컬럼 이름
  - `response`: 처리 반응 컬럼 이름
  - `aoi`: 관심 영역(AOI) ID 컬럼 이름
- **반환값:** 컬럼 이름을 포함하는 명명된 리스트

#### 함수: `fit_lmm_single_gene`
- **설명:** 복잡한 실험 설계를 가진 단일 유전자에 대한 lmer 모델을 피팅하는 내부 헬퍼 함수입니다.
- **파라미터:**
  - `gene_expr`: 발현 값의 숫자 벡터
  - `metadata`: 메타데이터 데이터 프레임
  - `config`: `create_analysis_config()`로부터의 분석 구성
  - `formula_str`: 선택적 사용자 정의 수식 문자열
  - `formula_components`: 고정, 상호작용, 랜덤 구성 요소를 포함하는 리스트
  - `use_config_names`: 일반 이름을 구성 이름에 매핑할지 여부
- **반환값:** 모델, 효과, anova, 수렴 상태를 포함하는 리스트

#### 함수: `summarize_lmm_results`
- **설명:** 여러 LMM의 결과를 결합하고 조정된 p-값을 계산하는 내부 헬퍼 함수입니다.
- **파라미터:**
  - `lmm_results`: `fit_lmm_single_gene`의 결과 리스트
  - `config`: 분석 구성
- **반환값:** 모든 모델 효과를 요약하는 깔끔한 데이터 프레임

#### 함수: `run_lmm_multiple_genes`
- **설명:** 여러 유전자에 대해 병렬로 LMM을 적용합니다. 복잡한 실험 설계(예: 환자, 약물, 시점, 반응)를 사용한 LMM 분석의 주요 작업 함수입니다.
- **파라미터:**
  - `seurat_obj`: Seurat 객체
  - `genes`: 분석할 유전자 이름의 문자 벡터
  - `config`: `create_analysis_config()`로부터의 분석 구성
  - `formula_str`: 선택적 사용자 정의 수식 문자열
  - `formula_components`: 모델 수식 구성 요소를 지정하는 리스트
  - `use_config_names`: 수식에서 구성 이름을 사용할지 여부
  - `n_cores`: 병렬 처리를 위한 CPU 코어 수
  - `verbose`: 진행 메시지를 출력할지 여부
- **반환값:** raw_results, summary, converged_genes, total_genes를 포함하는 리스트

#### 함수: `find_response_differential_genes`
- **설명:** LMM 결과로부터 약물에 따라 처리 반응이 다른 유전자를 식별합니다.
- **파라미터:**
  - `lmm_summary`: `run_lmm_multiple_genes`의 요약 데이터 프레임
  - `config`: 분석 구성
  - `drug_name`: 집중할 선택적 특정 약물 이름
  - `top_n`: 반환할 상위 유전자 수
- **반환값:** 효과 크기 순으로 정렬된 상위 유전자 데이터 프레임

#### 함수: `find_drug_specific_genes`
- **설명:** 특정 약물과 가장 강하게 연관된 유전자를 식별합니다.
- **파라미터:**
  - `lmm_summary`: `run_lmm_multiple_genes`의 요약 데이터 프레임
  - `config`: 분석 구성
  - `top_n`: 반환할 상위 유전자 수
- **반환값:** 효과 크기 순으로 정렬된 상위 유전자 데이터 프레임
