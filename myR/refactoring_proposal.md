# `myR` 패키지 리팩토링 제안

## 1. 현황 및 필요성

현재 `myR` 패키지는 `R/` 디렉토리 바로 아래와 그 하위 디렉토리에 다수의 R 스크립트가 혼재되어 있습니다. 이 구조는 작동은 하지만, 함수 간의 관계를 파악하고 코드를 탐색하기 어렵게 만듭니다. 특히 하위 디렉토리에 스크립트를 두는 것은 R 패키지의 표준적인 방식이 아니며, `devtools::load_all()` 사용 시 로딩 문제를 일으킵니다.

본 제안서는 함수들을 **맥락(예: scRNA-seq, 공간 전사체)**과 **기능(예: 데이터 준비, 분석, 시각화)**에 따라 그룹화하는 새로운 파일 구조를 제시합니다.

**목표:**
- 코드의 체계성과 유지보수성 향상
- 새로운 사용자가 패키지를 더 쉽게 이해할 수 있도록 개선
- R 패키지 개발 모범 사례 준수
- 모든 함수가 정확하게 로드되도록 보장

## 2. 제안 파일 구조

모든 R 소스 파일은 현재의 하위 디렉토리에서 최상위 `R/` 디렉토리로 이동해야 합니다. 파일 이름은 해당 파일의 내용을 잘 나타내도록 명명합니다.

### 핵심 및 유틸리티 함수 (범용)

이 함수들은 특정 분석 맥락에 종속되지 않으며, 다양한 종류의 분석에서 공통적으로 사용될 수 있습니다.

| 제안 파일명 | 내용 설명 | 포함될 함수 예시 |
|---|---|---|
| `utils-general.R` | 데이터 조작, 문자열 처리 등 기본적인 헬퍼 함수 | `sort_samples`, `convert_gene_ids` |
| `utils-seurat.R` | Seurat 객체와 상호작용하기 위한 헬퍼 함수 | `get_feature_vec`, `prepare_metadata_table`|
| `utils-plotting.R` | 범용 시각화 유틸리티 및 ggplot 테마 정의 | `save_plot`, custom ggplot themes |
| `data-preparation.R`| 분석을 위한 데이터를 준비하는 함수 | `prepare_count_matrix`, `convert_to_long_format` |

### scRNA-seq 분석 함수

표준적인 단일세포 RNA 시퀀싱 분석 워크플로우에 주로 사용되는 함수들입니다.

| 제안 파일명 | 내용 설명 | 포함될 함수 예시 |
|---|---|---|
| `scrnaseq-analysis.R` | DEG, 클러스터링 등 핵심 scRNA-seq 분석 단계 | `run_pseudobulk_deg`, `find_markers` |
| `scrnaseq-markers.R` | 마커 유전자 식별 및 처리 관련 함수 | `process_marker_list` |
| `scrnaseq-pathways.R`| 경로 및 유전자 세트 농축 분석 (GSEA, GO, KEGG) | `myGO`, `run_gsea_analysis` |
| `scrnaseq-signatures.R` | 유전자 시그니처 스코어링 및 발굴 | `AddMultipleModuleScores`, `score_signature` |
| `scrnaseq-trajectory.R` | 의사시간 및 세포궤적 추론 분석 (예: Slingshot) | `run_slingshot_from_seurat` |
| `scrnaseq-cci.R` | 세포 간 상호작용 분석 (예: NicheNet) | `run_nichenet_analysis` |

### 공간 전사체 분석 함수

공간 데이터 분석, 특히 GeoMx에 특화된 함수들입니다.

| 제안 파일명 | 내용 설명 | 포함될 함수 예시 |
|---|---|---|
| `spatial-geomx.R` | GeoMx 데이터 준비, 정규화 및 분석 함수 | `prepare_geomx_data`, `q3_normalize`, `find_deg_geomx` |
| `spatial-analysis.R` | 향후 추가될 수 있는 일반적인 공간 데이터 분석 함수 | (추후 개발을 위한 공간) |

### 시각화 함수

각 분석 맥락에 특화된 시각화 함수들입니다.

| 제안 파일명 | 내용 설명 | 포함될 함수 예시 |
|---|---|---|
| `plots-scrnaseq.R` | scRNA-seq 데이터 시각화 (예: 히트맵, 볼케이노 플롯) | `myhm_genes4`, `plot_volcano` |
| `plots-composition.R` | 클러스터/세포 유형 구성 비율 시각화 | `cmb` (proportional barplot), `plot_cluster_fractions` |
| `plots-spatial.R` | 공간 데이터 시각화 | (조직 이미지 위 플로팅 함수 등) |

### 사용되지 않는(Deprecated) 및 레거시 함수

| 제안 파일명 | 내용 설명 | 포함될 함수 예시 |
|---|---|---|
| `deprecated.R` | 더 이상 사용을 권장하지 않지만 재현성을 위해 남겨둔 함수 | `myhm_genesets2_legacy` |
| `legacy.R` | 리팩토링 이전 버전의 오래된 함수 | `pseudobulk_linear_fit_legacy` |


## 3. 실행 계획

1.  **새 파일 생성:** 제안에 따라 `R/` 디렉토리에 새로운 R 스크립트 파일들을 생성합니다.
2.  **함수 이동:** 기존 파일에 있던 함수 정의를 새로운 파일로 적절하게 이동시킵니다.
3.  **기존 파일 삭제:** 모든 함수가 이동되면, 기존 파일들과 비어있는 `R/` 내의 하위 디렉토리들을 삭제합니다.
4.  **문서 업데이트:** `devtools::document()`를 실행하여 변경 사항을 반영한 `NAMESPACE`와 문서를 재생성합니다.
5.  **테스트:** `devtools::check()`를 실행하여 패키지가 여전히 유효한지, 모든 함수가 정상적으로 익스포트되고 작동하는지 확인합니다.

이 리팩토링을 통해 훨씬 더 깔끔하고, 직관적이며, 견고한 R 패키지를 만들 수 있습니다.
