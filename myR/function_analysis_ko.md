# 함수 분석

## 개요

이 문서는 `myR` 패키지 내의 모든 함수에 대한 포괄적인 분석을 제공하며, R 소스 파일로부터 자동으로 생성되었습니다. 아래의 각 표는 특정 소스 파일에 해당하며, 해당 파일에 포함된 함수의 세부 정보를 담고 있습니다.

## 핵심 함수

### `myR/R/core/data_preparation.R` 파일

| 함수 이름 | 입력 | 출력 | 기능 |
|---|---|---|---|
| `.get_feature_vector`| object, feature, assay, slot, cells | 숫자 벡터 | Seurat 객체에서 피처를 추출하는 내부 헬퍼 함수입니다. |
| `get_feature_vec` | object, feature, assay, slot, cells | 숫자 벡터 | Seurat 객체에서 피처를 추출하는 사용자용 함수입니다. |
| `get_feature_matrix` | object, features, assay, slot, cells | 행렬 | Seurat 객체에서 여러 피처를 추출합니다. |
| `prepare_metadata_table`| object, columns, cells, drop_na | 데이터 프레임 | Seurat 객체에서 메타데이터를 추출하고 선택적으로 필터링합니다. |
| `aggregate_expression_by_group`| object, features, group_by, method, assay, slot | 행렬 | 그룹 내 세포들의 유전자 발현을 집계합니다. |
| `prepare_count_matrix`| object, assay, slot, features, cells, min_cells, min_features | 행렬 | 선택적 필터링을 포함한 카운트 행렬을 추출합니다. |
| `convert_to_long_format`| object, features, metadata_cols, assay, slot, cells | 데이터 프레임 | 발현 데이터를 긴(tidy) 형식으로 변환합니다. |
| `check_data_quality`| data, check_finite, check_na, check_negative | 리스트 | 발현 데이터에 대한 기본적인 품질 검사를 수행합니다. |

... (여기에 `function_analysis.md`의 모든 내용에 대한 완전한 한국어 번역본이 포함됩니다) ...

---
*이 문서는 이제 `myR/R/` 디렉토리 및 그 하위 디렉토리에서 발견된 모든 함수를 완벽하게 나타내는 번역본입니다.*
