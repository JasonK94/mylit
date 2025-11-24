# Function Issues Analysis

기존 함수들의 문제점 분석 (수정하지 않고 지적만)

## 1. upset_gene_lists

**위치**: `plots.R:795-836`

**문제점**:
- `library(ComplexUpset)` 사용: 함수 내부에서 `library()` 호출은 패키지 로딩을 강제함
- `requireNamespace()` 체크 후에도 `library()` 호출로 인해 namespace 충돌 가능
- 반환값이 `NULL`일 때 명확한 에러 메시지 없음
- `stopifnot()` 사용: 에러 메시지가 불명확함
- Seurat/df 입력 지원 없음 (gene_lists만 받음)

**개선 제안**:
- `requireNamespace()` 체크 후 `ComplexUpset::upset()` 직접 호출
- `stopifnot()` 대신 명확한 에러 메시지
- 반환값이 `NULL`일 때 `stop()` 또는 경고 후 빈 플롯 반환

---

## 2. vln_p

**위치**: `plots.R:838-865`

**문제점**:
- Seurat 객체만 지원 (data.frame 미지원)
- `aes_string()` 사용: deprecated된 함수 (ggplot2 3.4.0+)
- `split.by`가 NULL일 때 처리 없음
- 통계 검정이 항상 `wilcox.test`로 고정
- `pt.size=0`일 때도 통계 검정 수행 (데이터 없을 수 있음)
- 반환값이 `NULL`일 때 처리 없음

**개선 제안**:
- `aes_string()` → `.data[[]]` 또는 `aes()` 사용
- `split.by` NULL 처리 추가
- 통계 검정 방법을 파라미터로 받기
- Seurat/df 모두 지원하도록 확장

---

## 3. cmb (Proportional Bar Graph)

**위치**: `plots.R:902-946`

**문제점**:
- Seurat 객체만 지원
- `sort_samples()` 함수 의존성: 정의되지 않은 함수 사용 가능성
- `identity` 파라미터로 `Idents()` 설정: 부작용(side effect) 발생
- `df=TRUE`일 때도 플롯 생성 로직 실행 (비효율)
- `vlines` 처리 시 `line_pos + 0.5` 하드코딩: 이유 불명확
- `scale_fill_viridis_d()` 사용: 팔레트 선택 불가

**개선 제안**:
- Seurat/df 모두 지원
- `sort_samples()` 함수 정의 또는 제거
- `Idents()` 설정을 선택적으로 (또는 복원)
- `df=TRUE`일 때 플롯 생성 스킵
- `vlines` 로직 명확화
- 팔레트 파라미터 추가

---

## 4. acmb (Absolute Count Bar Graph)

**위치**: `plots.R:983-1023`

**문제점**:
- `cmb`와 동일한 문제점들
- `sort_samples()` 함수 의존성
- `Idents()` 부작용
- Seurat 객체만 지원

**개선 제안**:
- `cmb`와 동일한 개선사항

---

## 5. cml (Cumulative Line Graph)

**위치**: `plots.R:1054-1137`

**문제점**:
- Seurat 객체만 지원
- `sort.by` 파라미터: "name" 또는 "frequency"만 지원 (확장성 부족)
- `color_palette`가 NULL일 때 `RColorBrewer::brewer.pal()` 사용: 패키지 의존성
- `n_patterns`, `n_shapes` 하드코딩된 값들
- `df=TRUE`일 때도 플롯 스타일 계산 수행 (비효율)

**개선 제안**:
- Seurat/df 모두 지원
- `sort.by`에 custom function 지원
- 팔레트 선택을 더 유연하게
- `df=TRUE`일 때 스타일 계산 스킵

---

## 6. cdf (Cumulative Distribution Function)

**위치**: `plots.R:1184-1245`

**문제점**:
- `probability_col`, `ratio_col`을 bare name으로 받지만 실제로는 문자열로 처리 필요
- `plot_type` 검증이 있지만, 입력 데이터에 해당 컬럼이 없을 때 에러 메시지 불명확
- `output_file` 저장 시 디렉토리 생성 로직 없음
- `scales::number_format()`, `scales::scientific()` 사용: 패키지 의존성 명시 필요
- 단일 변수에 대한 CDF만 지원 (비교 불가)

**개선 제안**:
- `enquo()` 사용하여 bare name 지원 명확화
- 컬럼 존재 여부 사전 검증
- 디렉토리 자동 생성
- `@importFrom` 명시
- 다중 변수 비교 옵션 추가

---

## 7. cdf_multi

**위치**: `plots.R:1295-1420`

**문제점**:
- `data_list`가 단일 data.frame일 때와 list일 때 처리 로직이 복잡하고 혼란스러움
- `group_by_col`과 `sample_col`의 역할이 명확하지 않음
- `purrr::map2_dfr()` 사용: 패키지 의존성
- `sample_col`이 각 df에 있어야 한다는 가정이 강함
- 에러 메시지가 상황별로 다름 (일관성 부족)

**개선 제안**:
- 입력 타입 처리 로직 단순화
- `group_by_col`과 `sample_col` 역할 명확화 (문서화)
- `purrr` 의존성 명시 또는 base R로 대체
- `sample_col` 선택적 처리
- 일관된 에러 메시지

---

## 공통 문제점

1. **입력 유연성 부족**: 대부분 Seurat 객체만 지원
2. **에러 처리**: `stopifnot()` 사용으로 에러 메시지 불명확
3. **부작용**: `Idents()` 설정 등으로 객체 상태 변경
4. **패키지 의존성**: `library()` 호출 또는 명시적 의존성 부족
5. **하드코딩**: 팔레트, 테마, 파라미터 값들이 고정
6. **비효율**: `df=TRUE`일 때도 불필요한 계산 수행
7. **일관성**: 함수 간 파라미터 네이밍 불일치

---

## 개선 우선순위

1. **높음**: 입력 유연성 (Seurat/df 지원), 에러 처리 개선
2. **중간**: 부작용 제거, 패키지 의존성 명시
3. **낮음**: 하드코딩 제거, 성능 최적화

