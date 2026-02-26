# Plots Development Module Documentation

이 폴더는 시각화(Plotting) 관련 함수들의 문서를 포함합니다.

## 프로젝트 의의
- `myR` 전체 분석 모듈에서 필요한 공통 시각화를 재사용 가능한 형태로 묶어, scRNAseq/GeoMx 간 워크플로우의 마지막 단계를 책임집니다.
- (docs-main 컨텍스트) 통합 가이드가 요구하는 “데이터 준비 → 통계 분석 → 시각화” 흐름을 끊김 없이 이어 주며, 단일 worktree 환경에서 안전하게 실험/배포할 수 있게 해 줍니다.
- scatter, volcano, heatmap뿐 아니라 환자 빈도/UMAP 등 QC·EDA 필수 플롯을 하나의 네임스페이스로 통합하여, 각 모듈이 독자적으로 플롯을 관리할 때 발생하던 중복을 제거합니다.

## 최근 업데이트 (git log 요약)
- `c8afff6`: `cmb`/`acmb` 함수에 정렬·스택 순서 제어가 추가되어 집계 플롯 스타일을 세밀하게 조절할 수 있게 되었습니다.
- `0c8c5fd`: boxplot 축/레이블 버그 수정으로 QC용 플롯을 바로 리포트에 사용할 수 있습니다.
- `3efafdd`: `plot_box` 리팩터와 `.calculate_sort_values()` 추가로 `sort.by`, `split.by` 조합이 대형 메타데이터에서도 안정화되었습니다.
- (현재 PR) `plot_umap_density()` 도입으로 UMAP 상의 커널 밀도/등고선 비교, 그룹 차이 시각화 기능이 패키지 표준으로 편입되었습니다.

## 메인 가이드
*   [**PLOTS_INTEGRATED_GUIDE.md**](./PLOTS_INTEGRATED_GUIDE.md): 통합 가이드입니다. 주요 플로팅 함수 목록과 사용법(예정)을 설명합니다.
    - UMAP 밀도 플롯(`plot_umap_density`) 예제도 이 가이드에 포함해 전체 플롯 카탈로그를 한 번에 훑을 수 있습니다.

## 아카이브
과거 문서 및 개발 로그는 `archive/` 폴더에서 확인할 수 있습니다.
*   [function_issues.md](./archive/function_issues.md): 기존 함수의 문제점 분석
*   [TESTING.md](./archive/TESTING.md): 테스트 기록
*   [devlog.md](./archive/devlog.md): 개발 로그

