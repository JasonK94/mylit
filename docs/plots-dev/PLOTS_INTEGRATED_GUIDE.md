# Plotting Module Integrated Guide

이 문서는 Plotting 모듈의 통합 가이드입니다. (작성 중)

## 1. Introduction (소개)

### 목적
Seurat 객체 및 데이터프레임 기반의 다양한 시각화 함수를 제공하여 분석 결과를 효과적으로 표현합니다.

## 2. Functions (주요 함수)

현재 개발/리팩토링 중인 주요 함수들은 다음과 같습니다 (상세 분석은 `archive/function_issues.md` 참조):

*   **upset_gene_lists**: 유전자 리스트 간의 교집합 시각화 (Upset plot).
*   **vln_p**: Violin plot과 통계적 유의성(p-value) 표시.
*   **cmb**: Proportional Bar Graph (비율 막대 그래프).
*   **acmb**: Absolute Count Bar Graph (절대 수 막대 그래프).
*   **cml**: Cumulative Line Graph (누적 선 그래프).
*   **cdf**: Cumulative Distribution Function plot.
*   **cdf_multi**: 다중 변수 CDF plot.

## 3. Development Status (개발 현황)

*   입력 유연성(Seurat/data.frame 모두 지원) 확보 및 패키지 의존성 문제 해결을 위한 리팩토링이 진행 중입니다.

