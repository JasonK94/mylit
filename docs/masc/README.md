# MASC Analysis Pipeline

MASC (Mixed-effects Association testing for Single Cells) 파이프라인은 단일세포 데이터에서 클러스터 풍부도(abundance)와 공변량 간의 연관성을 테스트하는 도구입니다.

## 📌 개요

이 파이프라인은 Seurat 객체를 입력으로 받아, 각 클러스터별로 로지스틱 혼합효과 모델(GLMM)을 적합하고, 특정 조건(예: 질병 상태)에 따른 클러스터 비율 변화를 통계적으로 검정합니다.

- **원본 패키지**: [immunogenomics/masc](https://github.com/immunogenomics/masc)
- **주요 기능**:
  - Seurat 객체 지원
  - 데이터 전처리 및 자동 정제
  - 결과 캐싱 및 재현성 보장
  - 시각화 (OR Forest plot, P-value bar plot)

## 🚀 사용법

### 기본 실행

```r
library(Seurat)
# MASC 함수 로드
source("myR/R/masc.R")

# 파이프라인 실행
results <- run_masc_pipeline(
    seurat_obj = seurat_object,
    cluster_var = "cell_type",      # 클러스터 컬럼
    contrast_var = "condition",     # 비교할 변수 (예: Disease vs Control)
    random_effects = c("donor_id"), # 랜덤 효과 (필수)
    fixed_effects = c("sex", "age"),# 고정 효과 (선택)
    output_dir = "results/masc",
    save = TRUE
)

# 결과 확인
print(results$masc_results)
```

### 주요 파라미터

| 파라미터 | 설명 | 예시 |
|---|---|---|
| `seurat_obj` | 분석할 Seurat 객체 | `sobj` |
| `cluster_var` | 세포 유형이 저장된 메타데이터 컬럼 | `"cell_type"` |
| `contrast_var` | 주요 비교 변수 (Factor여야 함) | `"status"` |
| `random_effects` | 랜덤 효과 변수 (최소 1개 필수) | `c("patient_id")` |
| `fixed_effects` | 고정 효과 변수 (공변량) | `c("batch", "sex")` |
| `adjust_pvalue` | FDR 보정 여부 | `TRUE` |

## 🛠 개발 내역

### 2025-12-08: 초기 구현 및 안정화
- **핵심 함수 구현**: `run_masc_pipeline`, `.masc_run_analysis` 등
- **데이터 처리 강화**:
  - `hos_no` 등 숫자형 ID를 문자열로 자동 변환하여 모델 오류 방지
  - `hx_alcohol` 등 긴 문자열 변수를 단순화하는 전처리 로직 추가
  - `cli` 패키지 의존성 제거 및 표준 `cat`/`warning` 메시지로 전환
- **에러 핸들링**:
  - `glmer` 수렴 실패 시 경고 처리 및 진행
  - 최소 샘플 수 부족 클러스터 자동 필터링
  - `model.matrix` 생성 시 컬럼명 충돌 방지

### 2025-12-08: Stroke 데이터 적용
- `g3` 변수에 따른 `anno3big`, `anno3`, `anno.mo` 클러스터 비율 분석 완료.
- T 세포(`Tc`) 클러스터에서 유의미한 비율 증가 확인 (FDR < 0.05).

## 📁 파일 구조

```
myR/R/
└── masc.R              # 핵심 구현 파일

scripts/
├── run_masc_stroke.R   # Stroke 데이터 분석 스크립트
└── run_masc_simple.R   # 단순화된 테스트 스크립트
```

