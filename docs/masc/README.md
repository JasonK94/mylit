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

### CLI 실행 (추천)

`optparse` 기반 CLI: `scripts/masc/run_masc.R`

```bash
Rscript scripts/masc/run_masc.R \
  -i /data/user3/sobj/is2_IS_3_clustered.qs \
  -o /data/user3/sobj/masc/stroke_complex_cli \
  --cluster_var anno3 \
  --contrast_var g3 \
  --random_effects hos_no \
  --fixed_effects GEM,SET,age,sex,bmi,hx_smok,hx_alcohol \
  --prefix masc_anno3_complex
```

### Plot 저장 형식

- 기존: plot bundle을 `.qs`로 저장
- 현재: `--no_plot`이 아니면 **PNG/PDF 파일로 저장** (예: `..._plots_or_forest.png`, `..._plots_or_forest.pdf`)

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


### 변수 선택 및 모델 설정

#### 변수 타입 인식

MASC 파이프라인은 입력 변수들을 자동으로 분류하고 로그를 출력합니다:

```
Variable categories:
  * numeric: age
  * categorical: hos_no, GEM, sex
```

- **numeric**: 연속형 변수 (예: `age`, `bmi`). 강제 numeric 캐스팅 적용.
- **categorical**: 범주형 변수 (factor). 문자형/숫자형 ID는 자동으로 factor로 변환.

#### 중첩 구조 주의

**중요**: 이 데이터셋에서는 변수 간 완전 포함 관계가 있습니다:
- `hos_no` ⊂ `GEM` ⊂ `SET` (완전 중첩)

따라서 **동시에 포함하지 않도록 주의**:
- ✅ 권장: `g3 + age + sex + GEM + (1|hos_no)` (GEM과 hos_no는 동시 사용 가능, 단 GEM이 fixed, hos_no가 random)
- ❌ 피해야 할 조합: `GEM + SET` (SET이 GEM을 완전 포함)
- ⚠️ `bmi`는 결측값이 많아 기본적으로 제외하는 것을 권장합니다.

#### CLI에서 자동 필터링

`scripts/masc/run_masc.R`는 `--fixed_effects`에 `GEM,SET`이 함께 들어오면 자동으로 `SET`을 제거하고 경고를 출력합니다.

### renv 활성화

CLI 스크립트는 기본적으로 `/home/user3/GJC_KDW_250721/renv`를 활성화합니다.
다른 경로를 사용하려면 `--renv` 옵션을 지정:

```bash
Rscript scripts/masc/run_masc.R \
  --renv /path/to/your/renv \
  ...
```

