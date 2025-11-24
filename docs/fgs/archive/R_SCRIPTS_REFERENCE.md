# R 스크립트 파일 참조 가이드

이 문서는 FGS 프로젝트에서 R로 sourcing할 수 있는 주요 파일들을 정리합니다.

## 📁 실행 스크립트 (scripts/fgs/)

### 주요 실행 스크립트

1. **`scripts/fgs/benchmark_l2_methods.R`**
   - 각 L2 방법론별 소요시간 벤치마크 측정
   - 사용법:
     ```r
     cd /home/user3/GJC_KDW_250721
     Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/benchmark_l2_methods.R
     ```
   - 입력: `/data/user3/sobj/fgs/fgs2.qs`, `/data/user3/sobj/data_seurat_251104.qs`
   - 출력: `outputs/benchmark_l2/l2_benchmark_YYYYMMDD_HHMMSS.csv`

2. **`scripts/fgs/run_tml7_is5s_full.R`**
   - IS6 데이터셋 전체 FGS + TML7 파이프라인 실행
   - 사용법:
     ```r
     cd /home/user3/GJC_KDW_250721
     Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/run_tml7_is5s_full.R
     ```
   - 입력: `/data/user3/sobj/IS6_sex_added_0.1x_251110.qs`
   - 출력: `outputs/fgs_is5s/` (fgs, tml, cmgi 결과)

3. **`scripts/fgs/run_is5s_v54.R`** (레거시)
   - IS6 v5.4 실행 스크립트 (이전 버전)
   - `run_tml7_is5s_full.R`로 대체 권장

4. **`scripts/fgs/run_FGS_TML6_analysis.R`**
   - FGS + TML6 분석 실행 스크립트

## 📚 함수 정의 파일 (myR/R/)

### 핵심 함수 파일

1. **`myR/R/signature.R`** ⭐ **가장 중요**
   - `find_gene_signature_v5.4()`: Gene signature 찾기
   - `TML7()`: Meta-learner 학습
   - `compute_meta_gene_importance()`: Gene importance 계산
   - `add_meta_signature_score()`: Signature score를 Seurat에 추가
   - 사용법:
     ```r
     source('/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R')
     ```

2. **`myR/R/test.R`**
   - 테스트/레거시 함수들
   - 일부 중복 함수 포함 (주의)

### 유틸리티 함수 파일

3. **`myR/R/utils_*.R`**
   - `utils_my.R`: 일반 유틸리티
   - `utils_data.R`: 데이터 처리 유틸리티
   - `utils_demulti.R`: Demultiplexing 유틸리티
   - `utils_markers.R`: Marker 관련 유틸리티
   - `utils_aggregation.R`: Aggregation 유틸리티
   - `utils_validation.R`: Validation 유틸리티

### 분석 함수 파일

4. **`myR/R/analysis.R`**
   - 주요 분석 함수들

5. **`myR/R/analysis/*.R`**
   - `pseudobulk_deg.R`: Pseudobulk DEG 분석
   - `milo_pipeline.R`: MILO 분석 파이프라인
   - `nichenet_analysis.R`: NicheNet 분석
   - `trajectory_inference.R`: Trajectory 추론
   - `pathway_enrichment.R`: Pathway enrichment
   - 기타 분석 관련 함수들

### 시각화 함수 파일

6. **`myR/R/plots*.R`**
   - `plots.R`: 기본 plotting 함수
   - `plots_scatter.R`: Scatter plot
   - `plots_volcano.R`: Volcano plot
   - `plots_heatmap.R`: Heatmap
   - `plots_box.R`: Box plot

### 특수 목적 함수 파일

7. **`myR/R/lds*.R`**
   - `lds.R`: LDS (Limma-Dream-SVA) 분석
   - `lds_corrplot.R`: LDS correlation plot
   - `lds_08_heatmaps.R`: LDS heatmap

8. **`myR/R/GeoMx.R`**
   - GeoMx 데이터 처리 함수

9. **`myR/R/CCI.R`**
   - Cell-Cell Interaction 분석

10. **`myR/R/patient_dim_reduction.R`**
    - 환자 수준 차원 축소

## 🔧 환경 설정

모든 스크립트는 다음 환경에서 실행해야 합니다:

```r
# 1. 작업 디렉토리로 이동
cd /home/user3/GJC_KDW_250721

# 2. 환경 초기화 (패키지 로드)
source('start.R')

# 3. 함수 로드 (필요시)
source('/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R')

# 또는 devtools로 패키지 전체 로드
devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE)
```

## 📝 사용 예시

### 예시 1: 함수만 로드하기

```r
cd /home/user3/GJC_KDW_250721
R

# R 세션에서
source('start.R')
source('/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R')

# 함수 사용
fgs_result <- find_gene_signature_v5.4(data = is5s, target_var = "response", ...)
```

### 예시 2: 전체 스크립트 실행

```bash
cd /home/user3/GJC_KDW_250721
Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/benchmark_l2_methods.R
```

### 예시 3: devtools로 패키지 로드

```r
cd /home/user3/GJC_KDW_250721
R

# R 세션에서
source('start.R')
devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE)

# 모든 함수 사용 가능
```

## ⚠️ 주의사항

1. **패키지 충돌**: 반드시 `/home/user3/GJC_KDW_250721`에서 `start.R`을 먼저 실행
2. **환경 변수**: CPU 코어 제한 등 환경 변수는 스크립트 내에서 자동 설정됨
3. **레거시 파일**: `test.R`, `test_*.R` 파일들은 테스트용이므로 프로덕션에서는 사용 지양
4. **중복 함수**: `test.R`에 일부 중복 함수가 있어 `signature.R` 우선 사용 권장

## 📂 파일 구조 요약

```
fgs/
├── scripts/fgs/          # 실행 스크립트
│   ├── benchmark_l2_methods.R
│   ├── run_tml7_is5s_full.R
│   ├── run_is5s_v54.R
│   └── run_FGS_TML6_analysis.R
│
└── myR/R/                # 함수 정의
    ├── signature.R       # ⭐ 핵심 함수
    ├── test.R            # 테스트/레거시
    ├── utils_*.R         # 유틸리티
    ├── analysis.R        # 분석 함수
    ├── plots*.R          # 시각화
    └── analysis/         # 분석 서브모듈
        ├── pseudobulk_deg.R
        ├── milo_pipeline.R
        └── ...
```

## 🔗 관련 문서

- `docs/fgs/README.md`: 전체 프로젝트 개요
- `docs/fgs/TML6_IMPROVEMENTS_CONTEXT.md`: TML6/7 개선 내역

