# FGS (Find Gene Signature) & TML7 (Train Meta-Learner) Documentation

## 개요

이 디렉토리는 FGS (Find Gene Signature) 및 TML7 (Train Meta-Learner) 관련 문서와 스크립트를 관리합니다.

## 디렉토리 구조

```
fgs/
├── scripts/fgs/          # 실행 스크립트
│   ├── benchmark_l2_methods.R    # L2 방법론별 벤치마크
│   ├── run_tml7_is5s_full.R     # IS6 데이터 전체 파이프라인
│   └── run_is5s_v54.R           # IS6 v5.4 실행 (레거시)
└── docs/fgs/            # 문서
    ├── README.md                # 이 파일
    ├── TML6_IMPROVEMENTS_CONTEXT.md  # TML6/7 개선 작업 컨텍스트
    ├── FGS_TML6_ANALYSIS_CONTEXT.md  # FGS/TML6 분석 컨텍스트
    ├── COMPUTE_META_GENE_IMPORTANCE.md  # compute_meta_gene_importance 함수 문서
    ├── CPU_CONFIGURATION.md     # CPU 설정 가이드
    ├── PROGRESS_TRACKING.md     # 진행도 표시 기능
    └── L2_METHOD_TESTING.md     # L2 메서드 테스트 가이드
```

## 주요 기능

### TML7 개선 사항

1. **Group-wise Cross-Validation**
   - `cv_group_var` 파라미터로 환자 단위 누수 방지
   - 기본값: `"emrid"` (GeoMx 데이터) 또는 `"hos_no"` (IS6 데이터)
   - 그룹 메타데이터가 없으면 자동으로 cell-wise CV로 fallback

2. **확장된 L2 Meta-learner 방법론**
   - 기본: `glm`, `ranger`, `xgbTree`
   - 추가: `glmnet`, `svmRadial`, `mlp`, `mlpKerasDropout`, `nnet`, `earth`
   - 패키지 의존성 자동 검사 및 경고

3. **CPU 코어 제한**
   - 최대 16개 코어로 제한하여 과도한 리소스 사용 방지
   - `parallel_workers` 파라미터로 조정 가능

4. **compute_meta_gene_importance 개선**
   - `target_model` 파라미터 추가: 특정 모델의 importance 계산 가능
   - Graceful error handling: importance 추출 실패 시 `NULL` 반환 (경고)
   - 향상된 이름 매칭: `make.names()` 변환된 이름과 원본 이름 자동 매칭
   - 모델별 특화 처리: `ranger`, `earth` 등 모델별 importance 추출 로직 개선
   
   자세한 내용은 [COMPUTE_META_GENE_IMPORTANCE.md](./COMPUTE_META_GENE_IMPORTANCE.md) 참조

## 사용법

### 기본 초기화 (복사-붙여넣기)

```r
# FGS 환경 초기화
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')

# 필요한 패키지 로드
library(qs)
library(pROC)

# FGS 함수 로드
devtools::load_all('/home/user3/data_user3/git_repo/_wt/fgs/myR', quiet = TRUE)
```

### CPU 설정 커스터마이징

```r
# 병렬 처리 활성화 (최대 8코어, BLAS 4스레드)
Sys.setenv(
  FGS_MAX_CPU_CORES = "8",
  FGS_BLAS_THREADS = "4",
  FGS_DISABLE_PARALLEL = "FALSE"
)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

자세한 내용은 [CPU_CONFIGURATION.md](CPU_CONFIGURATION.md) 참조.

### 벤치마크 실행

각 L2 방법론의 소요시간을 측정:

```bash
cd /home/user3/GJC_KDW_250721 && Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/benchmark_l2_methods.R
```

결과는 `outputs/benchmark_l2/` 디렉토리에 저장됩니다.

### 전체 파이프라인 실행

IS6 데이터셋에 대한 전체 FGS + TML7 파이프라인:

```bash
cd /home/user3/GJC_KDW_250721 && Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/run_tml7_is5s_full.R
```

### 진행도 표시

FGS v5.4와 TML7은 자동으로 진행도를 표시합니다:
- 각 메서드별 예상 시간 (이전 실행 기록 기반)
- 실제 실행 시간
- 남은 시간 예측

자세한 내용은 [PROGRESS_TRACKING.md](PROGRESS_TRACKING.md) 참조.

## 데이터셋 정보

### IS6 (is5s)
- 경로: `/data/user3/sobj/IS6_sex_added_0.1x_251110.qs`
- Target: `g3` (1=NR, 2=R) → `response`로 변환
- Group variable: `hos_no` (23명 환자)
- Control: `hos_no` (confounder correction)

### GeoMx (data_seurat)
- 경로: `/data/user3/sobj/data_seurat_251104.qs`
- Target: `response` (R/NR)
- Group variable: `emrid` (환자 ID)
- 샘플 수: 114개 (작아서 테스트에 빠름)

## 주의사항

1. **공선성 문제**: `hos_no < GEM < set` 완전포함 관계이므로 동시 사용 시 주의
2. **CPU 사용량**: 새로 추가된 방법론(특히 `mlpKerasDropout`, `svmRadial`)은 시간이 오래 걸릴 수 있음
3. **패키지 의존성**: 일부 방법론은 추가 패키지 설치 필요 (자동으로 제외됨)

## 참고 문서

- `TML6_IMPROVEMENTS_CONTEXT.md`: TML6/7 개선 작업 상세 내역
- `FGS_TML6_ANALYSIS_CONTEXT.md`: FGS/TML6 분석 컨텍스트

