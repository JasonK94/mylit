# 다음 단계 가이드

## 완료된 작업 ✅

1. 워크트리 및 브랜치 생성
2. 컨텍스트 및 개발 계획 문서화
3. `run_monocle3_from_seurat()` 함수 개발
4. 데이터 정보 문서화 (feature, metadata)
5. 데이터 전처리 함수 개발 (`preprocess_pseudotime_data()`)
6. 테스트 스크립트 작성

## 즉시 실행할 작업

### 1. 데이터 전처리 실행

작업 디렉토리(`/home/user3/GJC_KDW_250721`)에서 실행:

```r
# 전처리 스크립트 실행
source("/home/user3/data_user3/git_repo/_wt/pseudotime/test_preprocessing.R")
```

또는 직접 함수 호출:

```r
source("/home/user3/data_user3/git_repo/_wt/pseudotime/preprocess_data.R")

# 다운샘플 데이터 전처리
sobj_ds <- preprocess_pseudotime_data(
  input_file = "/data/user3/sobj/IS_scvi_251107_ds2500.qs",
  output_file = "/data/user3/sobj/IS_scvi_251107_ds2500_preprocessed.qs"
)

# 원본 데이터 전처리
sobj_full <- preprocess_pseudotime_data(
  input_file = "/data/user3/sobj/IS_scvi_251107.qs",
  output_file = "/data/user3/sobj/IS_scvi_251107_preprocessed.qs"
)
```

**중요**: 전처리된 파일이 생성되면, 이후 모든 분석에서 이 파일을 사용하세요!

### 2. 기본 Pseudotime 분석 테스트

```r
# 기본 분석 테스트
source("/home/user3/data_user3/git_repo/_wt/pseudotime/test_pseudotime_basic.R")
```

이 스크립트는:
- 전처리된 데이터 로드
- 주요 feature 및 metadata 확인
- Slingshot trajectory inference 테스트
- Monocle3 trajectory inference 테스트
- 결과 저장

## 다음 개발 단계

### Phase 2: 고급 패턴 분석 함수 개발

1. **복잡한 패턴 탐지 함수**
   - `detect_expression_patterns()`
   - 증가/감소, oscillatory, bimodal 패턴

2. **Metadata 상관관계 동적 분석**
   - `analyze_metadata_correlation_dynamics()`
   - Pseudotime 구간별 상관관계 변화

3. **Branching Analysis**
   - `analyze_branching_dynamics()`
   - Branching point 탐지 및 분석

### Phase 3: 통합 워크플로우

1. **End-to-End 파이프라인**
   - `run_pseudotime_analysis_pipeline()`

2. **Batch Analysis**
   - `batch_pseudotime_analysis()`

## 파일 구조

```
_wt/pseudotime/
├── PSEUDOTIME_CONTEXT.md          # 작업 환경 및 데이터 정보
├── PSEUDOTIME_DEVELOPMENT_PLAN.md  # 상세 개발 계획
├── PSEUDOTIME_DEVELOPMENT_SUMMARY.md # 개발 요약
├── preprocess_data.R                # 데이터 전처리 함수
├── test_preprocessing.R             # 전처리 테스트
├── test_pseudotime_basic.R         # 기본 분석 테스트
└── myR/
    └── R/
        ├── pseudotime.R             # 주요 pseudotime 함수들
        └── preprocess_pseudotime_data.R # 전처리 함수 (패키지용)
```

## 참고 문서

- `PSEUDOTIME_CONTEXT.md`: 작업 환경, 데이터 정보, 개발 철학
- `PSEUDOTIME_DEVELOPMENT_PLAN.md`: 단계별 개발 계획
- `PSEUDOTIME_DEVELOPMENT_SUMMARY.md`: 완료된 작업 및 다음 단계

