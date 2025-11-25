# LDS 파이프라인 개발 로그

## 2025-01-XX: LDS 워크트리 생성 및 초기 테스트

### 목표
- Limma-Dream-SVA (LDS) 파이프라인 테스트 및 문서화
- 워크트리 기반 개발 환경 구축
- 표준화된 문서 구조 적용

### 완료된 작업

#### 1. 워크트리 및 브랜치 생성
- `../_wt/lds` 워크트리 생성
- `lds` 브랜치 생성 (main에서 분기)
- `docs/lds/`, `scripts/lds/` 디렉터리 구조 생성

#### 2. 문서 작성
- **`lds.md`**: LDS 파이프라인 상세 문서
  - 개요 및 핵심 개념 설명
  - SVA와 limma-dream 설명
  - 파라미터 및 사용 예시
  - 문제 해결 가이드
  
- **`TEST_INSTRUCTIONS.md`**: 테스트 실행 가이드
  - 환경 설정 방법
  - 데이터 준비 및 확인
  - 다양한 테스트 시나리오
  - 결과 분석 방법

- **`devlog.md`**: 개발 로그 (본 문서)

#### 3. 테스트 스크립트 작성
- **`scripts/lds/test_lds.R`**: 자동화된 테스트 스크립트
  - 데이터 로드 및 검증
  - 기본 LDS 테스트
  - SV 개수 지정 테스트
  - 결과 저장

#### 4. 표준화된 문서 구조 제안
다른 워크트리들과 일관된 구조:
```
docs/
  lds/
    lds.md              # 메인 문서 (기능 설명)
    TEST_INSTRUCTIONS.md # 테스트 가이드
    devlog.md           # 개발 로그
scripts/
  lds/
    test_lds.R          # 테스트 스크립트
```

### LDS 함수 구조 분석

#### 함수 위치
- `myR/R/test_working.R` (lines 1554-1757)

#### 주요 단계
1. 입력 데이터 검증 및 추출
2. 포뮬러 파싱 (고정 효과 추출)
3. DGEList 생성, 필터링, 정규화
4. SVA 실행 (voom 변환 후)
5. 최종 포뮬러 생성 (원본 + SV)
6. limma-dream 파이프라인 실행
7. 결과 반환

#### 의존성
- `limma` (dream, voomWithDreamWeights 포함)
- `edgeR` (DGEList, filterByExpr, calcNormFactors, voom)
- `lme4` (포뮬러 파싱)
- `BiocParallel` (병렬 처리)
- `sva` (SVA 실행)

### 다음 단계

1. **실제 데이터로 테스트 실행**
   - is5s 데이터로 기본 테스트
   - is5 전체 데이터로 확장 테스트
   - 다양한 포뮬러 테스트

2. **결과 검증**
   - topTable 결과 확인
   - SV 개수 및 분산 설명 비율 확인
   - 다른 방법과 비교 (MUSCAT, NEBULA)

3. **문서 보완**
   - 실제 테스트 결과 반영
   - 문제 해결 사례 추가
   - 성능 벤치마크 (선택적)

4. **표준화 제안**
   - 다른 워크트리들과 문서 구조 통일
   - 공통 템플릿 제안

## 2025-01-XX: LDS 함수 모듈화 및 개선

### 완료된 작업

#### 1. 함수 모듈화
- `LDS()` 함수를 8개의 독립적인 단계별 함수로 분리:
  - `lds_01_extract_data()`: 데이터 추출
  - `lds_01b_filter_na()`: NA 필터링
  - `lds_02_parse_formula()`: 포뮬러 파싱
  - `lds_03_preprocess_dge()`: DGEList 전처리
  - `lds_04_run_sva()`: SVA 실행
  - `lds_05_build_final_formula()`: 최종 포뮬러 생성
  - `lds_06_run_dream()`: Dream 실행
  - `lds_07_analyze_sva_correlation()`: SVA 상관관계 분석
  - `lds_08_create_heatmaps()`: Heatmap 생성

#### 2. 주요 개선사항

**2-1. NA 필터링 추가**
- `remove_na` 파라미터로 포뮬러 변수의 NA 값 필터링
- `counts_matrix`와 `meta.data`의 차원 일치 보장

**2-2. 필터링 유연성 향상**
- `filterByExpr` 옵션을 파라미터로 노출:
  - `min.count`, `min.total.count`, `min.prop`, `large.n`
- 기본값 완화 (scRNA-seq 데이터에 적합)

**2-3. p-value 계산 수정**
- `dream` 결과에서 t-statistics와 df를 사용하여 p-value 직접 계산
- `topTable()`에서 p-value 정상 출력 확인

**2-4. SV 개수 결정 개선**
- `sv_var_cutoff = NULL`이면 전체 유의미한 SV 사용
- 잔차 분산 기반 자동 결정 또는 사용자 지정

**2-5. SVA 상관관계 분석 및 시각화**
- `lds_07_analyze_sva_correlation()`: 상관관계 행렬 계산
- `lds_08_create_heatmaps()`: 3가지 Heatmap 자동 생성
  - 전체 상관관계 행렬 (metadata + SV)^2
  - 메타데이터 × SV
  - 상위 n개 메타데이터 × SV (p<0.05 기준)

**2-6. Helper 함수 추가**
- `lds_corrplot()`: 독립적으로 사용 가능한 상관관계 플롯 함수
- `lds_corrplot_summary()`: 상관관계 요약 테이블 생성
- SV가 열이 되도록 재배치

#### 3. 파일 구조

**모듈화된 함수들**:
- `myR/R/lds.R`: 메인 LDS 함수 및 단계별 함수
- `myR/R/lds_corrplot.R`: Helper 함수 (독립 사용 가능)
- `myR/R/lds_08_heatmaps.R`: Step 8 Heatmap 생성 함수

**차이점**:
- `lds_corrplot()`: 독립적으로 사용 가능한 helper 함수 (파이프라인 외부)
- `lds_08_create_heatmaps()`: LDS 파이프라인의 일부 (Step 8에서 자동 실행)

#### 4. SVA 원리 정리

**SVA가 설명하는 분산**:
- 전체 유전자의 총 분산이 아니라, **잔차 분산(residual variance)**을 설명
- 고정 효과 모델(`mod`)로 설명되지 않는 잔차에서 숨겨진 변동성 탐지
- 잔차 행렬의 SVD를 통해 주성분(SV) 추출

**유의미한 SV의 기준**:
- SVA 패키지의 `sva()` 함수가 permutation test와 통계적 검정으로 자동 결정
- 일반적으로 p-value < 0.05 기준
- 이후 `sv_var_cutoff`로 추가 필터링

#### 5. 테스트 결과

- p-value 정상 생성 확인 (291/291 유전자)
- 필터링 유연성 확인 (291개 유전자 통과, 0.6%)
- 3가지 Heatmap 자동 생성 확인
- `sv_var_cutoff = NULL` 시 전체 SV 사용 확인

## 표준화 제안

### 문서 구조 표준화

각 워크트리에서 다음 구조를 유지:

```
docs/
  <worktree_name>/
    <worktree_name>.md           # 메인 기능 문서
    TEST_INSTRUCTIONS.md          # 테스트 가이드
    devlog.md                     # 개발 로그
scripts/
  <worktree_name>/
    test_<worktree_name>.R        # 테스트 스크립트
```

### 문서 작성 가이드

1. **메인 문서 (`<name>.md`)**
   - 개요 및 핵심 개념
   - 주요 함수 설명
   - 파라미터 및 사용 예시
   - 주의사항 및 문제 해결

2. **테스트 가이드 (`TEST_INSTRUCTIONS.md`)**
   - 환경 설정
   - 데이터 준비
   - 테스트 실행 방법
   - 결과 확인 및 저장

3. **개발 로그 (`devlog.md`)**
   - 날짜별 작업 기록
   - 완료된 작업
   - 다음 단계
   - 이슈 및 해결 방법

### 공통 패턴

- **데이터 경로**: `/data/user3/sobj/` 사용
- **결과 저장**: `.qs` 형식 (큰 객체), `.csv` (표)
- **R 세션**: `/home/user3/GJC_KDW_250721`에서 시작
- **패키지 로드**: `devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")`

