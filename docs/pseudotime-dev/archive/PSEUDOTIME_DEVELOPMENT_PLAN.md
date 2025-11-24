# Pseudotime Analysis 개발 계획

## 데이터 정보

### 주요 Feature (유전자)
분석 시 다음 유전자들을 중시해야 함:
- `DDIT4`, `UTY`, `S100B`, `XIST`, `HLA-B`, `CCL4`, `HLA-C`, `TXNIP`

### Metadata Columns (매우 중요)
1. **nih_change** (numeric, 매우 중요): NIH 점수 변화량
2. **sex** (group, 매우 중요): M, F로 변환 필요 (전처리 스크립트 사용)
3. **g3** (group, 매우 중요): 1, 2, NA - factor로 변환 필요
4. **arrival_gcs_score** (numeric): 도착 시 GCS 점수
5. **age** (numeric): 나이
6. **ini_nih** (numeric): 초기 NIH 점수
7. **nih_end1** (numeric): 종료 시점 NIH 점수
8. **시간 변수들** (숫자 형식, 파싱 필요):
   - `icu_adm_dt`: ICU 입원 시간
   - `ia_start`: 시작 시간
   - `ia_end`: 종료 시간
   - 시간 차이 변수들도 자동 계산됨

**중요**: `preprocess_pseudotime_data()` 함수로 전처리 후 저장된 파일을 사용하세요!

## 현재 상태 분석

### 기존 구현된 기능
1. **Slingshot 통합** (`run_slingshot_from_seurat`)
   - Seurat 객체에서 Slingshot trajectory inference
   - SingleCellExperiment 객체로 결과 반환
   - 잘 구현되어 있음

2. **GAM 기반 동역학 분석** (`analyze_gene_dynamics`)
   - Monocle3 cell_data_set 객체 사용
   - 조건별 gene expression 패턴 분석
   - TV (Total Variation), DR (Dynamic Range) 계산
   - Interaction p-value 계산

3. **tradeSeq 통합** (`analyze_gene_dynamics_tradeSeq`)
   - Slingshot 결과를 tradeSeq로 분석
   - patternTest, conditionTest 등 지원

### 개선 및 추가 필요 사항

## 개발 계획

### Phase 1: 기본 기능 강화 (우선순위: 높음)

#### 1.1 Monocle3 Trajectory Inference 함수 추가
- **함수명**: `run_monocle3_from_seurat()`
- **기능**:
  - Seurat 객체를 Monocle3 cell_data_set로 변환
  - Preprocessing, dimensionality reduction, trajectory inference 자동화
  - Root cell 선택 옵션
  - Multiple trajectory 지원
- **출력**: Monocle3 cell_data_set 객체

#### 1.2 통합 Trajectory Inference 래퍼
- **함수명**: `run_trajectory_inference()`
- **기능**:
  - Slingshot, Monocle3 중 선택 가능
  - 자동 최적화 옵션
  - 결과 비교 및 검증
- **출력**: 통합된 결과 객체

### Phase 2: 고급 패턴 분석 (우선순위: 중간)

#### 2.1 복잡한 패턴 탐지
- **함수명**: `detect_expression_patterns()`
- **기능**:
  - 증가/감소 패턴
  - Oscillatory 패턴 (주기적 변화)
  - Bimodal 패턴
  - Branching point에서의 패턴 변화
- **방법론**:
  - GAM 기반 패턴 분류
  - Fourier transform 기반 주기 탐지
  - Change point detection

#### 2.2 Metadata 상관관계 동적 분석
- **함수명**: `analyze_metadata_correlation_dynamics()`
- **기능**:
  - Pseudotime 구간별 metadata와 gene expression 상관관계 계산
  - 상관관계 변화 패턴 탐지
  - 조건별 상관관계 비교
- **출력**: 구간별 상관계수 및 변화 통계

#### 2.3 Branching Analysis
- **함수명**: `analyze_branching_dynamics()`
- **기능**:
  - Branching point 탐지
  - Branch별 gene expression 차이 분석
  - Branch 선택 메커니즘 분석
- **방법론**: Monocle3 branching analysis + 통계 검정

### Phase 3: 통합 워크플로우 (우선순위: 중간)

#### 3.1 End-to-End 워크플로우
- **함수명**: `run_pseudotime_analysis_pipeline()`
- **기능**:
  - Seurat 객체 입력
  - Trajectory inference 자동 실행
  - Gene list에 대한 동역학 분석
  - 결과 통합 및 시각화
  - 결과 저장
- **파라미터**:
  - Trajectory method 선택
  - 분석할 유전자 리스트
  - 조건 변수
  - 출력 디렉토리

#### 3.2 Batch Analysis
- **함수명**: `batch_pseudotime_analysis()`
- **기능**:
  - 여러 조건/샘플에 대한 일괄 분석
  - 결과 비교 및 통합
  - 메타 분석

### Phase 4: 시각화 및 보고서 (우선순위: 낮음)

#### 4.1 고급 시각화
- Trajectory overlay plots
- Heatmap of gene dynamics
- Correlation heatmap across pseudotime
- Branching visualization

#### 4.2 자동 보고서 생성
- HTML 보고서 생성
- 주요 발견사항 요약
- 통계 테이블 및 플롯 포함

## 구현 순서

### Week 1: 기본 기능 강화
1. `run_monocle3_from_seurat()` 구현 및 테스트
2. `run_trajectory_inference()` 래퍼 구현
3. 기존 함수들 개선 및 버그 수정

### Week 2: 고급 분석 함수
1. `detect_expression_patterns()` 구현
2. `analyze_metadata_correlation_dynamics()` 구현
3. `analyze_branching_dynamics()` 구현

### Week 3: 통합 워크플로우
1. `run_pseudotime_analysis_pipeline()` 구현
2. `batch_pseudotime_analysis()` 구현
3. 통합 테스트

### Week 4: 시각화 및 문서화
1. 시각화 함수 개선
2. 문서 업데이트 (한글 우선)
3. 사용 예제 작성

## 테스트 전략

### 1단계: 다운샘플 데이터 (`IS_scvi_251107_ds2500.qs`)
- 각 함수별 단위 테스트
- 통합 워크플로우 테스트
- 오류 디버깅
- 결과 저장 확인

### 2단계: 원본 데이터 (`IS_scvi_251107.qs`)
- 동일한 테스트 재실행
- 성능 검증
- 결과 품질 검증

## 기술 스택

### 핵심 패키지
- **Monocle3**: Trajectory inference
- **Slingshot**: Trajectory inference
- **tradeSeq**: Differential expression along trajectories
- **mgcv**: GAM modeling
- **Seurat**: Single-cell analysis
- **SingleCellExperiment**: Data structure

### 추가 패키지 (필요시)
- **changepoint**: Change point detection
- **fourier**: Fourier analysis for oscillatory patterns
- **ggplot2**: Visualization
- **dplyr/tidyr**: Data manipulation

## 성공 기준

1. ✅ 모든 함수가 다운샘플 데이터에서 정상 작동
2. ✅ 원본 데이터에서도 정상 작동
3. ✅ 결과가 `/data/user3/sobj/`에 정상 저장
4. ✅ 문서가 한글 우선으로 업데이트됨
5. ✅ 코드가 `myR` 패키지에 통합됨

