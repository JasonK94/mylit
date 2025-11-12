# Pseudotime Analysis 개발 컨텍스트

## 작업 환경

### 워크트리 및 브랜치
- **워크트리 경로**: `/home/user3/data_user3/git_repo/_wt/pseudotime`
- **브랜치**: `pseudotime-dev`
- **메인 레포지토리**: `/home/user3/data_user3/git_repo/mylit`
- **목적**: 다른 브랜치들과의 충돌 방지

### 작업 디렉토리
- **프로젝트 루트**: `/home/user3/data_user3/git_repo/mylit`
- **작업 디렉토리**: `/home/user3/GJC_KDW_250721`
- **R 세션 초기화**: `st/start.R` 자동 실행
  - 작업 디렉토리에서 R을 실행하면 모든 라이브러리가 정상적으로 로드됨

### 테스트 데이터
- **다운샘플링 데이터**: `/data/user3/sobj/IS_scvi_251107_ds2500.qs`
  - 테스트 및 디버깅용
- **원본 데이터**: `/data/user3/sobj/IS_scvi_251107.qs`
  - 최종 검증용
- **전처리된 데이터**: `/data/user3/sobj/IS_scvi_251107_preprocessed.qs`
  - sex, time 변수가 전처리된 버전 (권장 사용)

### 데이터 Feature 및 Metadata

#### 주요 Feature (유전자)
다음 유전자들을 분석 시 중시해야 함:
- `DDIT4`, `UTY`, `S100B`, `XIST`, `HLA-B`, `CCL4`, `HLA-C`, `TXNIP`

#### Metadata Columns (매우 중요)
1. **nih_change** (numeric, 매우 중요)
   - NIH 점수 변화량

2. **sex** (group, 매우 중요)
   - 원본 데이터의 encoding이 이상할 수 있음
   - **반드시 M, F로 변환 필요**
   - 전처리 스크립트로 변환 후 저장된 파일 사용 권장

3. **g3** (group, 매우 중요)
   - 값: 1, 2, NA
   - numeric으로 인식될 수 있으므로 주의 필요
   - factor로 변환 권장

4. **arrival_gcs_score** (numeric)
   - 도착 시 GCS 점수

5. **age** (numeric)
   - 나이

6. **ini_nih** (numeric)
   - 초기 NIH 점수

7. **nih_end1** (numeric)
   - 종료 시점 NIH 점수

8. **시간 변수들** (숫자 형식, 파싱 필요할 수 있음)
   - **icu_adm_dt**: ICU 입원 시간
   - **ia_start**: 시작 시간
   - **ia_end**: 종료 시간
   - 이들 간의 차이 계산 가능
   - **전처리 스크립트로 파싱 후 저장된 파일 사용 권장**

### 결과 저장
- **저장 경로**: `/data/user3/sobj/`
- **저장 형식**: `qs::qsave()` 사용
- 모든 테스트 결과는 이 경로에 저장

## 개발 철학 및 목표

### 핵심 목표
1. **기본 분석**: Pseudotime에 따른 gene expression의 증가/감소 패턴 파악
2. **고급 분석**: Metadata와의 gene expression 상관관계 변화 분석
3. **복잡한 패턴**: 단순 증가/감소를 넘어선 복잡한 동적 패턴 탐지

### 주요 방법론
- **Monocle3**: 가장 널리 사용되는 trajectory inference 도구
- **Slingshot**: 유연하고 강력한 trajectory inference
- **최신 방법**: 성능이 검증된 최신 알고리즘 우선 사용

### 기존 코드 기반
- **파일**: `myR/R/pseudotime.R`
- **주요 함수**:
  - `run_slingshot_from_seurat()`: Seurat 객체에서 Slingshot 실행
  - `analyze_gene_dynamics()`: GAM을 이용한 gene expression 동역학 분석
  - `analyze_gene_dynamics_tradeSeq()`: tradeSeq를 이용한 분석
  - `process_gene_list_dynamics()`: 여러 유전자 일괄 처리

## 개발 프로세스

### 작업 순서
1. **다운샘플 데이터로 테스트**
   - 스크립트 실행
   - 오류 디버깅
   - 테스트 코드 정상 작동 확인
   - 데이터 생성 및 저장 확인

2. **원본 데이터로 검증**
   - 동일한 테스트를 원본 데이터에 적용
   - 성능 및 결과 검증

### 커밋 정책
- 작업 완료 후 변경사항 커밋
- 의미 있는 단위로 커밋

## 문서 업데이트

### 우선순위
- **한글 문서 우선**: `*_Korean.md` 파일들을 우선적으로 업데이트
- **영문 문서**: 필요시 동기화

### 주요 문서
- `myR/docs/CHANGELOG_Korean.md`: 변경 이력
- `myR/docs/DEVLOG_Korean.md`: 개발 로그
- `myR/docs/README_Korean.md`: 사용자 가이드

## 참고 사항

### 패키지 로드
- `devtools::load_all()`로 패키지 재로드 필요
- 모든 함수는 `myR` 패키지의 일부

### 데이터 형식
- **입력**: Seurat 객체 또는 SingleCellExperiment 객체
- **출력**: 분석 결과 리스트 또는 데이터프레임
- **시각화**: ggplot2 객체 및 저장된 플롯 파일

### 데이터 전처리
- **전처리 함수**: `preprocess_pseudotime_data()` (파일: `preprocess_data.R`)
- **사용법**:
  ```r
  # 전처리 실행
  sobj_preprocessed <- preprocess_pseudotime_data(
    input_file = "/data/user3/sobj/IS_scvi_251107_ds2500.qs",
    output_file = "/data/user3/sobj/IS_scvi_251107_ds2500_preprocessed.qs"
  )
  ```
- **전처리 내용**:
  - sex 변수: M, F로 변환 (factor)
  - g3 변수: factor로 변환 (1, 2)
  - 시간 변수: datetime으로 파싱 (icu_adm_dt, ia_start, ia_end)
  - 시간 차이 계산: icu_to_ia_start_hours, icu_to_ia_end_hours, ia_duration_hours
- **중요**: 전처리된 파일을 사용하여 분석을 진행하세요!

