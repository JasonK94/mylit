# DEG Consensus 모듈 개발 계획 (2025)

## 문서 개요

이 문서는 DEG Consensus 모듈의 향후 개발 계획을 정리한 것입니다. 현재 상태를 파악하고, dream과 MAST 방법론 추가, analysis 모듈 통합 등을 포함한 단계별 개발 로드맵을 제시합니다.

**대상 독자**: 다른 AI 에이전트 및 개발자  
**작성일**: 2025-01-XX  
**최종 업데이트**: 2025-01-XX

---

## 1. 현재 상태 (Current Status)

### 1.1 완료된 기능

#### 지원하는 DEG 방법론 (10개)
1. **muscat-edgeR**: `runMUSCAT2_v1(method="edgeR")`
2. **muscat-DESeq2**: `runMUSCAT2_v1(method="DESeq2")`
3. **muscat-limma-voom**: `runMUSCAT2_v1(method="limma-voom")`
4. **muscat-limma-trend**: `runMUSCAT2_v1(method="limma-trend")`
5. **limma-voom**: `runLIMMA_voom_v1()` (독립 구현)
6. **limma-trend**: `runLIMMA_trend_v1()` (독립 구현)
7. **edgeR-LRT**: `runEDGER_LRT_v1()` (독립 구현)
8. **edgeR-QLF**: `runEDGER_QLF_v1()` (독립 구현)
9. **DESeq2-Wald**: `runDESEQ2_Wald_v1()` (독립 구현)
10. **DESeq2-LRT**: `runDESEQ2_LRT_v1()` (독립 구현)

#### 핵심 기능
- ✅ `run_deg_consensus()`: 통합 실행 엔진
- ✅ `standardize_deg_results()`: 결과 표준화
- ✅ `compute_consensus_scores()`: Consensus 점수 계산
- ✅ `aggregate_cluster_deg_consensus()`: 클러스터별 결과 통합
- ✅ `covar_effects` 지원: Fixed effects (sex 등) 보정
- ✅ 작은 클러스터 처리: muscat 실패 시 fallback
- ✅ 중간 데이터 저장: 재분석을 위한 데이터 보존
- ✅ NEBULA 별도 실행: 전체 데이터에서 실행 (공선성 문제 회피)

### 1.2 부분 구현 / 미완성 기능

#### dream 방법론
- **상태**: 코드는 존재하나 (`deg_methods_dream.R`) 작동하지 않음
- **위치**: `myR/R/deg_consensus/deg_methods_dream.R`
- **문제점**: 
  - 원래 잘 작동하지 않았음 (사용자 피드백)
  - 현재 파이프라인에서 `include_dream <- FALSE`로 제외됨
  - 구현이 완전하지 않거나 데이터 구조와 맞지 않을 수 있음

#### MAST 방법론
- **상태**: 코드가 주석 처리됨 (`deg_methods_base.R`)
- **위치**: `myR/R/deg_consensus/deg_methods_base.R` (line 9-132)
- **문제점**:
  - MAST 패키지 호환성 문제 (MAST 1.28.0에서 `FromSeurat` 제거)
  - `runMAST_v1` 함수가 작동하지 않음
  - `analysis` 모듈에도 유사한 코드가 있으나 동일한 문제 존재

### 1.3 analysis 모듈과의 관계

#### 현재 구조
- **analysis 모듈**: `docs/analysis/` 폴더에 문서화
  - NEBULA, muscat을 이용한 Formula 1 분석
  - 단일세포 수준의 Mixed-Effects Model
  - 복잡한 실험 설계 지원 (교호작용, 공변량 보정)
  
- **deg-consensus 모듈**: `docs/deg-consensus-dev/` 폴더에 문서화
  - 여러 DEG 방법론 통합
  - Consensus Signature 도출
  - 클러스터별 결과 통합

#### 중복 및 통합 필요성
- **중복 함수**: `runNEBULA2_v1`, `runMUSCAT2_v1` 등이 양쪽 모듈에 존재
- **통합 이점**:
  - 코드 중복 제거
  - 일관된 API 제공
  - 유지보수 용이성 향상
  - 사용자 혼란 감소

---

## 2. 개발 목표 (Development Goals)

### 2.1 단기 목표 (Phase 1: 1-2주)

#### 2.1.1 dream 방법론 완성
- **목표**: dream 방법론을 정상 작동하도록 수정 및 테스트
- **작업 내용**:
  1. `deg_methods_dream.R` 코드 검토 및 디버깅
  2. `run_deg_consensus()`에 dream 통합
  3. 작은 데이터셋으로 테스트
  4. 오류 수정 및 문서화

- **예상 문제점**:
  - Random effects 처리 방식 확인 필요
  - VariancePartition 패키지 호환성
  - Pseudobulk 생성 방식 확인

#### 2.1.2 MAST 방법론 재구현
- **목표**: MAST 패키지 최신 버전과 호환되도록 재구현
- **작업 내용**:
  1. MAST 패키지 최신 문서 확인
  2. `FromMatrix` 직접 사용 또는 다른 변환 방법 구현
  3. `runMAST_v1` 함수 수정
  4. `run_deg_consensus()`에 MAST 통합
  5. 테스트 및 문서화

- **참고 자료**:
  - `docs/archive/general/NEXT_AGENT_PROMPT.md`: MAST 호환성 문제 언급
  - `myR/R/deg_consensus/deg_methods_base.R`: 주석 처리된 코드

### 2.2 중기 목표 (Phase 2: 2-4주)

#### 2.2.1 analysis 모듈 통합
- **목표**: analysis 모듈의 핵심 기능을 deg-consensus에 통합
- **작업 내용**:
  1. **코드 통합**:
     - `myR/R/analysis.R`의 함수들을 `myR/R/deg_consensus/`로 이동
     - 중복 함수 제거 및 통합
     - API 일관성 확보
     
  2. **문서 통합**:
     - `docs/analysis/` 내용을 `docs/deg-consensus-dev/`로 통합
     - Formula 1 분석 방법을 deg-consensus 가이드에 추가
     - 워크플로우 다이어그램 업데이트
     
  3. **테스트 및 검증**:
     - 기존 analysis 모듈 기능이 정상 작동하는지 확인
     - 통합 후 회귀 테스트 수행

- **통합 전략**:
  ```
  현재:
  - analysis 모듈: Formula 1 분석 (NEBULA 중심)
  - deg-consensus 모듈: Multi-method consensus
  
  통합 후:
  - deg-consensus 모듈: 
    * Multi-method consensus (기존 기능)
    * Formula-based analysis (analysis 모듈 기능 통합)
    * 통합된 API로 두 가지 접근 방식 모두 지원
  ```

#### 2.2.2 문서 통합 및 정리
- **목표**: 통합된 모듈에 대한 명확한 문서 제공
- **작업 내용**:
  1. `DEG_CONSENSUS_INTEGRATED_GUIDE_KR.md` 업데이트
  2. Formula 1 분석 방법 추가
  3. 방법론별 사용 가이드 정리
  4. 예제 코드 업데이트

### 2.3 장기 목표 (Phase 3: 1-2개월)

#### 2.3.1 추가 방법론 지원
- **목표**: 더 다양한 DEG 방법론 지원
- **후보 방법론**:
  - **scran**: Single-cell 분석 방법
  - **scde**: Single-cell differential expression
  - **BPSC**: Beta-Poisson model for single-cell
  - **SCOPE**: Single-cell analysis method

#### 2.3.2 성능 최적화
- **목표**: 대규모 데이터셋에서의 성능 향상
- **작업 내용**:
  1. 병렬 처리 최적화
  2. 메모리 사용량 최적화
  3. 캐싱 메커니즘 도입

---

## 3. 상세 개발 계획 (Detailed Development Plan)

### 3.1 Phase 1: dream 및 MAST 방법론 완성

#### 3.1.1 dream 방법론 개발

**현재 코드 위치**: `myR/R/deg_consensus/deg_methods_dream.R`

**개발 단계**:

1. **코드 검토** (1일)
   - 현재 구현 상태 확인
   - 오류 원인 파악
   - 필요한 수정 사항 정리

2. **디버깅 및 수정** (2-3일)
   - Random effects 처리 확인
   - VariancePartition 패키지 사용법 확인
   - Pseudobulk 생성 방식 확인
   - 오류 수정

3. **통합** (1일)
   - `run_deg_consensus()`에 dream 추가
   - Method handler 등록
   - 표준화 함수 추가

4. **테스트** (1-2일)
   - 작은 데이터셋으로 테스트
   - 오류 수정
   - 문서화

**예상 결과물**:
- `runDREAM_v1()` 함수 완성
- `run_deg_consensus()`에서 dream 사용 가능
- 테스트 리포트

#### 3.1.2 MAST 방법론 재구현

**현재 코드 위치**: `myR/R/deg_consensus/deg_methods_base.R` (주석 처리됨)

**개발 단계**:

1. **MAST 패키지 조사** (1일)
   - 최신 버전 문서 확인
   - `FromMatrix` 사용법 확인
   - 대안 방법 조사

2. **재구현** (2-3일)
   - `runMAST_v1()` 함수 재작성
   - Seurat → SingleCellAssay 변환 수정
   - 테스트

3. **통합** (1일)
   - `run_deg_consensus()`에 MAST 추가
   - Method handler 등록
   - 표준화 함수 추가

4. **테스트** (1-2일)
   - 작은 데이터셋으로 테스트
   - 오류 수정
   - 문서화

**예상 결과물**:
- `runMAST_v1()` 함수 완성
- `run_deg_consensus()`에서 MAST 사용 가능
- 테스트 리포트

### 3.2 Phase 2: analysis 모듈 통합

#### 3.2.1 코드 통합

**작업 순서**:

1. **함수 매핑 및 중복 확인** (1일)
   ```
   analysis 모듈 → deg-consensus 모듈
   - runNEBULA2_v1: 이미 존재 (통합 필요)
   - runMUSCAT2_v1: 이미 존재 (통합 필요)
   - runMAST_v1: 주석 처리됨 (재구현 필요)
   - Formula 1 관련 함수: 새로 추가 필요
   ```

2. **함수 이동 및 통합** (2-3일)
   - 중복 함수 제거
   - 공통 함수 추출
   - API 일관성 확보

3. **의존성 정리** (1일)
   - Import 문 정리
   - 패키지 의존성 확인

#### 3.2.2 문서 통합

**작업 순서**:

1. **현재 문서 분석** (1일)
   - `docs/analysis/ANALYSIS_INTEGRATED_GUIDE_KR.md` 내용 파악
   - `docs/deg-consensus-dev/DEG_CONSENSUS_INTEGRATED_GUIDE_KR.md`와 비교
   - 통합 전략 수립

2. **문서 통합** (2-3일)
   - Formula 1 분석 섹션 추가
   - 워크플로우 다이어그램 업데이트
   - 예제 코드 통합

3. **문서 정리** (1일)
   - 중복 내용 제거
   - 구조 개선
   - 링크 업데이트

#### 3.2.3 테스트 및 검증

**작업 순서**:

1. **회귀 테스트** (2-3일)
   - 기존 analysis 모듈 기능 테스트
   - 기존 deg-consensus 기능 테스트
   - 통합 후 기능 테스트

2. **통합 테스트** (2-3일)
   - 전체 파이프라인 테스트
   - 다양한 데이터셋으로 테스트
   - 성능 측정

---

## 4. 기술적 고려사항 (Technical Considerations)

### 4.1 dream 방법론

#### 4.1.1 Random Effects 처리
- **문제**: dream은 Random effects를 지원하지만, 현재 파이프라인은 클러스터별로 쪼갠 데이터를 사용
- **해결책**:
  - 클러스터별 실행 시: 클러스터 내에서 random effects 처리
  - 전체 데이터 실행 시: 클러스터를 random effect로 처리 가능

#### 4.1.2 VariancePartition 패키지
- **의존성**: `variancePartition` 패키지 필요
- **사용법**: `dream()` 함수를 통해 Linear Mixed Model 실행
- **주의사항**: Pseudobulk 데이터에 적용

### 4.2 MAST 방법론

#### 4.2.1 패키지 호환성
- **문제**: MAST 1.28.0에서 `FromSeurat` 제거
- **해결책**:
  ```r
  # 옵션 1: FromMatrix 직접 사용
  sca <- MAST::FromMatrix(
    exprsArray = counts_matrix,
    cData = cell_metadata,
    fData = gene_metadata,
    check_sanity = FALSE
  )
  
  # 옵션 2: SingleCellExperiment → SingleCellAssay 변환
  sce <- as.SingleCellExperiment(sobj)
  sca <- MAST::SceToSingleCellAssay(sce)
  ```

#### 4.2.2 Hurdle Model
- **특징**: MAST는 hurdle model을 사용 (발현 여부 + 발현량)
- **주의사항**: Single-cell 수준에서 작동 (pseudobulk 아님)

### 4.3 analysis 모듈 통합

#### 4.3.1 API 일관성
- **현재 문제**: 
  - analysis 모듈: Formula 기반 API
  - deg-consensus 모듈: 변수명 기반 API
  
- **해결책**:
  ```r
  # 통합된 API 예시
  run_deg_consensus(
    sobj = sobj,
    methods = c("limma-voom", "edgeR-LRT", "dream", "mast"),
    # 옵션 1: 변수명 기반 (기존 방식)
    group_id = "g3",
    covar_effects = "sex",
    # 옵션 2: Formula 기반 (analysis 모듈 방식)
    formula = ~ g3 + sex + (1|hos_no),
    contrast = "2 - 1"
  )
  ```

#### 4.3.2 함수 중복 해결
- **전략**: 
  1. 공통 함수는 `deg_methods_base.R`에 유지
  2. 모듈별 특화 함수는 각 모듈에 유지
  3. Wrapper 함수로 통합 API 제공

---

## 5. 테스트 계획 (Testing Plan)

### 5.1 단위 테스트

#### 5.1.1 dream 방법론
- **테스트 데이터**: 작은 클러스터 (예: 100-500 cells)
- **테스트 시나리오**:
  1. 기본 실행 테스트
  2. Random effects 포함 테스트
  3. Covariates 포함 테스트
  4. 오류 처리 테스트

#### 5.1.2 MAST 방법론
- **테스트 데이터**: 작은 클러스터 (예: 100-500 cells)
- **테스트 시나리오**:
  1. 기본 실행 테스트
  2. Formula 기반 테스트
  3. 오류 처리 테스트

### 5.2 통합 테스트

#### 5.2.1 전체 파이프라인 테스트
- **테스트 데이터**: 
  - Downsampled data: `/data/user3/sobj/IS6_sex_added_0.1x_251110.qs`
  - Full data: `/data/user3/sobj/IS6_sex_added_251110.qs`
  
- **테스트 시나리오**:
  1. 모든 방법론 포함 실행
  2. 클러스터별 실행
  3. Consensus 계산
  4. 결과 저장 및 로드

#### 5.2.2 analysis 모듈 통합 테스트
- **테스트 시나리오**:
  1. Formula 1 분석 실행
  2. Multi-method consensus 실행
  3. 두 방식 결과 비교

---

## 6. 문서화 계획 (Documentation Plan)

### 6.1 함수 문서화
- 모든 새로 추가/수정된 함수에 roxygen2 주석 추가
- 예제 코드 포함
- 매개변수 설명 상세화

### 6.2 사용자 가이드 업데이트
- `DEG_CONSENSUS_INTEGRATED_GUIDE_KR.md` 업데이트
- dream, MAST 방법론 사용법 추가
- Formula 기반 분석 방법 추가
- 예제 코드 업데이트

### 6.3 개발자 가이드
- 코드 구조 설명
- 새로운 방법론 추가 방법
- 테스트 방법

---

## 7. 일정 (Timeline)

### Phase 1: dream 및 MAST 완성 (2주)
- **Week 1**: dream 방법론 개발 및 테스트
- **Week 2**: MAST 방법론 재구현 및 테스트

### Phase 2: analysis 모듈 통합 (3-4주)
- **Week 1**: 코드 통합
- **Week 2**: 문서 통합
- **Week 3-4**: 테스트 및 검증

### Phase 3: 추가 개선 (선택적, 1-2개월)
- 추가 방법론 지원
- 성능 최적화

---

## 8. 리스크 및 대응 방안 (Risks and Mitigation)

### 8.1 dream 방법론
- **리스크**: 구현이 복잡하거나 패키지 호환성 문제
- **대응**: 단계별 개발, 작은 데이터셋으로 먼저 테스트

### 8.2 MAST 방법론
- **리스크**: 패키지 API 변경으로 인한 호환성 문제
- **대응**: 최신 문서 확인, 대안 방법 조사

### 8.3 analysis 모듈 통합
- **리스크**: 기존 기능 손상
- **대응**: 철저한 회귀 테스트, 단계별 통합

---

## 9. 참고 자료 (References)

### 9.1 프로젝트 문서
- `docs/deg-consensus-dev/DEG_CONSENSUS_INTEGRATED_GUIDE_KR.md`: 현재 가이드
- `docs/analysis/ANALYSIS_INTEGRATED_GUIDE_KR.md`: analysis 모듈 가이드
- `docs-main/guide_KR.md`: 프로젝트 전체 가이드
- `docs/deg-consensus-dev/archive/DEVELOPMENT_PLAN.md`: 이전 개발 계획

### 9.2 코드 위치
- `myR/R/deg_consensus/`: deg-consensus 모듈 코드
- `myR/R/analysis.R`: analysis 모듈 코드 (통합 예정)
- `scripts/deg-consensus-dev/`: 테스트 스크립트

### 9.3 외부 자료
- dream: `variancePartition` 패키지 문서
- MAST: MAST 패키지 최신 문서
- NEBULA: NEBULA 패키지 문서

---

## 10. 다음 단계 (Next Steps)

### 즉시 시작 가능한 작업
1. **dream 방법론 코드 검토**
   - `myR/R/deg_consensus/deg_methods_dream.R` 파일 읽기
   - 오류 원인 파악
   - 수정 계획 수립

2. **MAST 패키지 조사**
   - 최신 버전 문서 확인
   - `FromMatrix` 사용법 확인
   - 테스트 코드 작성

3. **analysis 모듈 코드 분석**
   - `myR/R/analysis.R` 파일 분석
   - 중복 함수 확인
   - 통합 전략 수립

---

## 부록: 대화 컨텍스트 요약

### 최근 완료된 작업 (2025-01-XX)
1. **`covar_effects` 지원**: Fixed effects (sex 등) 보정 기능 추가
2. **작은 클러스터 처리**: muscat 실패 시 fallback 추가
3. **건너뛴 클러스터 추적**: 이유와 함께 기록
4. **중간 데이터 저장**: 재분석을 위한 데이터 보존
5. **NEBULA 별도 실행**: 전체 데이터에서 실행하여 공선성 문제 회피
6. **문서화**: 방법론별 이슈와 해결책 정리

### 사용자 피드백
- "dream과 mast도 개발하면 좋겠는데"
- "analysis는 흡수 통합하는 게 낫지 않을까"
- "개발 계획을 나에게 제안해봐. docs 폴더에서 나와의 대화 컨텍스트, 그리고 이 프로젝트 자체의 컨텍스트(.md파일 등으로 정리된)를 포함하여 .md 파일로 만들어. 다른 에이전트와 개발 진행하게."

### 현재 상태 요약
- **지원 방법론**: 10개 (muscat 4개, limma 2개, edgeR 2개, DESeq2 2개)
- **부분 구현**: dream (코드 존재하나 작동 안 함), MAST (주석 처리됨)
- **통합 필요**: analysis 모듈과의 통합

---

**문서 작성자**: AI Assistant  
**검토 필요**: 사용자 검토 후 수정  
**업데이트 주기**: 개발 진행에 따라 업데이트

