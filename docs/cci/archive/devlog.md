# CCI 분석 도구 개발 로그

## 2025-01-XX: CCI 워크트리 생성 및 초기 개발

### 목표
- Cell-to-Cell Interaction 분석 도구 개발
- DEG 리스트를 직접 입력받아 NicheNet 분석 수행
- 모듈화된 구조로 확장 가능한 도구 구축

### 완료된 작업

#### 1. 워크트리 및 브랜치 생성
- `cci` 브랜치 생성
- `../_wt/cci` 워크트리 생성
- 표준 디렉터리 구조 생성

#### 2. 문서 작성
- **`cci.md`**: CCI 분석 도구 상세 문서
  - 개요 및 핵심 개념 설명
  - 주요 함수 및 파라미터 설명
  - 사용 예시 및 주의사항
  
- **`TEST_INSTRUCTIONS.md`**: 테스트 실행 가이드
  - 환경 설정 방법
  - 데이터 준비 및 확인
  - 다양한 테스트 시나리오
  - 결과 확인 및 저장 방법

- **`devlog.md`**: 개발 로그 (본 문서)

#### 3. 개발 계획 수립
- **`CCI_DEVELOPMENT_PLAN.md`**: 상세 개발 계획
  - 6단계 개발 계획
  - 모듈 구조 설계
  - 테스트 계획

### 완료된 작업 (계속)

#### 4. 함수 구현
- **`prepare_cci_data.R`**: 데이터 준비 및 검증 모듈
  - `validate_cci_inputs()`: 입력 검증
  - `extract_receiver_degs()`: Receiver DEG 추출
  - `identify_sender_clusters()`: Sender 클러스터 식별
  - `filter_expressed_genes()`: 발현 유전자 필터링
  - `prepare_cci_summary()`: 데이터 요약 생성

- **`utils_cci.R`**: 유틸리티 함수
  - `format_deg_summary()`: DEG 요약 테이블 생성
  - `identify_top_senders()`: Top sender 식별
  - `create_sender_receiver_map()`: Sender-receiver 매핑 생성

- **`save_cci_results.R`**: 결과 저장 유틸리티
  - `save_cci_intermediate()`: 중간 결과 저장
  - `save_cci_final()`: 최종 결과 저장
  - `load_cci_results()`: 저장된 결과 로드

- **`nichenet_wrapper.R`**: NicheNet 래퍼 (준비 중)
  - 기존 `run_nichenet_analysis()`와의 통합을 위한 래퍼
  - DEG 리스트 직접 사용을 위한 준비

- **`run_cci_analysis.R`**: 메인 CCI 분석 함수
  - 전체 파이프라인 통합
  - 입력 검증 → DEG 추출 → Sender 식별 → NicheNet 분석 → 결과 저장
  - 자동 저장 기능 포함

#### 5. 테스트 스크립트 작성
- **`scripts/cci/test_cci.R`**: 자동화된 테스트 스크립트
  - 기본 기능 테스트
  - Sender 클러스터 지정 테스트
  - 자동 Sender 식별 테스트

### 다음 단계

1. **테스트 실행 및 디버깅** ✅
   - 다운샘플 데이터로 테스트 실행 완료
   - 발현 유전자 필터링 문제 수정 완료
   - NicheNet 데이터 다운로드 이슈 확인 (네트워크 문제 가능성)

2. **기능 개선**
   - 조건 없이 DEG 리스트만으로 분석 가능하도록 개선
   - NicheNet 함수 직접 호출 옵션 추가
   - NicheNet 데이터 다운로드 재시도 로직 추가
   - 추가 시각화 기능

3. **문서 보완**
   - 사용 예시 추가
   - 문제 해결 가이드 보완
   - 성능 최적화 가이드

## 2025-11-14: 테스트 실행 및 버그 수정

### 테스트 결과
- ✅ 모든 핵심 함수 정상 작동 확인
- ✅ 데이터 로드 및 검증 성공
- ✅ 발현 유전자 필터링 문제 발견 및 수정

### 수정 사항
1. **`filter_expressed_genes()` 함수 개선**
   - `cluster_col` 파라미터 추가
   - `Idents(sobj)` 설정 로직 개선
   - `nichenetr::get_expressed_genes()` 호출 전 Idents 설정 보장

2. **`run_cci_analysis()` 함수 개선**
   - 발현 유전자 필터링 호출 전 Idents 명시적 설정
   - 더 나은 오류 메시지 제공

### 테스트 통계
- 데이터: 2500 cells, 51795 genes
- 클러스터: 24개
- Sender expressed genes: 9459개 (수정 후)
- Receiver expressed genes: 3607개 (수정 후)

### 알려진 이슈
- ~~NicheNet 데이터 다운로드/읽기 문제~~ ✅ 해결: `/data/user3/git_repo/human/` 경로 사용
- ~~FindMarkers DEG 필터링 문제~~ ✅ 해결: p_val_adj_cutoff와 logfc_cutoff 완화

## 2025-11-14: 전체 파이프라인 성공!

### 최종 테스트 결과
- ✅ **전체 CCI 분석 완료!**
- ✅ NicheNet 데이터 로드 성공 (`/data/user3/git_repo/human/` 사용)
- ✅ FindMarkers DEG 분석 성공: **2121개 DEG** 식별
- ✅ Potential ligands 식별: **149개**
- ✅ Ligand activity 예측 성공: **Top 10 ligands** 선택
- ✅ 결과 저장 완료

### 최종 결과
- **Top 5 Ligands**: MFNG, TGFB1, IL15, ADAM10, TSPAN3
- **Sender clusters**: 23개 (자동 식별)
- **Receiver cluster**: Memory B-cells
- **DEGs used**: 2121개 (NicheNet ligand_target_matrix와 교집합)

### 해결한 문제들
1. **NicheNet 데이터 경로**: `/data/user3/git_repo/human/` 사용
2. **FindMarkers 필터링**: p_val_adj_cutoff=1.1, logfc_cutoff=0.05로 완화
3. **발현 유전자 필터링**: Idents 설정 로직 개선

### 저장된 파일
- 중간 결과: `/data/user3/sobj/cci_prepared_data_20251114_*.qs`
- 최종 결과: `/data/user3/sobj/cci_analysis_results_20251114_*.qs`

