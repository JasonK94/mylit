# Cell-to-Cell Interaction (CCI) Tool Development Plan

## 개요
scRNAseq 데이터에서 Cell-to-Cell Interaction 분석을 수행하는 통합 도구 개발. NicheNet을 중심으로 하되, 확장 가능한 모듈 구조로 설계.

## 목표
1. **최종 함수**: `sobj`, 분석 대상 클러스터, `deg list(dataframe)` 입력 → CCI 분석 결과 생성
2. **모듈화**: 각 단계별 독립적인 함수로 구성
3. **자동 저장**: 각 단계마다 `/data/user3/sobj`에 중간 결과 저장
4. **문서화**: 표준화된 문서 구조 준수

## 프로젝트 구조

### 디렉터리 구조
```
_wt/cci/
  docs/
    cci/
      cci.md                    # 메인 기능 문서
      TEST_INSTRUCTIONS.md      # 테스트 가이드
      devlog.md                 # 개발 로그
  scripts/
    cci/
      test_cci.R                # 테스트 스크립트
  myR/
    R/
      cci/
        run_cci_analysis.R      # 메인 CCI 분석 함수
        prepare_cci_data.R      # 데이터 준비 및 검증
        nichenet_wrapper.R      # NicheNet 분석 래퍼
        save_cci_results.R       # 결과 저장 유틸리티
        utils_cci.R             # CCI 유틸리티 함수
```

## 개발 단계

### Phase 1: 환경 설정 및 기본 구조
1. **워크트리 및 브랜치 생성**
   - `git worktree add ../_wt/cci cci` (cci 브랜치 생성)
   - 표준 디렉터리 구조 생성

2. **문서 작성**
   - `cci.md`: 기능 상세 설명
   - `TEST_INSTRUCTIONS.md`: 테스트 가이드
   - `devlog.md`: 개발 로그

### Phase 2: 데이터 준비 모듈
**파일**: `myR/R/cci/prepare_cci_data.R`

**함수들**:
- `validate_cci_inputs()`: 입력 검증 (sobj, cluster_col, deg_df)
- `prepare_receiver_deg()`: DEG 리스트에서 receiver DEG 추출
- `identify_sender_receiver()`: 클러스터에서 sender/receiver 식별
- `filter_expressed_genes()`: 발현 유전자 필터링

**출력 저장**:
- `/data/user3/sobj/cci_prepared_data_{timestamp}.qs`

### Phase 3: NicheNet 래퍼 모듈
**파일**: `myR/R/cci/nichenet_wrapper.R`

**함수들**:
- `run_nichenet_for_cci()`: CCI 분석용 NicheNet 실행
  - 기존 `run_nichenet_analysis()` 활용
  - DEG 리스트를 직접 입력으로 받도록 수정
  - sender/receiver 자동 식별

**출력 저장**:
- `/data/user3/sobj/cci_nichenet_results_{timestamp}.qs`

### Phase 4: 메인 CCI 분석 함수
**파일**: `myR/R/cci/run_cci_analysis.R`

**함수**: `run_cci_analysis()`

**입력**:
- `sobj`: Seurat 객체
- `cluster_col`: 클러스터 컬럼명 (예: "anno3.scvi")
- `deg_df`: DEG 결과 데이터프레임 (gene, cluster, logFC, p_val_adj 등)
- `sender_clusters`: sender 클러스터 벡터 (선택적)
- `receiver_cluster`: receiver 클러스터 (선택적)
- `condition_col`: 조건 컬럼 (예: "g3")
- `condition_oi`: 관심 조건 (예: "2")
- `condition_ref`: 참조 조건 (예: "1")
- `species`: "human" 또는 "mouse" (기본: "human")
- `output_dir`: 출력 디렉터리 (기본: NULL)

**처리 단계**:
1. 입력 검증
2. DEG 리스트에서 receiver 클러스터의 DEG 추출
3. Sender 클러스터 자동 식별 (또는 사용자 지정)
4. NicheNet 데이터 로드
5. NicheNet 분석 실행
6. 결과 정리 및 저장

**출력**:
- 리스트 객체:
  - `nichenet_results`: NicheNet 분석 결과
  - `sender_receiver_map`: sender-receiver 매핑
  - `deg_summary`: 사용된 DEG 요약
  - `output_path`: 저장 경로

**출력 저장**:
- `/data/user3/sobj/cci_analysis_results_{timestamp}.qs`

### Phase 5: 결과 저장 유틸리티
**파일**: `myR/R/cci/save_cci_results.R`

**함수들**:
- `save_cci_intermediate()`: 중간 결과 저장
- `save_cci_final()`: 최종 결과 저장
- `load_cci_results()`: 저장된 결과 로드

### Phase 6: 유틸리티 함수
**파일**: `myR/R/cci/utils_cci.R`

**함수들**:
- `extract_receiver_degs()`: DEG 리스트에서 receiver DEG 추출
- `identify_top_senders()`: DEG 기반으로 top sender 식별
- `format_cci_summary()`: 결과 요약 테이블 생성

## 데이터 스펙

### 입력 데이터
- **테스트 데이터**: `/data/user3/sobj/IS6_sex_added_251110_ds2500.qs`
- **원본 데이터**: `/data/user3/sobj/IS6_sex_added_251110.qs`

### 메타데이터 구조 (예상)
- `hos_no`: 샘플 ID
- `g3`: 타겟 그룹 변수 (1, 2, "NA")
- `GEM`: 배치 변수
- `sex`: 생물학적 공변량
- `anno3.scvi`: 클러스터 주석

### DEG 데이터프레임 구조 (예상)
```r
# 필수 컬럼:
# - gene: 유전자명
# - cluster: 클러스터 ID
# - logFC 또는 avg_log2FC: 로그 폴드 체인지
# - p_val_adj: 조정된 p-value
# 선택적 컬럼:
# - p_val: 원시 p-value
# - FDR: FDR 값
```

## 테스트 계획

### Test 1: 기본 기능 테스트
- 다운샘플 데이터 사용
- 단일 sender-receiver 쌍
- 기본 파라미터

### Test 2: 다중 sender 테스트
- 여러 sender 클러스터
- 단일 receiver

### Test 3: 전체 데이터 테스트
- 원본 데이터 사용
- 모든 클러스터 조합

### Test 4: DEG 리스트 직접 입력 테스트
- 외부 DEG 결과 사용
- 다양한 DEG 형식 지원

## 구현 우선순위

1. **High Priority**
   - 워크트리 및 브랜치 생성
   - 기본 문서 작성
   - `run_cci_analysis()` 메인 함수
   - 기본 테스트

2. **Medium Priority**
   - 모듈화 (데이터 준비, NicheNet 래퍼 분리)
   - 자동 저장 기능
   - 결과 요약 함수

3. **Low Priority**
   - 추가 시각화
   - 다른 CCI 도구 통합 (CellChat 등)
   - 배치 처리 기능

## 의존성

### 필수 패키지
- `Seurat`
- `nichenetr`
- `dplyr`
- `qs` (데이터 저장)

### 선택적 패키지
- `circlize` (Circos plot)
- `DiagrammeR` (Signaling path)
- `ComplexHeatmap` (고급 시각화)

## 예상 타임라인

- **Phase 1-2**: 1일 (환경 설정, 데이터 준비)
- **Phase 3-4**: 2일 (NicheNet 래퍼, 메인 함수)
- **Phase 5-6**: 1일 (유틸리티, 저장 기능)
- **테스트 및 디버깅**: 2-3일

**총 예상 기간**: 6-7일

## 다음 단계

1. 이 계획 검수 및 승인
2. 워크트리 생성
3. 기본 구조 및 문서 작성
4. 단계별 구현 및 테스트
5. 디버깅 및 최적화

