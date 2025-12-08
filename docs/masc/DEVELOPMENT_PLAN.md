# MASC 파이프라인 개발 계획서

**작성일**: 2025-12-08  
**버전**: 1.0  
**목적**: MASC (Mixed-effects Association testing for Single Cells) 파이프라인 구현 및 통합

---

## 📚 개요

### MASC란?
MASC (Mixed-effects Association testing for Single Cells)는 단일세포 데이터에서 클러스터 abundance와 공변량(예: 질병 상태, 치료 그룹) 간의 연관성을 테스트하는 혼합효과 모델링 방법입니다.

**주요 특징**:
- 로지스틱 혼합효과 모델 (logistic mixed-effects model) 사용
- 클러스터별로 null 모델과 full 모델 비교 (LRT)
- 랜덤 효과(예: 개체, 배치) 및 고정 효과(예: 성별, 나이) 고려
- 각 클러스터별 Odds Ratio (OR) 및 신뢰구간 계산

### 참고 자료
- **원본 패키지**: https://github.com/immunogenomics/masc
- **원본 논문**: nihms-1009396.pdf
- **사용 예시**: [TAURUS paper - LR analysis](https://github.com/DendrouLab/TAURUS_paper/blob/main/3_cell_cell_interactions/LR_analysis_CD_multiNicheNet_longitudinal_subbucket.r)

---

## 🎯 개발 목표

### Phase 1: 핵심 기능 구현
1. **MASC 핵심 함수**
   - 원본 MASC 함수를 기반으로 Seurat 객체 지원 추가
   - 메타데이터 자동 추출 및 검증
   - 에러 핸들링 강화

2. **파이프라인 함수**
   - `run_masc_pipeline()`: 전체 워크플로우 통합
   - Seurat 객체 입력 지원
   - 결과 저장 및 캐싱 기능
   - 시각화 옵션

3. **헬퍼 함수**
   - 데이터 준비 및 검증
   - 결과 포맷팅 및 시각화
   - 모델 저장/로드 유틸리티

### Phase 2: 고급 기능
1. **다중 비교 보정** (FDR, Bonferroni)
2. **시각화 함수** (클러스터별 OR 플롯, p-value 히트맵)
3. **배치 처리** (여러 contrast 변수 동시 분석)

---

## 📋 함수 명세

### 1. 메인 파이프라인 함수

#### `run_masc_pipeline()`
Seurat 객체 또는 데이터프레임에서 MASC 분석을 수행하는 통합 파이프라인

**입력**:
- `seurat_obj`: Seurat 객체 (또는 `seurat_qs_path`로 파일 경로)
- `cluster_var`: 클러스터 할당 컬럼명
- `contrast_var`: 테스트할 주요 공변량 (예: "status", "treatment")
- `random_effects`: 랜덤 효과 변수 벡터
- `fixed_effects`: 고정 효과 변수 벡터
- `save`: 중간 결과 저장 여부
- `force_run`: 캐시 무시하고 재실행 여부
- `plotting`: 시각화 생성 여부

**출력**:
- `masc_results`: MASC 결과 데이터프레임
- `models`: (선택적) 각 클러스터별 모델 객체
- `plots`: (선택적) 시각화 객체 리스트

---

### 2. 핵심 MASC 함수

#### `run_masc_analysis()`
원본 MASC 함수를 래핑하여 Seurat 객체 지원 추가

**기능**:
- 데이터프레임 또는 Seurat 객체 입력 처리
- 클러스터별 로지스틱 GLMM 피팅
- Likelihood Ratio Test (LRT) 수행
- OR 및 95% CI 계산

---

### 3. 헬퍼 함수

- `.masc_prepare_data()`: 데이터 준비 및 검증
- `.masc_format_results()`: 결과 포맷팅
- `.masc_plot_results()`: 시각화 생성
- `.masc_save_models()`: 모델 저장/로드

---

## 🔧 구현 세부사항

### 의존성 패키지
- `lme4`: 혼합효과 모델링
- `Seurat`: 단일세포 데이터 처리
- `qs`: 빠른 직렬화
- `cli`: 진행 상황 메시지
- `ggplot2`: 시각화

### 에러 핸들링
1. 입력 검증: 필수 변수 존재 확인, factor 타입 확인
2. 모델 수렴 문제: 실패 시 경고 및 NA 반환
3. 데이터 부족: 최소 셀 수 확인

---

## ✅ 구현 체크리스트

### Phase 1: 핵심 기능
- [ ] `run_masc_analysis()`: 기본 MASC 분석 함수
- [ ] `.masc_prepare_data()`: 데이터 준비 헬퍼
- [ ] `run_masc_pipeline()`: 파이프라인 통합
- [ ] 입력 검증 및 에러 핸들링

### Phase 2: 고급 기능
- [ ] FDR 보정 옵션
- [ ] 시각화 함수
- [ ] 모델 저장/로드 유틸리티

---

## 📝 참고사항

1. **원본 MASC 패키지와의 차이점**:
   - Seurat 객체 직접 지원
   - 파이프라인 형태로 통합
   - 캐싱 및 재현성 강화

2. **기존 코드와의 일관성**:
   - `milo_pipeline.R`의 구조 및 컨벤션 준수
   - roxygen2 문서화 형식 통일
