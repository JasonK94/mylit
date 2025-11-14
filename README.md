# Multi-Model DEG Consensus Engine

여러 DEG 방법론을 통합하여 robust consensus DEG signature를 생성하는 시스템입니다.

## 개요

이 프로젝트는 g3==2 vs g3==1 그룹 간 DEG 분석을 여러 방법론으로 수행하고, 방법론 수준에서의 클러스터링을 통해 consensus DEG signature를 생성합니다.

## 주요 기능

1. **여러 DEG 방법론 통합 실행**
   - muscat (edgeR, DESeq2, limma-voom, limma-trend)
   - limma 계열 (voom, trend)
   - edgeR 독립 (LRT, QLF)
   - DESeq2 독립 (Wald, LRT)
   - nebula

2. **결과 표준화 및 통합**
   - 다양한 결과 형식을 자동으로 표준화
   - 유전자 × 방법론 행렬 구성

3. **Consensus 분석**
   - 방법론 간 일치도 계산
   - 방법론 수준 클러스터링
   - Consensus DEG signature 생성

## 빠른 시작

```r
# 1. 함수 로드
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 2. 데이터 로드
is5s <- qs::qread("/data/user3/sobj/IS6_sex_added_251110_ds2500.qs")

# 3. 전체 파이프라인 실행
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_full_pipeline.R")
```

## 파일 구조

```
_wt/deg-consensus/
├── myR/R/deg_consensus/
│   ├── run_deg_consensus.R          # 메인 통합 함수
│   ├── deg_methods_limma.R          # limma 계열
│   ├── deg_methods_edger.R          # edgeR 계열
│   ├── deg_methods_deseq2.R         # DESeq2 계열
│   ├── deg_standardize.R            # 결과 표준화
│   └── deg_consensus_analysis.R     # Consensus 분석
├── scripts/
│   ├── test_phase2_limma.R          # limma 테스트
│   └── test_full_pipeline.R         # 전체 파이프라인 테스트
└── docs/deg_consensus/
    ├── DEVELOPMENT_PLAN.md           # 개발 계획
    ├── FINAL_SUMMARY.md             # 최종 요약
    └── COMPLETION_STATUS.md         # 완료 상태
```

## 지원하는 방법론

- muscat-edgeR, muscat-DESeq2, muscat-limma-voom, muscat-limma-trend
- limma-voom, limma-trend
- edgeR-LRT, edgeR-QLF
- DESeq2-Wald, DESeq2-LRT
- nebula

## 문서

자세한 내용은 `docs/deg_consensus/` 폴더를 참조하세요.

