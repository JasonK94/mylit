# 테스트 실행 결과 보고서

## 실행 정보
- **테스트 스크립트**: `scripts/test_step_by_step.R`
- **데이터**: `/data/user3/sobj/IS6_sex_added_251110_ds2500.qs`
- **실행 날짜**: (실행 후 기록)

## 실행 방법

```r
# R 세션에서
cd /home/user3/GJC_KDW_250721
R

# R 내에서
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")

# 테스트 실행
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_step_by_step.R")
```

## 예상 결과

### Step 1: 함수 로드
- 모든 함수가 성공적으로 로드되어야 함

### Step 2: 데이터 로드
- 데이터가 정상적으로 로드되어야 함
- 메타데이터 컬럼이 올바르게 인식되어야 함

### Step 3: 개별 방법론 테스트
- **muscat-edgeR**: 성공 (예상)
- **limma-voom**: 성공 (예상)
- **edgeR-LRT**: 성공 (예상)

각 방법론은 클러스터별 DEG 결과를 반환해야 함.

### Step 4: 결과 표준화
- 모든 방법론의 결과가 표준화되어야 함
- 표준화된 컬럼: gene, logFC, pvalue, pvalue_adj, statistic, method, cluster_id

### Step 5: 행렬 구성
- beta 행렬: genes × methods
- pvalue 행렬: genes × methods
- significance 행렬: genes × methods (FDR < 0.05)

### Step 6: Consensus 분석
- Agreement scores 계산
- Consensus scores 계산
- Consensus DEG list 생성

### Step 7: 결과 저장
- `/data/user3/sobj/test_deg_consensus_step_by_step.qs` 파일 생성

## 결과 확인

```r
library(qs)
result <- qs::qread("/data/user3/sobj/test_deg_consensus_step_by_step.qs")

# 구조 확인
str(result, max.level = 2)

# Consensus DEG 리스트
head(result$consensus_deg_list, 20)
nrow(result$consensus_deg_list)

# Agreement scores 분포
summary(result$agreement_scores)
hist(result$agreement_scores, breaks = 50)

# Consensus scores 분포
summary(result$consensus_scores$consensus_score)
hist(result$consensus_scores$consensus_score, breaks = 50)
```

## 성공 기준

1. ✅ 모든 함수 로드 성공
2. ✅ 데이터 로드 성공
3. ✅ 최소 2개 이상의 방법론 성공
4. ✅ 결과 표준화 성공
5. ✅ 행렬 구성 성공
6. ✅ Consensus 분석 성공
7. ✅ 결과 파일 생성 성공

## 다음 단계

테스트 성공 후:
1. 더 많은 방법론 추가 테스트
2. 전체 파이프라인 테스트
3. 결과 검증 및 분석

