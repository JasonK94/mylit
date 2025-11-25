# 실행 가이드

## 빠른 시작

### 1. R 세션 시작
```bash
cd /home/user3/GJC_KDW_250721
R
```

### 2. 함수 로드
```r
# 필수 함수 로드
devtools::load_all("/home/user3/data_user3/git_repo/mylit/myR")
source("/home/user3/data_user3/git_repo/mylit/myR/R/test_analysis.R")

# DEG consensus 함수 로드
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_limma.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_edger.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_methods_deseq2.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_standardize.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/deg_consensus_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus/run_deg_consensus.R")
```

### 3. 테스트 실행
```r
# 단계별 테스트 (권장)
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_step_by_step.R")

# 또는 전체 파이프라인 테스트
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_full_pipeline.R")
```

## 결과 확인

```r
library(qs)

# 결과 로드
result <- qs::qread("/data/user3/sobj/test_deg_consensus_step_by_step.qs")

# 구조 확인
str(result, max.level = 2)

# Consensus DEG 리스트 확인
head(result$consensus_deg_list, 20)
summary(result$consensus_deg_list$agreement)
summary(result$consensus_deg_list$consensus_score)
```

## 문제 해결

### 함수가 로드되지 않음
- 모든 source() 명령이 성공했는지 확인
- 파일 경로가 올바른지 확인

### 방법론 실행 실패
- 개별 방법론을 먼저 테스트
- 오류 메시지 확인
- 필요한 패키지가 설치되어 있는지 확인

### 표준화 실패
- 결과 형식 확인
- 컬럼명 확인

