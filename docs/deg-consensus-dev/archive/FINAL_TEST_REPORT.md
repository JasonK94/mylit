# 최종 테스트 실행 보고서

## 테스트 실행 방법

### 방법 1: 단계별 테스트 (권장)
```r
# R 세션에서
cd /home/user3/GJC_KDW_250721
R

# R 내에서
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_step_by_step.R")
```

### 방법 2: 전체 파이프라인 테스트
```r
source("/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/test_full_pipeline.R")
```

### 방법 3: 셸 스크립트 실행
```bash
/home/user3/data_user3/git_repo/_wt/deg-consensus/scripts/deg-consensus/run_full_test.sh
```

## 예상 출력 파일

테스트 성공 시 다음 파일들이 생성됩니다:

1. **단계별 테스트 결과**
   - `/data/user3/sobj/test_deg_consensus_step_by_step.qs`

2. **전체 파이프라인 결과**
   - `/data/user3/sobj/test_deg_consensus_full_pipeline.qs`
   - `/data/user3/sobj/test_deg_consensus_phase3.qs`
   - `/data/user3/sobj/test_deg_consensus_final_result.qs`

## 결과 확인 방법

```r
library(qs)

# 결과 로드
result <- qs::qread("/data/user3/sobj/test_deg_consensus_step_by_step.qs")

# 개별 방법론 결과 확인
str(result$individual_results, max.level = 1)

# 표준화된 결과 확인
str(result$standardized_results, max.level = 1)

# 행렬 확인
dim(result$deg_matrices$beta)
dim(result$deg_matrices$significance)

# Consensus DEG 리스트 확인
head(result$consensus_deg_list)
nrow(result$consensus_deg_list)
```

## 주요 검증 사항

1. **함수 로드**: 모든 함수가 정상적으로 로드되는가?
2. **데이터 로드**: 데이터가 정상적으로 로드되는가?
3. **방법론 실행**: 각 방법론이 정상적으로 실행되는가?
4. **결과 표준화**: 결과가 올바르게 표준화되는가?
5. **행렬 구성**: 행렬이 올바르게 구성되는가?
6. **Consensus 분석**: Consensus 분석이 정상적으로 수행되는가?
7. **결과 저장**: 결과가 올바르게 저장되는가?

## 알려진 제한사항

1. **함수 의존성**: 모든 함수 파일을 source해야 함
2. **패키지 의존성**: Seurat, muscat, limma, edgeR, DESeq2 등 필요
3. **메모리**: 큰 데이터셋의 경우 메모리 사용량 주의

## 디버깅 팁

1. **단계별 실행**: test_step_by_step.R을 사용하여 문제 지점 파악
2. **traceback()**: 오류 발생 시 traceback()으로 스택 확인
3. **개별 테스트**: 각 방법론을 개별적으로 테스트
4. **로그 확인**: 각 단계의 메시지 출력 확인

