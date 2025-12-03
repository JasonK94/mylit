# Demultiplexing Utils 비교 분석

## 기존 단순 버전 vs 현재 복잡 버전

### 1. `get_best_two` 함수

**기존:**
```r
get_best_two <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  best <- sorted_indices[1]
  second_best <- sorted_indices[2]
  return(c(best, second_best))
}
```

**현재:**
```r
get_best_two <- function(probs) {
  sorted_idx <- order(probs, decreasing = TRUE)
  list(
    best = probs[sorted_idx[1]],
    best_idx = sorted_idx[1],
    second = probs[sorted_idx[2]],
    second_idx = sorted_idx[2]
  )
}
```

**차이점:**
- 기존: 인덱스만 반환
- 현재: 값과 인덱스 모두 반환 (더 유용)

**필요성:** ✅ 필요 - 확률 값도 필요함

### 2. `get_barcode_mapping` 함수

**기존:**
- 단순히 Best_Sample, Best_Probability, Second_Best_Sample, Second_Best_Probability, Probability_Ratio 반환
- threshold 기반 필터링 없음

**현재:**
- singlet_threshold, doublet_threshold 기반 필터링
- return_probs 옵션
- Negative 셀 처리

**필요성:** ✅ 필요 - threshold 기반 필터링이 중요함

### 3. `generate_sample_names` 함수

**기존:**
```r
generate_sample_names=function(vector){
  singlets=vector
  doublets=combn(vector, 2, FUN=function(x) paste(x,collapse="+"))
  return(c(singlets,doublets))
}
```

**현재:**
- 여러 버전 존재 (format 옵션 등)
- 더 복잡함

**필요성:** ⚠️ 부분적 - 샘플 개수 역산 기능 추가 필요

### 4. 추가된 함수들

- `filter_demulti_results`: 필터링 유틸리티
- `summarize_demulti_results`: 요약 통계
- `parse_doublet_name`: doublet 파싱

**필요성:** ⚠️ 선택적 - 유용하지만 필수는 아님

## 결론

대부분의 복잡화는 **필요함**. 다만 `generate_sample_names`는 샘플 개수 역산 기능이 필요함.

