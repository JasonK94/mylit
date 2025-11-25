# compute_meta_gene_importance 함수 개선사항

## 개요

`compute_meta_gene_importance` 함수는 TML7의 결과물에서 메타 모델의 변수 중요도(variable importance)를 추출하여 각 L1 signature의 중요도를 계산하고, 이를 gene-level importance로 변환하는 함수입니다.

## 주요 개선사항

### 1. `target_model` 파라미터 추가

기존에는 `best_model`만 사용할 수 있었으나, 이제 `trained_models`에서 원하는 모델을 선택할 수 있습니다.

```r
# Best model 사용 (기본)
cmgi <- compute_meta_gene_importance(tmla)

# 특정 모델 지정
cmgi_glm <- compute_meta_gene_importance(tmla, target_model = "glm")
cmgi_glmnet <- compute_meta_gene_importance(tmla, target_model = "glmnet")
cmgi_ranger <- compute_meta_gene_importance(tmla, target_model = "ranger")
```

**사용 사례:**
- `ranger`나 `earth` 같은 모델이 best model이지만 importance가 없는 경우
- 다른 모델의 importance를 비교하고 싶은 경우
- 특정 모델의 성능을 분석하고 싶은 경우

### 2. Graceful Error Handling

이전에는 importance 추출 실패 시 `stop()`으로 에러가 발생했으나, 이제 `warning()` 후 `NULL`을 반환합니다.

**개선된 동작:**
- `ranger` 모델에서 importance가 없는 경우 → `NULL` 반환 (경고만)
- `earth` 모델에서 `caret::varImp()` 실패 시 → `NULL` 반환 (경고만)
- 기타 모델에서 importance 추출 실패 시 → `NULL` 반환 (경고만)

**장점:**
- List 처리 시 일부 모델 실패해도 전체 작업이 중단되지 않음
- `tryCatch()`로 쉽게 처리 가능

**사용 예시:**
```r
# List 처리 시 NULL 반환되는 경우 필터링
cmgi_list <- lapply(tml_list, function(x) {
  tryCatch(
    compute_meta_gene_importance(x),
    error = function(e) {
      warning("Error computing importance: ", conditionMessage(e))
      NULL
    },
    warning = function(w) {
      # Warning이 있어도 NULL 반환하는 경우가 있음
      NULL
    }
  )
})
cmgi_list <- cmgi_list[!sapply(cmgi_list, is.null)]
```

### 3. 향상된 이름 매칭 로직

`l2_train`의 컬럼 이름은 `make.names()`로 변환될 수 있지만, `l1_signatures`는 원본 이름을 유지합니다. 이제 두 가지 경우를 모두 처리합니다:

1. **직접 매칭**: `sig_names`가 `l1_signatures`의 원본 이름과 직접 일치하는 경우
2. **변환 매칭**: `sig_names`가 `make.names()`로 변환된 이름인 경우, 원본 이름으로 매핑

**처리 로직:**
```r
# 1. 직접 매칭 시도
direct_match <- sig_names %in% l1_sig_names
if (all(direct_match)) {
  sig_names_for_l1 <- sig_names  # 원본 이름 그대로 사용
} else {
  # 2. make.names() 변환 매칭
  name_mapping <- setNames(l1_sig_names, make.names(l1_sig_names))
  original_sig_names <- name_mapping[sig_names]
  sig_names_for_l1 <- as.character(original_sig_names[!is.na(original_sig_names)])
}
```

**장점:**
- 기존 TML 객체와 새로 생성된 TML 객체 모두 지원
- `make.names()` 변환 여부와 관계없이 올바른 매칭

### 4. 모델별 Importance 처리 개선

#### GLM 모델
- `stats::coef()`로 회귀 계수 추출
- `(Intercept)` 제외
- 매칭되지 않으면 `NULL` 반환

#### Ranger 모델
1. `ranger::importance()` 직접 시도
2. 실패 시 `caret::varImp()` 시도
3. 둘 다 실패하면 `NULL` 반환 (경고)

**참고:** `ranger` 모델에서 importance를 얻으려면 학습 시 `importance = "permutation"` 또는 `importance = "impurity"`를 설정해야 합니다. TML7에서는 `importance = "permutation"`으로 설정하지만, 기존 객체는 설정이 없을 수 있습니다.

#### Earth 모델
- `caret::varImp()` 시도
- 실패 시 `NULL` 반환 (경고)

**참고:** `earth` 모델은 MARS (Multivariate Adaptive Regression Splines) 모델로, 기본적으로 `caret::varImp()`를 지원하지 않습니다.

#### 기타 모델 (glmnet, xgbTree, nnet 등)
- `caret::varImp()` 사용
- 실패 시 `NULL` 반환 (경고)

## 함수 시그니처

```r
compute_meta_gene_importance(
  meta_result,      # TML7() 결과 객체
  normalize = TRUE, # 중요도를 최대값으로 정규화 여부
  target_model = NULL  # 사용할 모델 이름 (NULL이면 best_model)
)
```

## 반환값

- **성공 시**: Named list
  - 각 signature 이름을 키로 하고, gene-level importance를 값으로 가짐
  - 각 element는 named numeric vector (gene name → importance score)
- **실패 시**: `NULL` (경고 메시지 출력)

## 사용 예시

### 기본 사용
```r
# TML7 결과
tmla <- TML7(
  l1_signatures = fgsa,
  holdout_data = data_seurat,
  target_var = "g3",
  l2_methods = c("glm", "ranger", "xgbTree")
)

# Best model 사용
cmgi <- compute_meta_gene_importance(tmla)

# 특정 모델 지정
cmgi_glm <- compute_meta_gene_importance(tmla, target_model = "glm")
```

### List 처리
```r
# 여러 TML 객체 처리
tml_list <- list(
  "CD4+ T-cells" = tmla1,
  "CD8+ Cytotoxic T-cells" = tmla2,
  ...
)

# Best model 사용 (일부는 NULL 반환될 수 있음)
cmgi_list <- lapply(tml_list, function(x) {
  tryCatch(
    compute_meta_gene_importance(x),
    error = function(e) NULL,
    warning = function(w) NULL
  )
})

# NULL 제거
cmgi_list <- cmgi_list[!sapply(cmgi_list, is.null)]

# 특정 모델로 재시도 (best model이 importance 없는 경우)
failed_clusters <- names(tml_list)[!names(tml_list) %in% names(cmgi_list)]
cmgi_retry <- lapply(failed_clusters, function(cluster) {
  compute_meta_gene_importance(tml_list[[cluster]], target_model = "glm")
})
names(cmgi_retry) <- failed_clusters
```

### 결과 활용
```r
# 특정 signature의 gene importance
cmgi[["lasso"]]

# 특정 gene의 중요도 (모든 signature에 대해)
gene_importance <- sapply(cmgi, function(x) x[["CD3D"]])
gene_importance <- gene_importance[!is.na(gene_importance)]

# 시각화
library(ggplot2)
imp_df <- data.frame(
  signature = names(cmgi[["lasso"]]),
  importance = cmgi[["lasso"]]
)
imp_df <- imp_df[order(abs(imp_df$importance), decreasing = TRUE)[1:20], ]
ggplot(imp_df, aes(x = reorder(signature, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip()
```

## 모델별 Importance 지원 상태

| 모델 | Importance 제공 방식 | 기본 제공 | 비고 |
|------|---------------------|----------|------|
| **glm** | 회귀 계수 (`coef()`) | ✅ 항상 제공 | 가장 안정적 |
| **glmnet** | `caret::varImp()` | ✅ 항상 제공 | 안정적 |
| **ranger** | `ranger::importance()` 또는 `caret::varImp()` | ⚠️ 설정 필요 | 학습 시 `importance` 파라미터 필요 |
| **xgbTree** | `caret::varImp()` | ✅ 항상 제공 | 안정적 |
| **nnet** | `caret::varImp()` | ✅ 제공 | 신경망 가중치 기반 |
| **earth** | `caret::varImp()` | ❌ 제공 안 됨 | MARS 모델 특성상 제한적 |
| **svmRadial** | `caret::varImp()` | ⚠️ 제한적 | SVM 특성상 제한적 |

## 문제 해결

### Q: "Model was trained without importance calculation" 에러
**A:** 기존 TML 객체의 `ranger` 모델은 `importance` 파라미터가 설정되지 않았을 수 있습니다. `target_model`로 다른 모델을 사용하거나 재학습하세요:
```r
# 다른 모델 사용
cmgi <- compute_meta_gene_importance(tmla, target_model = "glm")

# 또는 재학습
tmla_new <- TML7(..., l2_methods = c("glm", "ranger", "xgbTree"))
```

### Q: "Some signatures not found in importance" 경고
**A:** 일부 signature가 zero variance로 제거되었거나, 모델이 해당 signature를 사용하지 않은 경우입니다. 정상적인 동작이며, 사용 가능한 signature만 반환됩니다.

### Q: `earth` 모델에서 `NULL` 반환
**A:** `earth` 모델은 `caret::varImp()`를 지원하지 않습니다. `target_model`로 다른 모델을 사용하거나, TML7에서 `earth`를 제외하고 학습하세요.

## 관련 문서

- [TML7 함수 문서](../README.md)
- [L2 Method Testing](./L2_METHOD_TESTING.md)
- [FGS/TML6 Analysis Context](./FGS_TML6_ANALYSIS_CONTEXT.md)

