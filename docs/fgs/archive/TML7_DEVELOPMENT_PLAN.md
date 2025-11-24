# TML7 Gene Importance 개선 개발 계획

## 배경

현재 `compute_meta_gene_importance` 함수는 여러 L2 모델(ranger, glm, xgbTree 등)에서 얻은 계수들을 단순 합산하여 gene-level importance를 계산합니다. 각 모델의 계수 의미가 다르고(로그 오즈, permutation importance, gain-based importance 등), 이들의 단순 합산은 정보 손실을 초래할 수 있습니다.

**중요 제약사항**: AMSC(단순 합산)만으로도 성능이 좋기 때문에, 너무 간단한 모델은 실제로 사용되지 않을 가능성이 높습니다.

## 개발 목표

1. Model-specific normalization으로 계수 의미 차이 해결
2. Weighted ensemble로 여러 모델의 정보 통합
3. 비선형 모델을 위한 Gradient/SHAP-based importance
4. 최종 성능 평가를 통한 방법론 선택

## 개발 단계

### Phase 1: Model-specific Normalization (우선순위: 높음)

**목적**: 각 모델 타입에 맞는 정규화 방법 적용

**구현 방법**:
- GLM: 계수는 이미 log-odds scale이므로 그대로 사용 또는 softmax 변환
- Ranger: Permutation importance를 probability scale로 변환
- XGBoost: Gain-based importance를 probability scale로 변환
- 기타: 일반적인 정규화 적용

**코드 위치**: `compute_meta_gene_importance` 함수 내부

**평가 지표**:
- Gene importance의 해석 가능성
- 기존 방법 대비 성능 변화
- 계산 비용

**예상 결과**: Model-specific normalization만으로도 계수 의미 차이 문제 해결 가능

---

### Phase 2: Weighted Ensemble (우선순위: 중간)

**목적**: 여러 모델의 importance를 CV 성능으로 가중하여 통합

**구현 방법**:
1. 각 L2 모델의 CV 성능 측정
2. 성능을 기반으로 가중치 계산 (softmax 또는 정규화)
3. 각 모델의 importance를 정규화 후 가중 평균
4. 가중 평균된 importance를 사용하여 gene contribution 계산

**장점**:
- Best model에만 의존하지 않음
- 여러 모델의 정보를 통합
- Robust한 importance 추정

**단점**:
- 계산 비용 증가
- 모든 모델의 importance 추출 필요

**평가 지표**:
- AMSC만 사용한 경우와 성능 비교
- 여러 모델 통합의 이점 평가
- 계산 시간 측정

**예상 결과**: Best model만 사용하는 것보다 더 robust할 수 있음. 단, 성능 향상이 명확해야 사용.

---

### Phase 3: Gradient-based Importance (우선순위: 낮음, 비선형 모델 전용)

**목적**: 비선형 모델(xgbTree, nnet, mlp 등)의 해석 가능성 향상

**구현 방법**:
- XGBoost: SHAP values 사용 (가능하면)
- Neural Network: Gradient-based importance (Integrated Gradients)
- Fallback: Permutation importance

**장점**:
- 비선형 모델의 feature importance를 더 정확하게 추정
- SHAP values는 해석 가능성이 높음

**단점**:
- 계산 비용이 매우 높음
- SHAP 패키지 의존성 추가
- 모든 모델에 적용 가능하지 않음

**평가 지표**:
- 기존 방법 대비 해석 가능성 향상
- 계산 시간
- SHAP vs Permutation importance 비교

**예상 결과**: 비선형 모델의 경우 유용할 수 있지만, 계산 비용이 높아 선택적으로만 사용.

---

### Phase 4: Signature-specific Aggregation (우선순위: 낮음)

**목적**: Signature 타입에 따라 다른 집계 방법 적용

**구현 방법**:
- 선형 모델 기반 signature: 합산 (현재 방식)
- Tree-based signature: 절댓값 합산 또는 최대값 선택
- 통계 모델 기반 signature: 가중 평균 (신뢰도 기반)

**장점**:
- Signature 타입의 특성을 반영
- 더 정확한 importance 추정

**단점**:
- 구현 복잡도 증가
- Signature 타입 분류 로직 필요

**평가 지표**:
- 기존 방법 대비 성능
- Signature 타입별 contribution의 의미

**예상 결과**: 약간의 개선은 있을 수 있지만, 구현 복잡도 대비 이득이 클지 않을 수 있음.

---

## 최종 평가 기준

각 방법론은 다음 기준으로 평가됩니다:

1. **성능 향상**: AMSC만 사용한 경우 대비 성능 개선 (AUC, Accuracy 등)
2. **해석 가능성**: Gene importance의 의미가 명확한가
3. **계산 비용**: 구현 및 실행 시간이 허용 가능한가
4. **Robustness**: 다양한 데이터셋에서 일관된 성능인가

**중요**: AMSC만으로도 성능이 좋기 때문에, 제안한 방법론이 **명확한 성능 향상**을 보여주지 않으면 사용하지 않습니다.

---

## 구현 순서

1. **Phase 1 (Model-specific Normalization)** 구현 및 평가
   - 가장 간단하고 효과적일 것으로 예상
   - 계산 비용 증가 최소

2. **Phase 2 (Weighted Ensemble)** 구현 및 평가
   - Phase 1 결과에 따라 필요성 판단
   - 성능 향상이 명확할 때만 사용

3. **Phase 3 (Gradient-based)** 구현 및 평가
   - 비선형 모델이 best model일 때만 고려
   - 계산 비용이 높으므로 선택적 사용

4. **Phase 4 (Signature-specific Aggregation)** 구현 및 평가
   - 다른 방법들이 효과적이지 않을 때만 고려
   - 구현 복잡도 대비 이득이 명확할 때만 사용

---

## 파일 구조

```
docs/tml/
├── DEVELOPMENT_PLAN.md          # 이 파일
├── MODEL_SPECIFIC_NORMALIZATION.md  # Phase 1 상세 계획
├── WEIGHTED_ENSEMBLE.md         # Phase 2 상세 계획
├── GRADIENT_BASED_IMPORTANCE.md # Phase 3 상세 계획
└── EVALUATION_RESULTS.md        # 평가 결과 기록
```

---

## 참고사항

- 모든 구현은 기존 `compute_meta_gene_importance` 함수와 호환되어야 함
- `target_model` 파라미터를 통해 특정 모델 선택 가능해야 함
- 각 방법론은 독립적으로 사용 가능해야 함 (조합 가능)
- 성능 평가는 실제 데이터셋에서 진행
- 문서화는 각 단계마다 업데이트

---

## 다음 단계

1. Phase 1 구현 시작
2. 평가 결과에 따라 다음 Phase 결정
3. 최종 방법론 선택 및 문서화
