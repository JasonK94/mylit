# TML7 Gene Importance Improvement Documentation

## 개요

이 디렉토리는 TML7 (Train Meta-Learner v7)의 `compute_meta_gene_importance` 함수 개선을 위한 개발 계획과 구현을 관리합니다.

## 목적

현재 구현에서 여러 L2 모델(ranger, glm, xgbTree 등)의 계수를 단순 합산하여 gene-level importance를 계산하는 방식의 한계를 개선하기 위함입니다.

**핵심 제약사항**: AMSC(단순 합산)만으로도 성능이 좋기 때문에, 제안된 방법론이 명확한 성능 향상을 보여주지 않으면 사용하지 않습니다.

## 주요 문서

- **[DEVELOPMENT_PLAN.md](./DEVELOPMENT_PLAN.md)**: 전체 개발 계획 및 단계별 접근 방법
- 각 Phase별 상세 문서 (향후 추가 예정)

## 개발 단계 요약

### Phase 1: Model-specific Normalization
각 모델 타입에 맞는 정규화 방법 적용 (우선순위: 높음)

### Phase 2: Weighted Ensemble
여러 모델의 importance를 CV 성능으로 가중하여 통합 (우선순위: 중간)

### Phase 3: Gradient-based Importance
비선형 모델을 위한 SHAP/Gradient-based importance (우선순위: 낮음)

### Phase 4: Signature-specific Aggregation
Signature 타입에 따라 다른 집계 방법 적용 (우선순위: 낮음)

## 평가 기준

1. **성능 향상**: AMSC만 사용한 경우 대비 성능 개선
2. **해석 가능성**: Gene importance의 의미가 명확한가
3. **계산 비용**: 구현 및 실행 시간이 허용 가능한가
4. **Robustness**: 다양한 데이터셋에서 일관된 성능인가

## 현재 상태

- ✅ 개발 계획 수립 완료
- ⏳ Phase 1 구현 준비 중

## 참고사항

- 모든 구현은 기존 `compute_meta_gene_importance` 함수와 호환되어야 함
- `target_model` 파라미터를 통해 특정 모델 선택 가능해야 함
- 각 방법론은 독립적으로 사용 가능해야 함 (조합 가능)

