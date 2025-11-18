# 질문에 대한 빠른 답변

## 1. GEM을 covar_effects에 포함하지 않는 것이 좋을까?

**네, 권장합니다!**

**이유:**
- GEM과 g3 간 완전 분리 문제 해결
- 설계 행렬이 특이하지 않음
- 분석이 안정적으로 수행됨
- 모델 단순화 및 해석 용이

**권장 설정:**
```r
neb1 <- runNEBULA(
  is5s, 
  fixed_effects = c("g3"), 
  covar_effects = NULL,  # GEM 제외
  patient_col = "hos_no", 
  offset = "nCount_RNA"
)
```

## 2. hos_no와 g3 간에도 완전 분리가 있는데 어떻게 되는 거지?

**문제 없습니다!**

**이유:**
- hos_no는 **random effect** (patient-level)
- g3는 **fixed effect**
- Random effect와 fixed effect 간 완전 분리는 문제가 되지 않음

**왜 문제가 되지 않는가?**
1. Random effects는 설계 행렬 X에 직접 포함되지 않음
2. Random effects는 평균 0인 정규분포를 따름 (b_i ~ N(0, σ²_b))
3. Patient 간 비교를 통해 fixed effects를 추정할 수 있음

**예시:**
```
Patient 1 (g3=1): log(μ_1j) = log(offset_1j) + β_g3=1 + b_1
Patient 2 (g3=2): log(μ_2j) = log(offset_2j) + β_g3=2 + b_2

Patient 간 비교를 통해 β_g3=2 - β_g3=1을 추정할 수 있음!
```

## 3. NEBULA의 가정, 모델 복잡성 등에 대해 상세히 설명

**상세한 설명은 `NEBULA_MODEL_EXPLANATION.md`를 참조하세요.**

**핵심 요약:**

### 모델 구조
```
Y_ij ~ Negative Binomial(μ_ij, φ)
log(μ_ij) = log(offset_ij) + X_ij^T β + b_i

where:
  X_ij: 설계 행렬 (fixed effects)
  β: fixed effects 계수
  b_i: random effect (patient-level)
      b_i ~ N(0, σ²_b)
```

### 주요 가정
1. **Negative Binomial 분포**: Count 데이터
2. **Log-linear 모델**: 선형 관계
3. **Random Effects**: 정규분포, patient 간 독립성
4. **Fixed Effects**: 설계 행렬이 full rank (특이하지 않음)

### 모델 복잡성
1. **Fixed Effects 수**: 더 많으면 더 복잡
2. **Random Effects 수**: 더 많은 patient → 더 많은 random effects
3. **유전자 수**: 계산 비용 선형 증가

### 완전 분리 문제
- **Fixed Effects 간 완전 분리**: 문제! (GEM과 g3)
- **Fixed Effects와 Random Effects 간 완전 분리**: 문제 없음! (hos_no와 g3)

## 요약

1. **GEM은 포함하지 않는 것을 권장**
2. **hos_no와 g3 간 완전 분리는 문제 없음**
3. **상세한 설명은 `NEBULA_MODEL_EXPLANATION.md` 참조**

## 참고 문서

- `NEBULA_MODEL_EXPLANATION.md`: NEBULA 모델 상세 설명
- `NEBULA_ERROR_EXPLANATION.md`: 오류 해결 방법
- `QUICK_ANSWER.md`: 빠른 답변 (본 문서)
