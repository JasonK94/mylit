# TML6 개선 작업 컨텍스트

**작업 브랜치**: `fgs`  
**작업 디렉토리**: `/home/user3/data_user3/git_repo/_wt/fgs`  
**작업 날짜**: 2025-01-XX  
**상태**: 진행 중 (환자 단위 CV 및 L2 메서드 확장 대기 중)

---

## 1. 완료된 작업

### 1.1 TML6 함수 버그 수정

#### 문제점
- `.score_signature()` 함수에서 `scale()` 사용 시 zero-variance signature가 있으면 NA가 발생
- 모든 signature가 zero-variance일 때 "No usable data remains after cleaning" 에러 발생
- xgboost의 deprecated 경고(`ntree_limit` → `iteration_range`)가 너무 많이 출력됨

#### 해결 방법
- `.score_signature()` 함수에서 `scale()` 전에 분산 체크 추가
  - 분산이 0이면 normalize를 건너뛰고 0으로 설정
  - `scale()` 결과가 NA/Inf인 경우 원래 점수 사용
- xgboost 경고 억제 시도:
  - `xgb.set.config(verbosity = 0)` 호출 (가능한 경우)
  - `suppressWarnings()`로 R-level 경고 억제
- 에러 메시지 개선: 더 자세한 디버깅 정보 제공

**파일**: `myR/R/signature.R` (라인 ~1149-1180)

### 1.2 compute_meta_gene_importance 함수 버그 수정

#### 문제점
- `signature_importance[[sig]]` 사용 시 "subscript out of bounds" 에러
- `caret::varImp`의 rownames와 `sig_names`가 일치하지 않을 수 있음

#### 해결 방법
- `signature_importance[[sig]]` → `signature_importance[sig]`로 변경 (named vector이므로 `[` 사용)
- `intersect()`로 실제 존재하는 signature만 추출
- 존재하지 않는 signature는 경고 후 건너뛰기
- NA/Inf importance 값 체크 및 처리

**파일**: `myR/R/signature.R` (라인 ~1521-1583)

### 1.3 add_meta_signature_score 함수 추가

#### 목적
TML6와 `compute_meta_gene_importance()`로 만든 gene weights를 사용하여 signature score를 계산하고 Seurat 객체에 AddMetaData로 추가

#### 기능
- Meta-learner의 gene-level importance를 가중치로 사용
- 가중합 계산: `sum(w * expr) / sum(|w|)`
- 선택적 z-score normalization
- Signed/absolute contribution 선택 가능
- Sparse/dense matrix 모두 지원

**파일**: `myR/R/signature.R` (라인 ~1616-1761)

### 1.4 문서화 및 커밋

- `DEVLOG_Korean.md`에 변경 사항 기록
- 커밋 해시: `edf439a`
- 커밋 메시지: "Fix TML6 zero-variance handling and compute_meta_gene_importance subscript error"

---

## 2. 진행 예정 작업

### 2.1 환자 단위 누수 방지용 group-wise CV 추가

#### 문제점
- 현재 TML6는 cell/AOI 단위로 랜덤 k-fold CV를 수행
- 같은 환자의 `pre`와 `post`가 서로 다른 fold로 나뉘어 train/test에 동시에 들어갈 수 있음
- 이로 인해 메타러너의 AUC가 과도하게 높게 나올 수 있음 (예: 0.9993)

#### 해결 방법
- **새 파라미터**: `cv_group_var = "emrid"` (기본값)
- Seurat 입력 + `meta.data`에 해당 컬럼이 있으면 group-wise CV 활성화
- `caret::trainControl`의 `index`/`indexOut`을 사용하여 같은 그룹은 같은 fold에 묶음
- `cv_group_var = NULL`로 설정하면 기존 방식(셀 단위 CV) 유지

#### 구현 위치
- `TML6` 함수 시그니처에 `cv_group_var = "emrid"` 파라미터 추가
- 타깃 준비/정리 후 `cv_group` 변수 추출 및 `keep`/`row_ok`와 동기화
- `trainControl` 생성 전에 group-wise CV index/indexOut 구성

**예상 코드 위치**: `myR/R/signature.R` (라인 ~1320-1379 근처)

### 2.2 L2 meta-learner 메서드 확장

#### 추가할 메서드
- `glmnet`: 정규화된 로지스틱 회귀 (L1/L2/elastic-net)
- `svmRadial`: RBF 커널 SVM (비선형 결정 경계)
- `mlp`: RSNNS 기반 다층 퍼셉트론
- `mlpKerasDropout`: Keras 기반 드롭아웃 MLP (고급 사용자용)

#### 구현 방법
- `TML6` 함수에서 패키지 가용성 체크 추가
- 패키지가 없으면 경고 후 해당 메서드 자동 제거
- caret이 이미 지원하므로 별도 구현 불필요

**예상 코드 위치**: `myR/R/signature.R` (라인 ~1381-1402 근처, xgbTree 체크 블록 앞)

---

## 3. 테스트 계획

### 3.1 단계별 테스트 순서

1. **기본 테스트**: glm + 셀 단위 CV
   ```r
   tml_glm <- TML6(fgs2, data_seurat, "response",
                   l2_methods = "glm", cv_group_var = NULL)
   ```

2. **환자 단위 CV 테스트**: glm + emrid 기반 group CV
   ```r
   tml_glm_grp <- TML6(fgs2, data_seurat, "response",
                       l2_methods = "glm", cv_group_var = "emrid")
   ```

3. **L2 메서드 확장 테스트**: 모든 메서드 포함
   ```r
   tml_deep <- TML6(fgs2, data_seurat, "response",
                    l2_methods = c("glm","glmnet","ranger","xgbTree",
                                   "svmRadial","mlp","mlpKerasDropout"),
                    cv_group_var = "emrid")
   ```

### 3.2 예상 결과

- 환자 단위 CV를 사용하면 AUC가 약간 낮아질 수 있음 (과적합 완화)
- 패키지가 없는 메서드는 자동으로 제거되고 경고 출력
- `model_comparison`으로 여러 메타러너의 CV 성능 비교 가능

---

## 4. 현재 코드 상태

### 4.1 TML6 함수 시그니처 (현재)
```r
TML6 <- function(
  l1_signatures,
  holdout_data,
  target_var,
  l2_methods = c("glm","ranger","xgbTree"),
  k_folds   = 5,
  metric    = c("AUC","ROC","Accuracy","Kappa"),
  fgs_seed  = 42,
  layer     = "data",
  allow_parallel = FALSE,
  parallel_workers = NULL
)
```

### 4.2 예상 변경 후 시그니처
```r
TML6 <- function(
  l1_signatures,
  holdout_data,
  target_var,
  l2_methods = c("glm","ranger","xgbTree"),
  k_folds   = 5,
  metric    = c("AUC","ROC","Accuracy","Kappa"),
  fgs_seed  = 42,
  layer     = "data",
  allow_parallel = FALSE,
  parallel_workers = NULL,
  cv_group_var = "emrid"   # ★ 신규 파라미터
)
```

---

## 5. 주요 파일 경로

- **메인 함수 파일**: `/home/user3/data_user3/git_repo/_wt/fgs/myR/R/signature.R`
- **문서**: `/home/user3/data_user3/git_repo/_wt/fgs/myR/docs/DEVLOG_Korean.md`
- **테스트 데이터**: `/data/user3/sobj/data_seurat_251104.qs` (참고용)
- **워크트리**: `/home/user3/data_user3/git_repo/_wt/fgs`

---

## 6. 기술적 세부사항

### 6.1 Group-wise CV 구현 방법

- `caret::trainControl`의 `index`와 `indexOut` 파라미터 사용
- 각 fold마다:
  - `index[[fold]]`: 해당 fold의 training set 인덱스 (같은 그룹의 모든 셀 포함)
  - `indexOut[[fold]]`: 해당 fold의 test set 인덱스 (다른 그룹의 모든 셀 포함)
- 그룹 수가 `k_folds`보다 작으면 자동으로 standard CV로 fallback

### 6.2 L2 메서드 패키지 의존성

- `glmnet`: `glmnet` 패키지
- `svmRadial`: `kernlab` 패키지
- `mlp`: `RSNNS` 패키지
- `mlpKerasDropout`: `keras` 패키지 (+ `tensorflow` 백엔드)

---

## 7. 다음 단계

1. ✅ **완료**: TML6 zero-variance 버그 수정
2. ✅ **완료**: compute_meta_gene_importance 버그 수정
3. ✅ **완료**: add_meta_signature_score 함수 추가
4. ⏳ **대기**: 환자 단위 CV (`cv_group_var`) 구현
5. ⏳ **대기**: L2 메서드 확장 (glmnet, svmRadial, mlp, mlpKerasDropout)
6. ⏳ **대기**: 테스트 및 검증
7. ⏳ **대기**: 최종 커밋 및 문서화

---

## 8. 참고사항

- **데이터 특성**: GeoMx dataset, 114~121개 AOI, `emrid` (환자 ID), `treatment` (pre/post), `response` (R/NR)
- **과적합 우려**: 현재 메타러너 AUC가 0.9993으로 매우 높음 → 환자 단위 CV로 재평가 필요
- **성능 vs 해석**: `pred` (TML6 예측값)는 성능이 좋지만, `meta_signature_score`는 해석 가능한 선형 요약

---

**마지막 업데이트**: 2025-01-XX  
**작성자**: Auto (GPT-5 Codex)

